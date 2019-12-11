#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "metadata.hpp"
#include "parser.hpp"
#include "background.hpp"
#include "chealpix.h"
#include <map>
#include <gsl/gsl_spline.h>

// g++ raytracer.cpp -o raytracer -std=c++11 -O3 -fopenmp -DCOLORTERMINAL -DLIGHTCONE_THICKNESS=3 -I/astro/adamek/local/include -L/astro/adamek/local/lib -lgsl -lgslcblas -lchealpix -lcfitsio

#define HEALPIX_UNDEF -1.6375e30
#ifndef READAHEAD
#define READAHEAD 16
#endif

using namespace std;

struct healpix_header
{
	uint32_t Nside;
	uint32_t Npix;
	uint32_t precision;
	uint32_t Ngrid;
	double direction[3];
	double distance;
	double boxsize;
	char fill[256 - 4 * 4 - 5 * 8]; /* fills to 256 Bytes */
	uint32_t header_blocksize;
	uint32_t data_blocksize; /* to get rid of two fread commands per header/data pair */
};

struct cycle_info
{
	int num;
	double tau;
	double a;
	double d_in;
	double d_01;
	double d_12;
	double d_out;
	double opening[3];
};

struct photon_data
{
	double pos[3];
	double n[3];
	double e1[3];
	double e2[3];
	double lnk0;
	double delay;
	double DA;
	double DAprime;
	double sigma[2];
	double ellip[2];
	double omega;
	int32_t stop;
	int32_t padding;
};

struct metric_data
{
	healpix_header hdr;
	float * pixel;
};

struct metric_container
{
	std::map<int,metric_data> healpix_data;
	cycle_info cinfo;
	char name[4];
	char * dir;
	char * basename;
	void init(cycle_info & c, char * d, char * b, const char * n)
	{
		cinfo = c;
		strcpy(name, n);
		dir = d;
		basename = b;
	}
	void truncate(double min_dist)
	{
		std::map<int,metric_data>::iterator it;

		for (it = healpix_data.begin(); it != healpix_data.end(); it++)
		{
			if (it->second.hdr.distance >= min_dist) break;
			else free(it->second.pixel);
		}

		if (it != healpix_data.begin()) healpix_data.erase(healpix_data.begin(), it);
	}
	void clear()
	{
		for (std::map<int,metric_data>::iterator it = healpix_data.begin(); it != healpix_data.end(); it++)
			free(it->second.pixel);
		healpix_data.clear();
	}
};

int loadHealpixData(metric_container * field, double min_dist, double max_dist);

inline double linear_interpolation(double field1, double field2, double weight)
{
	return field1 * (1.-weight) + field2 * weight;
}

void rotate_vector(double * input, double * output, double * coordsys);

void counterrotate_vector(double * input, double * output, double * coordsys);

double get_direction(double * input, double * output, double * coordsys)
{
	double dist = sqrt(input[0]*input[0] + input[1]*input[1] + input[2]*input[2]);
	if (dist > 0)
	{
		rotate_vector(input, output, coordsys);
		output[0] /= dist;
		output[1] /= dist;
		output[2] /= dist;
		if (output[0] == 0 && output[1] == 0) output[3] = 0;
		else output[3] = atan2(output[1], output[0]);
		if (output[3] < 0) output[3] += 2 * M_PI;
		return dist;
	}
	else
	{
		output[0] = 0;
		output[1] = 0;
		output[2] = 1;
		output[3] = 0;
		return 0;
	}
}

std::map<int,metric_data>::iterator get_iterator(metric_container * field, const double dist, const int dir);

bool interpolation(metric_container * field_prev, metric_container * field_next, const double weight, const double dist, const double dr, double * v, double * result);

bool time_derivative(metric_container * field_prev, metric_container * field_next, const double dtau, const double dist, const double dr, double * v, double * result);

bool gradient(metric_container * field_prev, metric_container * field_next, const double weight, const double dist, const double dr, double * v, double * result, int debug = 0);

bool time_derivative_gradient(metric_container * field_prev, metric_container * field_next, const double dtau, const double dist, const double dr, double * v, double * result, int debug = 0);

bool tidal_matrix(metric_container * field_prev, metric_container * field_next, const double weight, const double dist, const double dr, double * v, double * result, int debug = 0);

int main(int argc, char **argv)
{
	int numtasks = 1;
	int taskindex = 0;
	char * settingsfile = NULL;

	parameter * params = NULL;
	int numparam;

	long Nsource = 0;
	long Nsource_tot = 0;
	long Nsource_skip = 0;
	int num_passes = 2;
	double integration_step = 1.0;
	char parameter_filename[PARAM_MAX_LENGTH];
	char output_dir[PARAM_MAX_LENGTH];
	char dummy[PARAM_MAX_LENGTH];
	char catalog_filename[PARAM_MAX_LENGTH];
	char output_filename[PARAM_MAX_LENGTH];
	char restart_filename[PARAM_MAX_LENGTH];

	metadata sim;
	icsettings ic;
	cosmology cosmo;

	bool use_lightcone[MAX_OUTPUTS];
	char filename[1024];
	int num_usedlightcone = 0;
	int lightcone_fields = MASK_PHI;
	std::map<int,cycle_info> lightcone_info;
	std::map<int,cycle_info>::iterator shell;
	cycle_info cinfo;
	double tau_obs = -1.;
	double z_obs;
	double dummy_read[LIGHTCONE_THICKNESS+5];
	double * vertex = NULL;
	double fourpiG;

	FILE * infile;
	FILE * outfile;

	size_t num_read1, num_read2;

	int index_x=-1, index_y=-1, index_z=-1, index_vx=-1, index_vy=-1, index_vz=-1, index = 0;
	double sourcedata[6];
	char line[1024];
	char format[1024];
	char * token;
	char c;
	long backtrack;
	int passes_done = 0;

	double * source_pos;
	double * source_vel;
	double * source_direction;
	double meancorr, maxcorr;
	int count;
	double * tauvec;
	double * avec;
	gsl_interp_accel * acc;
	gsl_spline * aspline;

	double kobs = 0.;
	double tau;
	double weight;
	double dr = -1.;
	double coordsys[4];
	photon_data * photon;
	long num_live;
	long stepcount;

	metric_container phi_next;
	metric_container phi_prev;

	metric_container B1_next;
	metric_container B1_prev;
	metric_container B2_next;
	metric_container B2_prev;
	metric_container B3_next;
	metric_container B3_prev;

	cycle_info cycle_next, cycle_prev;
	double min_dist = 0;
	double max_dist = 0;
	double buffer_dist;
	double dist, source_dist;

	double phi;
	double B[3] = {0, 0, 0};
	double gradB1[3] = {0, 0, 0};
	double gradB2[3] = {0, 0, 0};
	double gradB3[3] = {0, 0, 0};
	double gradnB[3] = {0, 0, 0};
	double grade1B[3] = {0, 0, 0};
	double grade2B[3] = {0, 0, 0};
	double gradBprime[3] = {0, 0, 0};
	double tidalB[6];
	double nngradB, tidalBterm;
	double gradphi[3];
	double tidalphi[6];
	double phiprime;
	double ngradphi, nntidalphi;
	double e1gradphi, e2gradphi;
	double eetidalphi[2];
	double eetidalBterm[2];
	double dn_dtau[3];
	double dlnk0_dtau;
	double exp2phi;
	double v[4];
	double tmp[3];

	for (int i = 1; i < argc; i++ ) {
		if ( argv[i][0] != '-' )
			continue;
		switch(argv[i][1]) {
			case 'i':
				taskindex = atoi(argv[++i]);
				break;
			case 't':
				numtasks =  atoi(argv[++i]);
				break;
			case 's':
				settingsfile = argv[++i];
			}
	}

	if (settingsfile == NULL)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no settings file specified!" << endl;
		return -1;
	}

	numparam = loadParameterFile(settingsfile, params);

	if (!parseParameter(params, numparam, "gevolution parameter file", parameter_filename))
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": gevolution parameter file not specified!" << endl;
		return -1;
	}

	if (!parseParameter(params, numparam, "gevolution output directory", output_dir))
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": gevolution output directory not specified!" << endl;
		return -1;
	}

	parseParameter(params, numparam, "integration time step", integration_step);

	if (!parseParameter(params, numparam, "output file", output_filename))
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": output file not specified!" << endl;
		return -1;
	}

	if (!parseParameter(params, numparam, "catalog file", catalog_filename))
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": catalog file not specified!" << endl;
		return -1;
	}

	parseParameter(params, numparam, "passes", num_passes);

	parseParameter(params, numparam, "metric dr", dr);

	if (!parseParameter(params, numparam, "restart file", restart_filename))
	{
		cout << COLORTEXT_YELLOW << " warning" << COLORTEXT_RESET << ": restart file not specified! using default: raytracer_restart.bin" << endl;
		sprintf(restart_filename, "raytracer_restart.bin");
	}

	free(params);

	if (taskindex == 0)
		cout << COLORTEXT_WHITE << " LCARS tools: raytracer" << endl << COLORTEXT_RESET << " compiled with " << COLORTEXT_RED << "LIGHTCONE_THICKNESS=" << LIGHTCONE_THICKNESS << endl << COLORTEXT_RESET << endl << " opening settings file of simulation: " << parameter_filename << endl << " parser output:" << endl << endl;

	numparam = loadParameterFile(parameter_filename, params);
	parseMetadata(params, numparam, sim, cosmo, ic);

	free(params);

	if (dr < 0) dr = 1. / (double) sim.numpts;
	integration_step /= -sim.numpts;

	fourpiG = 1.5 * sim.boxsize * sim.boxsize / C_SPEED_OF_LIGHT / C_SPEED_OF_LIGHT;

	if (taskindex == 0)
	{
		cout << " number of light cones: " << sim.num_lightcone << endl << endl;

		for (int i = 0; i < sim.num_lightcone; i++)
		{
			cout << " light cone " << i << " parameters:" << endl << "    vertex = " << sim.lightcone[i].vertex[0] << ", " << sim.lightcone[i].vertex[1] << ", " << sim.lightcone[i].vertex[2] << endl;
			cout << "    redshift of observation = " << sim.lightcone[i].z << endl;
			cout << "    direction = " << sim.lightcone[i].direction[0] << ", " << sim.lightcone[i].direction[1] << ", " << sim.lightcone[i].direction[2] << endl;
			cout << "    opening half-angle = " << ((sim.lightcone[i].opening > -1.) ? acos(sim.lightcone[i].opening) * 180. / M_PI : 180.) /*sim.lightcone[i].opening * 180. / M_PI*/ << " degrees" << endl;
			cout << "    distance interval = " << sim.lightcone[i].distance[0] << ", " << sim.lightcone[i].distance[1] << endl << endl;
		}
	}

	for (int i = 0; i < sim.num_lightcone; i++)
	{
		if (sim.out_lightcone[i] & MASK_PHI)
  		{
			use_lightcone[i] = true;
			if (sim.out_lightcone[i] & MASK_CHI)
				lightcone_fields |= MASK_CHI;
			if (sim.out_lightcone[i] & MASK_B)
				lightcone_fields |= MASK_B;
		}
		else use_lightcone[i] = false;
	}

	if (taskindex == 0)
		cout << " using light cone(s) ";
	for (int i = 0; i < sim.num_lightcone; i++)
	{
		if (use_lightcone[i])
		{
			if (taskindex == 0)
				cout << i << " ";
			num_usedlightcone++;
			if (sim.lightcone[i].direction[0] == 0 && sim.lightcone[i].direction[1] == 0)
			{
				coordsys[0] = 1;
				coordsys[1] = 0;
				coordsys[3] = 0;
			}
			else
			{
				coordsys[2] = atan2(sim.lightcone[i].direction[1], sim.lightcone[i].direction[0]);
				coordsys[0] = cos(coordsys[2]);
				coordsys[1] = sin(coordsys[2]);
				coordsys[3] = sqrt(1. - sim.lightcone[i].direction[2]*sim.lightcone[i].direction[2]);
			}
			coordsys[2] = sim.lightcone[i].direction[2];
		}
	}
	if (taskindex == 0)
		cout << endl << endl;

	if (num_usedlightcone < 1)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no usable light cones!" << endl;
		return -1;
	}

	cout << " reading light cone information files...: " << endl;
	for (int i = 0; i < sim.num_lightcone; i++)
	{
		if (use_lightcone[i])
		{
			if (tau_obs < 0.)
			{
				vertex = sim.lightcone[i].vertex;
				z_obs = sim.lightcone[i].z;
				tau_obs = particleHorizon(1. / (1. + z_obs), fourpiG, cosmo);
			}
			else if (sim.lightcone[i].z != z_obs || sim.lightcone[i].vertex[0] != vertex[0] || sim.lightcone[i].vertex[1] != vertex[1] || sim.lightcone[i].vertex[2] != vertex[2])
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": not all selected light cones refer to the same observation event!" << endl << " incompatible lightcone " << i << " will be omitted!" << endl << endl;
				use_lightcone[i] = false;
				continue;
			}

			if (sim.num_lightcone > 1)
				sprintf(filename, "%s%s%d_info.bin", output_dir, sim.basename_lightcone, i);
			else
				sprintf(filename, "%s%s_info.bin", output_dir, sim.basename_lightcone);

			infile = fopen(filename, "rb");

			if (infile == NULL)
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to open light cone information file " << filename << endl;
				return -1;
			}

			while (true)
			{
				num_read1 = fread((void *) &cinfo.num, sizeof(int), 1, infile);
				num_read2 = fread((void *) dummy_read, sizeof(double), LIGHTCONE_THICKNESS+5, infile);

				if (feof(infile) || ferror(infile)) break;
				else if (num_read1 != 1 || num_read2 != LIGHTCONE_THICKNESS+5)
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": read error in light cone information file " << filename << " (first pass)" << endl;
					fclose(infile);
					return -1;
				}

				cinfo.tau = dummy_read[0];
				cinfo.a = dummy_read[1];
				cinfo.d_in = dummy_read[4];
				cinfo.d_01 = dummy_read[5];
				cinfo.d_12 = dummy_read[6];
				cinfo.d_out = dummy_read[7];
				cinfo.opening[0] = (sim.lightcone[i].distance[0] <= cinfo.d_01 || sim.lightcone[i].distance[1] > cinfo.d_01 || cinfo.d_01 <= 0) ? 2. : sim.lightcone[i].opening;
				cinfo.opening[1] = (sim.lightcone[i].distance[0] <= cinfo.d_12 || sim.lightcone[i].distance[1] > cinfo.d_12 || cinfo.d_12 <= 0) ? 2. : sim.lightcone[i].opening;
				cinfo.opening[2] = (sim.lightcone[i].distance[0] <= cinfo.d_out || sim.lightcone[i].distance[1] > cinfo.d_out || cinfo.d_out <= 0) ? 2. : sim.lightcone[i].opening;

				if ((shell = lightcone_info.find(cinfo.num)) != lightcone_info.end())
				{
					for (int l = 0; l < 3; l++)
					{
						if (cinfo.opening[l] < shell->second.opening[l])
						{
							shell->second.opening[l] = cinfo.opening[l];
						}
					}
				}
				else
					lightcone_info.insert(std::pair<int,cycle_info>(cinfo.num, cinfo));
			}

			fclose(infile);
		}
	}

	tauvec = (double *) malloc(lightcone_info.size() * sizeof(double));
	avec = (double *) malloc(lightcone_info.size() * sizeof(double));

	//cout << " scale factor table:" << endl;

	shell = lightcone_info.begin();
	for (int i = 0; shell != lightcone_info.end(); shell++, i++)
	{
		tauvec[i] = shell->second.tau;
		avec[i] = shell->second.a;

		//cout << " " << shell->first << ": " << tauvec[i] << ", " << avec[i] << endl;
	}

	acc = gsl_interp_accel_alloc();
	aspline = gsl_spline_alloc(gsl_interp_cspline, lightcone_info.size());
	gsl_spline_init(aspline, tauvec, avec, lightcone_info.size());

	free(tauvec);
	free(avec);

	if (taskindex == 0)
		cout << endl << " reading source catalog from " << catalog_filename << " ..." << endl;

	infile = fopen(catalog_filename, "r");

	if(infile == NULL || fgets(line, 1024, infile) == NULL)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": could not read from catalog file!" << endl;
		return -1;
	}

	if ((token = strchr(line, '\n')) != NULL) (*token) = '\0';

	token = strtok(line, " \t;,");

	for (int marker = 0; token != NULL; marker += 4)
	{
		if (strcmp(token, "x") == 0 || strcmp(token, "X") == 0)
		{
			sprintf(format+marker, " %%lf");
			index_x = index++;
		}
		else if (strcmp(token, "y") == 0 || strcmp(token, "Y") == 0)
		{
			sprintf(format+marker, " %%lf");
			index_y = index++;
		}
		else if (strcmp(token, "z") == 0 || strcmp(token, "Z") == 0)
		{
			sprintf(format+marker, " %%lf");
			index_z = index++;
		}
		else if (strcmp(token, "vx") == 0 || strcmp(token, "VX") == 0)
		{
			sprintf(format+marker, " %%lf");
			index_vx = index++;
		}
		else if (strcmp(token, "vy") == 0 || strcmp(token, "VY") == 0)
		{
			sprintf(format+marker, " %%lf");
			index_vy = index++;
		}
		else if (strcmp(token, "vz") == 0 || strcmp(token, "VZ") == 0)
		{
			sprintf(format+marker, " %%lf");
			index_vz = index++;
		}
		else sprintf(format+marker, " %%*f");

		token = strtok(NULL, " \t;,");
	}

	if (index_x < 0 || index_y < 0 || index_z < 0 || index_vx < 0 || index_vy < 0 || index_vz < 0)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": first line of catalog file does not contain expected header information!" << endl;
		return -1;
	}

#ifdef DEBUG
	if (taskindex == 0)
		cout << " DEBUG: format string for catalog reads" << endl << format << endl;
#endif

	backtrack = ftell(infile);

	while (!feof(infile) && !ferror(infile))
	{
		if (fscanf(infile, "%c%*[^\n]\n", &c) != 1) break;

		if (c == '#')
			backtrack = ftell(infile);
		else
			Nsource_tot++;
	}

	if (Nsource_tot % numtasks == 0)
	{
		Nsource = Nsource_tot/numtasks;
		Nsource_skip = taskindex * Nsource;
	}
	else if (taskindex == numtasks-1)
	{
		Nsource = Nsource_tot - (numtasks-1) * (1 + Nsource_tot/numtasks);
		Nsource_skip = (numtasks-1) * (1 + Nsource_tot/numtasks);
	}
	else
	{
		Nsource = 1 + Nsource_tot/numtasks;
		Nsource_skip = taskindex * (1 + Nsource_tot/numtasks);
	}

	source_pos = (double *) malloc(3l*Nsource*sizeof(double));
	source_vel = (double *) malloc(3l*Nsource*sizeof(double));

	fseek(infile, backtrack, SEEK_SET);

	for (int n = 0; n < Nsource_skip; n++)
	{
		if((c = fscanf(infile, format, sourcedata, sourcedata+1, sourcedata+2, sourcedata+3, sourcedata+4, sourcedata+5)) != 6)
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": could not read source data from catalog for source number " << n << ", number of items read = " << c << endl;
			return -1;
		}
	}

	for (int n = 0; n < Nsource; n++)
	{
		if((c = fscanf(infile, format, sourcedata, sourcedata+1, sourcedata+2, sourcedata+3, sourcedata+4, sourcedata+5)) != 6)
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": could not read source data from catalog for source number " << n << ", number of items read = " << c << endl;
			return -1;
		}
		source_pos[3*n] = sourcedata[index_x] / sim.boxsize; // assuming Mpc/h
		source_pos[3*n+1] = sourcedata[index_y] / sim.boxsize;
		source_pos[3*n+2] = sourcedata[index_z] / sim.boxsize;
		source_vel[3*n] = sourcedata[index_vx] / (100.0 * C_SPEED_OF_LIGHT); // assuming km/s
		source_vel[3*n+1] = sourcedata[index_vy] / (100.0 * C_SPEED_OF_LIGHT);
		source_vel[3*n+2] = sourcedata[index_vz] / (100.0 * C_SPEED_OF_LIGHT);
	}

	fclose(infile);

	cout << " task " << taskindex << " done! " << Nsource << " of " << Nsource_tot << " sources read from catalog." << endl << endl;

	source_direction = (double *) malloc(2l*Nsource*sizeof(double));

	photon = (photon_data *) malloc(Nsource*sizeof(photon_data));

	for (int n = 0; n < Nsource; n++)
	{
		photon[n].stop = 0;
		photon[n].delay = 0;
	}

	if (numtasks > 1)
	{
		sprintf(filename, "%s.%d", output_filename, taskindex);
		infile = fopen(filename, "r");	
	}
	else
		infile = fopen(output_filename, "r");

	if (infile == NULL)
	{
#pragma omp parallel for
		for (int n = 0; n < Nsource; n++)
		{
			source_direction[2*n] = acos(source_pos[3*n+2] / sqrt(source_pos[3*n]*source_pos[3*n] + source_pos[3*n+1]*source_pos[3*n+1] + source_pos[3*n+2]*source_pos[3*n+2]));
			source_direction[2*n+1] = atan2(source_pos[3*n+1], source_pos[3*n]);
		}
	}
	else
	{
		cout << " task " << taskindex << " reading source angles from previous pass..." << endl << endl;
		
		if (fgets(line, 1024, infile) == NULL || fscanf(infile, "# output for pass %d\n", &passes_done) != 1)
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": could not read header from previous output!" << endl;
			fclose(infile);
			return -1;
		}
		fgets(line, 1024, infile);

		for (int n = 0; n < Nsource; n++)
		{
			if((c = fscanf(infile, " %lf %lf %lf %*f %*e %*e %*e %*f %*f %*f %*f %*f", source_direction+(2*n), source_direction+(1+2*n), dummy_read)) != 3)
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": could not read source data from catalog for source number " << n << ", number of items read = " << c << endl;
				return -1;
			}

			if (dummy_read[0] < 0.) photon[n].stop = -1;
		}

		fclose(infile);
	}

	if (numtasks > 1)
	{
		sprintf(filename, "%s.%d", restart_filename, taskindex);
		infile = fopen(filename, "rb");	
	}
	else
		infile = fopen(restart_filename, "rb");

	if (infile != NULL)
	{
		cout << " task " << taskindex << " restarting from restart point..." << endl;
		if (fread((void *) &tau, sizeof(double), 1, infile) != 1 || fread((void *) &stepcount, sizeof(long), 1, infile) != 1 || fread((void *) &buffer_dist, sizeof(double), 1, infile) != 1 || fread((void *) &min_dist, sizeof(double), 1, infile) != 1 || fread((void *) &max_dist, sizeof(double), 1, infile) != 1)
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": could not read header from restart file!" << endl;
			fclose(infile);
			return -1;
		}
		if (fread((void *) photon, sizeof(photon_data), Nsource, infile) != Nsource)
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": could not read photon data from restart file!" << endl;
			fclose(infile);
			return -1;
		}
		if (fread((void *) source_pos, sizeof(double), 3l*Nsource, infile) != 3l*Nsource)
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": could not read source data from restart file!" << endl;
			fclose(infile);
			return -1;
		}
		fclose(infile);

		cout << " restart at step " << stepcount << ", tau = " << tau << ", (min,max,buffer)_dist = " << min_dist << ", " << max_dist << ", " << buffer_dist << endl << endl;
	}
	else
	{
		tau = tau_obs;
		stepcount = 1;
		buffer_dist = READAHEAD * dr - integration_step;
		min_dist = 0;
		max_dist = 0;
	}

	for (int pass = passes_done+1; pass < num_passes+1; pass++)
	{
		for (shell = lightcone_info.begin(); shell != lightcone_info.end() && shell->second.tau < tau; shell++);
		if (shell == lightcone_info.end())
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": tau out of range covered in light cone information files!" << endl;
			return -1;
		}

		cycle_next = shell->second;

		shell--;

		cycle_prev = shell->second;

		weight = (tau - cycle_prev.tau) / (cycle_next.tau - cycle_prev.tau);

		phi_prev.init(cycle_prev, output_dir, sim.basename_lightcone, "phi");
		phi_next.init(cycle_next, output_dir, sim.basename_lightcone, "phi");

		if (loadHealpixData(&phi_prev, min_dist - 1.7321 * dr, buffer_dist) < 0)
		{
			cout << COLORTEXT_YELLOW << " warning" << COLORTEXT_RESET << ": read-ahead failed on phi_prev at initialisation!" << endl;
			buffer_dist = max_dist;
		}
		if (loadHealpixData(&phi_next, min_dist - 1.7321 * dr, buffer_dist) < 0)
		{
			cout << COLORTEXT_YELLOW << " warning" << COLORTEXT_RESET << ": read-ahead failed on phi_next at initialisation!" << endl;
			buffer_dist = max_dist;
		}
		
		if (lightcone_fields & MASK_B)
		{
			B1_prev.init(cycle_prev, output_dir, sim.basename_lightcone, "B1");
			B1_next.init(cycle_next, output_dir, sim.basename_lightcone, "B1");

			B2_prev.init(cycle_prev, output_dir, sim.basename_lightcone, "B2");
			B2_next.init(cycle_next, output_dir, sim.basename_lightcone, "B2");

			B3_prev.init(cycle_prev, output_dir, sim.basename_lightcone, "B3");
			B3_next.init(cycle_next, output_dir, sim.basename_lightcone, "B3");

			if (loadHealpixData(&B1_prev, min_dist - 1.7321 * dr, buffer_dist) < 0)
			{
				cout << COLORTEXT_YELLOW << " warning" << COLORTEXT_RESET << ": read-ahead failed on B1_prev at initialisation!" << endl;
				buffer_dist = max_dist;
			}
			if (loadHealpixData(&B1_next, min_dist - 1.7321 * dr, buffer_dist) < 0)
			{
				cout << COLORTEXT_YELLOW << " warning" << COLORTEXT_RESET << ": read-ahead failed on B1_next at initialisation!" << endl;
				buffer_dist = max_dist;
			}

			if (loadHealpixData(&B2_prev, min_dist - 1.7321 * dr, buffer_dist) < 0)
			{
				cout << COLORTEXT_YELLOW << " warning" << COLORTEXT_RESET << ": read-ahead failed on B2_prev at initialisation!" << endl;
				buffer_dist = max_dist;
			}
			if (loadHealpixData(&B2_next, min_dist - 1.7321 * dr, buffer_dist) < 0)
			{
				cout << COLORTEXT_YELLOW << " warning" << COLORTEXT_RESET << ": read-ahead failed on B2_next at initialisation!" << endl;
				buffer_dist = max_dist;
			}

			if (loadHealpixData(&B3_prev, min_dist - 1.7321 * dr, buffer_dist) < 0)
			{
				cout << COLORTEXT_YELLOW << " warning" << COLORTEXT_RESET << ": read-ahead failed on B3_prev at initialisation!" << endl;
				buffer_dist = max_dist;
			}
			if (loadHealpixData(&B3_next, min_dist - 1.7321 * dr, buffer_dist) < 0)
			{
				cout << COLORTEXT_YELLOW << " warning" << COLORTEXT_RESET << ": read-ahead failed on B3_next at initialisation!" << endl;
				buffer_dist = max_dist;
			}
		}

		if (tau == tau_obs)
		{
			cout << " task " << taskindex << " setting photon initial conditions..." << endl;

			if (!interpolation(&phi_prev, &phi_next, weight, 0, dr, v, &phi))
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to load metric data at light cone vertex!" << endl;
				return -1;
			}

			exp2phi = exp(2 * phi);

			kobs = exp(-3. * phi);

			kobs *= gsl_spline_eval (aspline, tau_obs, acc);

#pragma omp parallel for
			for (int n = 0; n < Nsource; n++)
			{
				if (photon[n].stop < 0) continue;
				photon[n].pos[0] = 0;
				photon[n].pos[1] = 0;
				photon[n].pos[2] = 0;
				photon[n].n[0] = -sin(source_direction[2*n]) * cos(source_direction[2*n+1]);
				photon[n].n[1] = -sin(source_direction[2*n]) * sin(source_direction[2*n+1]);
				photon[n].n[2] = -cos(source_direction[2*n]);
				photon[n].e2[0] = sin(source_direction[2*n+1]);
				photon[n].e2[1] = -cos(source_direction[2*n+1]);
				photon[n].e2[2] = 0;
				photon[n].e1[0] = cos(source_direction[2*n]) * cos(source_direction[2*n+1]);
				photon[n].e1[1] = cos(source_direction[2*n]) * sin(source_direction[2*n+1]);
				photon[n].e1[2] = -sin(source_direction[2*n]);
				/*photon[n].e1[0] = sin(source_direction[2*n+1]);
				photon[n].e1[1] = -cos(source_direction[2*n+1]);
				photon[n].e1[2] = 0;
				photon[n].e2[0] = -cos(source_direction[2*n]) * cos(source_direction[2*n+1]);
				photon[n].e2[1] = -cos(source_direction[2*n]) * sin(source_direction[2*n+1]);
				photon[n].e2[2] = sin(source_direction[2*n]);*/
				photon[n].lnk0 = log(kobs);
				photon[n].DA = 0;
				photon[n].DAprime = -exp2phi;
				photon[n].sigma[0] = 0;
				photon[n].sigma[1] = 0;
				photon[n].ellip[0] = 0;
				photon[n].ellip[1] = 0;
				photon[n].omega = 0;
				photon[n].stop = 0;
			}
		}

		cout << " task " << taskindex << " starting time integration on pass " << pass << endl << endl;
		
		do
		{
			num_live = 0;

			if (shell->second.tau > tau)
			{
				while (shell->second.tau > tau && shell != lightcone_info.begin())
					shell--;

				if (shell->second.tau > tau) break; // leaving time domain

				cycle_next = cycle_prev;
				cycle_prev = shell->second;

				phi_next.healpix_data.swap(phi_prev.healpix_data);
				phi_next.cinfo = cycle_next;

				phi_prev.clear();
				phi_prev.cinfo = cycle_prev;

				if (loadHealpixData(&phi_prev, min_dist - 1.7321 * dr, buffer_dist) < 0)
				{
					cout << COLORTEXT_YELLOW << " warning" << COLORTEXT_RESET << ": read-ahead failed on phi_prev at cycle transition! step " << stepcount << ", available distance range = ";
					if (phi_prev.healpix_data.empty())
						cout << "none";
					else
						cout << phi_prev.healpix_data.begin()->second.hdr.distance << "-" << phi_prev.healpix_data.rbegin()->second.hdr.distance;
					cout << ", required distance range = " << min_dist-dr << "-" << max_dist+dr << endl;
				}

				if (lightcone_fields & MASK_B)
				{
					B1_next.healpix_data.swap(B1_prev.healpix_data);
					B1_next.cinfo = cycle_next;

					B1_prev.clear();
					B1_prev.cinfo = cycle_prev;

					if (loadHealpixData(&B1_prev, min_dist - 1.7321 * dr, buffer_dist) < 0)
					{
						cout << COLORTEXT_YELLOW << " warning" << COLORTEXT_RESET << ": read-ahead failed on B1_prev at cycle transition! step " << stepcount << ", available distance range = ";
						if (B1_prev.healpix_data.empty())
							cout << "none";
						else
							cout << B1_prev.healpix_data.begin()->second.hdr.distance << "-" << B1_prev.healpix_data.rbegin()->second.hdr.distance;
						cout << ", required distance range = " << min_dist-dr << "-" << max_dist+dr << endl;
					}

					B2_next.healpix_data.swap(B2_prev.healpix_data);
					B2_next.cinfo = cycle_next;

					B2_prev.clear();
					B2_prev.cinfo = cycle_prev;

					if (loadHealpixData(&B2_prev, min_dist - 1.7321 * dr, buffer_dist) < 0)
					{
						cout << COLORTEXT_YELLOW << " warning" << COLORTEXT_RESET << ": read-ahead failed on B2_prev at cycle transition! step " << stepcount << ", available distance range = ";
						if (B2_prev.healpix_data.empty())
							cout << "none";
						else
							cout << B2_prev.healpix_data.begin()->second.hdr.distance << "-" << B2_prev.healpix_data.rbegin()->second.hdr.distance;
						cout << ", required distance range = " << min_dist-dr << "-" << max_dist+dr << endl;
					}

					B3_next.healpix_data.swap(B3_prev.healpix_data);
					B3_next.cinfo = cycle_next;

					B3_prev.clear();
					B3_prev.cinfo = cycle_prev;

					if (loadHealpixData(&B3_prev, min_dist - 1.7321 * dr, buffer_dist) < 0)
					{
						cout << COLORTEXT_YELLOW << " warning" << COLORTEXT_RESET << ": read-ahead failed on B3_prev at cycle transition! step " << stepcount << ", available distance range = ";
						if (B3_prev.healpix_data.empty())
							cout << "none";
						else
							cout << B3_prev.healpix_data.begin()->second.hdr.distance << "-" << B3_prev.healpix_data.rbegin()->second.hdr.distance;
						cout << ", required distance range = " << min_dist-dr << "-" << max_dist+dr << endl;
					}
				}
			}

			weight = (tau - cycle_prev.tau) / (cycle_next.tau - cycle_prev.tau);

			phi_prev.truncate(min_dist - 2.7321 * dr);
			phi_next.truncate(min_dist - 2.7321 * dr);

			if (lightcone_fields & MASK_B)
			{
				B1_prev.truncate(min_dist - 2.7321 * dr);
				B1_next.truncate(min_dist - 2.7321 * dr);
				B2_prev.truncate(min_dist - 2.7321 * dr);
				B2_next.truncate(min_dist - 2.7321 * dr);
				B3_prev.truncate(min_dist - 2.7321 * dr);
				B3_next.truncate(min_dist - 2.7321 * dr);
			}

			if (buffer_dist <= max_dist + 1.7321 * dr)
			{
				buffer_dist = max_dist + (READAHEAD + 1.7321) * dr;
				if (loadHealpixData(&phi_prev, min_dist - 1.7321 * dr, buffer_dist) < 0)
				{
					cout << COLORTEXT_YELLOW << " warning" << COLORTEXT_RESET << ": read-ahead failed on phi_prev! step " << stepcount << ", available distance range = ";
					if (phi_prev.healpix_data.empty())
						cout << "none";
					else
						cout << phi_prev.healpix_data.begin()->second.hdr.distance << "-" << phi_prev.healpix_data.rbegin()->second.hdr.distance;
					cout << ", required distance range = " << min_dist-dr << "-" << max_dist+dr << endl;
					buffer_dist = max_dist + 1.7321 * dr;
				}
				if (loadHealpixData(&phi_next, min_dist - 1.7321 * dr, buffer_dist) < 0)
				{
					cout << COLORTEXT_YELLOW << " warning" << COLORTEXT_RESET << ": read-ahead failed on phi_next! step " << stepcount << ", available distance range = ";
					if (phi_next.healpix_data.empty())
						cout << "none";
					else
						cout << phi_next.healpix_data.begin()->second.hdr.distance << "-" << phi_next.healpix_data.rbegin()->second.hdr.distance;
					cout << ", required distance range = " << min_dist-dr << "-" << max_dist+dr << endl;
					buffer_dist = max_dist + 1.7321 * dr;
				}

				if (lightcone_fields & MASK_B)
				{
					if (loadHealpixData(&B1_prev, min_dist - 1.7321 * dr, buffer_dist) < 0)
					{
						cout << COLORTEXT_YELLOW << " warning" << COLORTEXT_RESET << ": read-ahead failed on B1_prev! step " << stepcount << ", available distance range = ";
						if (B1_prev.healpix_data.empty())
							cout << "none";
						else
							cout << B1_prev.healpix_data.begin()->second.hdr.distance << "-" << B1_prev.healpix_data.rbegin()->second.hdr.distance;
						cout << ", required distance range = " << min_dist-dr << "-" << max_dist+dr << endl;
						buffer_dist = max_dist + 1.7321 * dr;
					}
					if (loadHealpixData(&B1_next, min_dist - 1.7321 * dr, buffer_dist) < 0)
					{	
						cout << COLORTEXT_YELLOW << " warning" << COLORTEXT_RESET << ": read-ahead failed on B1_next! step " << stepcount << ", available distance range = ";
						if (B1_next.healpix_data.empty())
							cout << "none";
						else
							cout << B1_next.healpix_data.begin()->second.hdr.distance << "-" << B1_next.healpix_data.rbegin()->second.hdr.distance;
						cout << ", required distance range = " << min_dist-dr << "-" << max_dist+dr << endl;
						buffer_dist = max_dist + 1.7321 * dr;
					}
					if (loadHealpixData(&B2_prev, min_dist - 1.7321 * dr, buffer_dist) < 0)
					{
						cout << COLORTEXT_YELLOW << " warning" << COLORTEXT_RESET << ": read-ahead failed on B2_prev! step " << stepcount << ", available distance range = ";
						if (B2_prev.healpix_data.empty())
							cout << "none";
						else
							cout << B2_prev.healpix_data.begin()->second.hdr.distance << "-" << B2_prev.healpix_data.rbegin()->second.hdr.distance;
						cout << ", required distance range = " << min_dist-dr << "-" << max_dist+dr << endl;
						buffer_dist = max_dist + 1.7321 * dr;
					}
					if (loadHealpixData(&B2_next, min_dist - 1.7321 * dr, buffer_dist) < 0)
					{	
						cout << COLORTEXT_YELLOW << " warning" << COLORTEXT_RESET << ": read-ahead failed on B2_next! step " << stepcount << ", available distance range = ";
						if (B2_next.healpix_data.empty())
							cout << "none";
						else
							cout << B2_next.healpix_data.begin()->second.hdr.distance << "-" << B2_next.healpix_data.rbegin()->second.hdr.distance;
						cout << ", required distance range = " << min_dist-dr << "-" << max_dist+dr << endl;
						buffer_dist = max_dist + 1.7321 * dr;
					}
					if (loadHealpixData(&B3_prev, min_dist - 1.7321 * dr, buffer_dist) < 0)
					{
						cout << COLORTEXT_YELLOW << " warning" << COLORTEXT_RESET << ": read-ahead failed on B3_prev! step " << stepcount << ", available distance range = ";
						if (B3_prev.healpix_data.empty())
							cout << "none";
						else
							cout << B3_prev.healpix_data.begin()->second.hdr.distance << "-" << B3_prev.healpix_data.rbegin()->second.hdr.distance;
						cout << ", required distance range = " << min_dist-dr << "-" << max_dist+dr << endl;
						buffer_dist = max_dist + 1.7321 * dr;
					}
					if (loadHealpixData(&B3_next, min_dist - 1.7321 * dr, buffer_dist) < 0)
					{	
						cout << COLORTEXT_YELLOW << " warning" << COLORTEXT_RESET << ": read-ahead failed on B3_next! step " << stepcount << ", available distance range = ";
						if (B3_next.healpix_data.empty())
							cout << "none";
						else
							cout << B3_next.healpix_data.begin()->second.hdr.distance << "-" << B3_next.healpix_data.rbegin()->second.hdr.distance;
						cout << ", required distance range = " << min_dist-dr << "-" << max_dist+dr << endl;
						buffer_dist = max_dist + 1.7321 * dr;
					}
				}
			}

			min_dist = 1e30;
			max_dist = 0;

#pragma omp parallel for private(phi,phiprime,gradphi,tidalphi,ngradphi,nntidalphi,e1gradphi,e2gradphi,eetidalphi,exp2phi,B,gradB1,gradB2,gradB3,gradnB,grade1B,grade2B,nngradB,gradBprime,tidalB,tidalBterm,eetidalBterm,v,tmp,dist,source_dist,dn_dtau,dlnk0_dtau,filename,outfile), reduction(+:num_live), reduction(min:min_dist), reduction(max:max_dist)
			for (int n = 0; n < Nsource; n++)
			{
				if (photon[n].stop) continue;

				dist = get_direction(photon[n].pos, v, coordsys);

				if (!interpolation(&phi_prev, &phi_next, weight, dist, dr, v, &phi) || !time_derivative(&phi_prev, &phi_next, cycle_next.tau - cycle_prev.tau, dist, dr, v, &phiprime) || !gradient(&phi_prev, &phi_next, weight, dist, dr, v, gradphi) || !tidal_matrix(&phi_prev, &phi_next, weight, dist, dr, v, tidalphi, (n+Nsource_skip) % 500000 ? 0 : 1))
				{
					if ((n+Nsource_skip) % 500000 == 0)
					{
						cout << " photon " << n+Nsource_skip << " leaving data region at (" << photon[n].pos[0] << ", " << photon[n].pos[1] << ", " << photon[n].pos[2] << ") with n = (" << photon[n].n[0] << ", " << photon[n].n[1] << ", " << photon[n].n[2] << ")" << endl;
					}
					photon[n].stop = -1;
					continue;
				}

				//if (lightcone_fields & MASK_B)
				//{
				if (lightcone_fields & MASK_B == 0 || !interpolation(&B1_prev, &B1_next, weight, dist, dr, v, B) || !gradient(&B1_prev, &B1_next, weight, dist, dr, v, gradB1) || !interpolation(&B2_prev, &B2_next, weight, dist, dr, v, B+1) || !gradient(&B2_prev, &B2_next, weight, dist, dr, v, gradB2) || !interpolation(&B3_prev, &B3_next, weight, dist, dr, v, B+2) || !gradient(&B3_prev, &B3_next, weight, dist, dr, v, gradB3) || !time_derivative_gradient(&B1_prev, &B1_next, cycle_next.tau - cycle_prev.tau, dist, dr, v, gradBprime) || !tidal_matrix(&B1_prev, &B1_next, weight, dist, dr, v, tidalB, 0))
				/*	{
						photon[n].stop = -1;
						continue;
					}
				}
				else */
				{
					B[0]      = 0.; B[1]      = 0.; B[2] = 0.;
					gradB1[0] = 0.; gradB1[1] = 0.; gradB1[2] = 0.;
					gradB2[0] = 0.; gradB2[1] = 0.; gradB2[2] = 0.;
					gradB3[0] = 0.; gradB3[1] = 0.; gradB3[2] = 0.;
					grade1B[0]    = 0.; grade1B[1]    = 0.; grade1B[2]    = 0.;
					grade2B[0]    = 0.; grade2B[1]    = 0.; grade2B[2]    = 0.;
					gradBprime[0] = 0.; gradBprime[1] = 0.; gradBprime[2] = 0.;
					tidalB[0] = 0.; tidalB[1] = 0.; tidalB[2] = 0.;
					tidalB[3] = 0.; tidalB[4] = 0.; tidalB[5] = 0.;
					tidalBterm = 0.; eetidalBterm[0] = 0.; eetidalBterm[1] = 0.;
					gradnB[0] = 0.; gradnB[1] = 0.; gradnB[2] = 0.;
				}
				
				exp2phi = exp(2 * phi);

				tidalphi[0] += 2 * gradphi[0] * gradphi[0];
				tidalphi[1] += 2 * gradphi[0] * gradphi[1];
				tidalphi[2] += 2 * gradphi[0] * gradphi[2];
				tidalphi[3] += 2 * gradphi[1] * gradphi[1];
				tidalphi[4] += 2 * gradphi[1] * gradphi[2];
				tidalphi[5] += 2 * gradphi[2] * gradphi[2];

				/*photon[n].e1[0] = atan2(photon[n].n[1], photon[n].n[0]); // FIXME
				photon[n].e2[0] = -sin(photon[n].e1[0]);
				photon[n].e2[1] = cos(photon[n].e1[0]);
				photon[n].e2[2] = 0;
				photon[n].e1[0] = photon[n].n[2] * photon[n].e2[1];
				photon[n].e1[1] = -photon[n].n[2] * photon[n].e2[0];
				//photon[n].e1[2] = -sin(acos(photon[n].n[2]));
				photon[n].e1[2] = -sqrt(1.-photon[n].n[2]*photon[n].n[2]);*/

				rotate_vector(photon[n].e1, tmp, coordsys);
				rotate_vector(photon[n].e2, dn_dtau, coordsys); // dn_dtau is a temporary array here...

				/*for (int i = 0; i < 3; i++) // FIXME
				{
					photon[n].e1[i] = tmp[i];
					photon[n].e2[i] = dn_dtau[i];
				}*/

				e1gradphi = tmp[0] * gradphi[0] + tmp[1] * gradphi[1] + tmp[2] * gradphi[2];
				e2gradphi = dn_dtau[0] * gradphi[0] + dn_dtau[1] * gradphi[1] + dn_dtau[2] * gradphi[2];

				eetidalphi[0] = (tmp[0]*tmp[0] - dn_dtau[0]*dn_dtau[0])*tidalphi[0] + 2*(tmp[0]*tmp[1] - dn_dtau[0]*dn_dtau[1])*tidalphi[1] + 2*(tmp[0]*tmp[2] - dn_dtau[0]*dn_dtau[2])*tidalphi[2] + (tmp[1]*tmp[1] - dn_dtau[1]*dn_dtau[1])*tidalphi[3] + 2*(tmp[2]*tmp[1] - dn_dtau[2]*dn_dtau[1])*tidalphi[4] + (tmp[2]*tmp[2] - dn_dtau[2]*dn_dtau[2])*tidalphi[5];
				eetidalphi[1] = 2*tmp[0]*dn_dtau[0]*tidalphi[0] + 2*(tmp[0]*dn_dtau[1] + tmp[1]*dn_dtau[0])*tidalphi[1] + 2*(tmp[0]*dn_dtau[2] + tmp[2]*dn_dtau[0])*tidalphi[2] + 2*tmp[1]*dn_dtau[1]*tidalphi[3] + 2*(tmp[2]*dn_dtau[1] + tmp[1]*dn_dtau[2])*tidalphi[4] + 2*tmp[2]*dn_dtau[2]*tidalphi[5];


				if (lightcone_fields & MASK_B)
				{
					rotate_vector(photon[n].n, gradnB, coordsys); // gradnB is temporary array here!

					tidalBterm = photon[n].n[0] * (tidalB[0] + tidalB[3] + tidalB[5] + gradBprime[0] * gradnB[0] + gradBprime[1] * gradnB[1] + gradBprime[2] * gradnB[2]);
					eetidalBterm[0] = (tmp[0]*tmp[0] - dn_dtau[0]*dn_dtau[0])*(photon[n].n[0]*tidalB[0] - gradBprime[0] - tidalB[0]*gradnB[0] - tidalB[1]*gradnB[1] - tidalB[2]*gradnB[2]) + (tmp[0]*tmp[1] - dn_dtau[0]*dn_dtau[1])*(2*photon[n].n[0]*tidalB[1] - gradBprime[1] - tidalB[1]*gradnB[0] - tidalB[3]*gradnB[1] - tidalB[4]*gradnB[2]) + (tmp[0]*tmp[2] - dn_dtau[0]*dn_dtau[2])*(2*photon[n].n[0]*tidalB[2] - gradBprime[2] - tidalB[2]*gradnB[0] - tidalB[4]*gradnB[1] - tidalB[5]*gradnB[2]) + (tmp[1]*tmp[1] - dn_dtau[1]*dn_dtau[1])*(photon[n].n[0]*tidalB[3]) + (tmp[2]*tmp[1] - dn_dtau[2]*dn_dtau[1])*(2*photon[n].n[0]*tidalB[4]) + (tmp[2]*tmp[2] - dn_dtau[2]*dn_dtau[2])*(photon[n].n[0]*tidalB[5]);
					eetidalBterm[1] = 2*tmp[0]*dn_dtau[0]*(photon[n].n[0]*tidalB[0] - gradBprime[0] - tidalB[0]*gradnB[0] - tidalB[1]*gradnB[1] - tidalB[2]*gradnB[2]) + (tmp[0]*dn_dtau[1] + tmp[1]*dn_dtau[0])*(2*photon[n].n[0]*tidalB[1] - gradBprime[1] - tidalB[1]*gradnB[0] - tidalB[3]*gradnB[1] - tidalB[4]*gradnB[2]) + (tmp[0]*dn_dtau[2] + tmp[2]*dn_dtau[0])*(2*photon[n].n[0]*tidalB[2] - gradBprime[2] - tidalB[2]*gradnB[0] - tidalB[4]*gradnB[1] - tidalB[5]*gradnB[2]) + 2*tmp[1]*dn_dtau[1]*(photon[n].n[0]*tidalB[3]) + (tmp[2]*dn_dtau[1] + tmp[1]*dn_dtau[2])*(2*photon[n].n[0]*tidalB[4]) + 2*tmp[2]*dn_dtau[2]*(photon[n].n[0]*tidalB[5]);

					if(!time_derivative_gradient(&B2_prev, &B2_next, cycle_next.tau - cycle_prev.tau, dist, dr, v, gradBprime) || !tidal_matrix(&B2_prev, &B2_next, weight, dist, dr, v, tidalB, 0))
					{
						photon[n].stop = -1;
						continue;
					}

					tidalBterm += photon[n].n[1] * (tidalB[0] + tidalB[3] + tidalB[5] + gradBprime[0] * gradnB[0] + gradBprime[1] * gradnB[1] + gradBprime[2] * gradnB[2]);
					eetidalBterm[0] += (tmp[0]*tmp[0] - dn_dtau[0]*dn_dtau[0])*(photon[n].n[1]*tidalB[0]) + (tmp[0]*tmp[1] - dn_dtau[0]*dn_dtau[1])*(2*photon[n].n[1]*tidalB[1] - gradBprime[0] - tidalB[0]*gradnB[0] - tidalB[1]*gradnB[1] - tidalB[2]*gradnB[2]) + (tmp[0]*tmp[2] - dn_dtau[0]*dn_dtau[2])*(2*photon[n].n[1]*tidalB[2]) + (tmp[1]*tmp[1] - dn_dtau[1]*dn_dtau[1])*(photon[n].n[1]*tidalB[3] - gradBprime[1] - tidalB[1]*gradnB[0] - tidalB[3]*gradnB[1] - tidalB[4]*gradnB[2]) + (tmp[2]*tmp[1] - dn_dtau[2]*dn_dtau[1])*(2*photon[n].n[1]*tidalB[4] - gradBprime[2] - tidalB[2]*gradnB[0] - tidalB[4]*gradnB[1] - tidalB[5]*gradnB[2]) + (tmp[2]*tmp[2] - dn_dtau[2]*dn_dtau[2])*(photon[n].n[1]*tidalB[5]);
					eetidalBterm[1] += 2*tmp[0]*dn_dtau[0]*(photon[n].n[1]*tidalB[0]) + (tmp[0]*dn_dtau[1] + tmp[1]*dn_dtau[0])*(2*photon[n].n[1]*tidalB[1] - gradBprime[0] - tidalB[0]*gradnB[0] - tidalB[1]*gradnB[1] - tidalB[2]*gradnB[2]) + (tmp[0]*dn_dtau[2] + tmp[2]*dn_dtau[0])*(2*photon[n].n[1]*tidalB[2]) + 2*tmp[1]*dn_dtau[1]*(photon[n].n[1]*tidalB[3] - gradBprime[1] - tidalB[1]*gradnB[0] - tidalB[3]*gradnB[1] - tidalB[4]*gradnB[2]) + (tmp[2]*dn_dtau[1] + tmp[1]*dn_dtau[2])*(2*photon[n].n[1]*tidalB[4] - gradBprime[2] - tidalB[2]*gradnB[0] - tidalB[4]*gradnB[1] - tidalB[5]*gradnB[2]) + 2*tmp[2]*dn_dtau[2]*(photon[n].n[1]*tidalB[5]);

					if(!time_derivative_gradient(&B3_prev, &B3_next, cycle_next.tau - cycle_prev.tau, dist, dr, v, gradBprime) || !tidal_matrix(&B3_prev, &B3_next, weight, dist, dr, v, tidalB, 0))
					{
						photon[n].stop = -1;
						continue;
					}

					tidalBterm += photon[n].n[2] * (tidalB[0] + tidalB[3] + tidalB[5] + gradBprime[0] * gradnB[0] + gradBprime[1] * gradnB[1] + gradBprime[2] * gradnB[2]);
					eetidalBterm[0] += (tmp[0]*tmp[0] - dn_dtau[0]*dn_dtau[0])*(photon[n].n[2]*tidalB[0]) + (tmp[0]*tmp[1] - dn_dtau[0]*dn_dtau[1])*(2*photon[n].n[2]*tidalB[1]) + (tmp[0]*tmp[2] - dn_dtau[0]*dn_dtau[2])*(2*photon[n].n[2]*tidalB[2] - gradBprime[0] - tidalB[0]*gradnB[0] - tidalB[1]*gradnB[1] - tidalB[2]*gradnB[2]) + (tmp[1]*tmp[1] - dn_dtau[1]*dn_dtau[1])*(photon[n].n[2]*tidalB[3]) + (tmp[2]*tmp[1] - dn_dtau[2]*dn_dtau[1])*(2*photon[n].n[2]*tidalB[4] - gradBprime[1] - tidalB[1]*gradnB[0] - tidalB[3]*gradnB[1] - tidalB[4]*gradnB[2]) + (tmp[2]*tmp[2] - dn_dtau[2]*dn_dtau[2])*(photon[n].n[2]*tidalB[5] - gradBprime[2] - tidalB[2]*gradnB[0] - tidalB[4]*gradnB[1] - tidalB[5]*gradnB[2]);
					eetidalBterm[1] += 2*tmp[0]*dn_dtau[0]*(photon[n].n[2]*tidalB[0]) + (tmp[0]*dn_dtau[1] + tmp[1]*dn_dtau[0])*(2*photon[n].n[2]*tidalB[1]) + (tmp[0]*dn_dtau[2] + tmp[2]*dn_dtau[0])*(2*photon[n].n[2]*tidalB[2] - gradBprime[0] - tidalB[0]*gradnB[0] - tidalB[1]*gradnB[1] - tidalB[2]*gradnB[2]) + 2*tmp[1]*dn_dtau[1]*(photon[n].n[2]*tidalB[3]) + (tmp[2]*dn_dtau[1] + tmp[1]*dn_dtau[2])*(2*photon[n].n[2]*tidalB[4] - gradBprime[1] - tidalB[1]*gradnB[0] - tidalB[3]*gradnB[1] - tidalB[4]*gradnB[2]) + 2*tmp[2]*dn_dtau[2]*(photon[n].n[2]*tidalB[5] - gradBprime[2] - tidalB[2]*gradnB[0] - tidalB[4]*gradnB[1] - tidalB[5]*gradnB[2]);

					for (int i = 0; i < 3; i++)
						gradnB[i] = -photon[n].e1[0] * gradB1[i] - photon[n].e1[1] * gradB2[i] - photon[n].e1[2] * gradB3[i];

					counterrotate_vector(gradnB, grade1B, coordsys);

					for (int i = 0; i < 3; i++)
						gradnB[i] = -photon[n].e2[0] * gradB1[i] - photon[n].e2[1] * gradB2[i] - photon[n].e2[2] * gradB3[i];

					counterrotate_vector(gradnB, grade2B, coordsys);
	
					for (int i = 0; i < 3; i++)
						gradnB[i] = photon[n].n[0] * gradB1[i] + photon[n].n[1] * gradB2[i] + photon[n].n[2] * gradB3[i];

					for (int i = 0; i < 3; i++)
					{
						grade1B[0] += tmp[i] * (photon[n].n[0] * gradnB[i] + gradB1[i]);
						grade1B[1] += tmp[i] * (photon[n].n[1] * gradnB[i] + gradB2[i]);
						grade1B[2] += tmp[i] * (photon[n].n[2] * gradnB[i] + gradB3[i]);
						grade2B[0] += dn_dtau[i] * (photon[n].n[0] * gradnB[i] + gradB1[i]);
						grade2B[1] += dn_dtau[i] * (photon[n].n[1] * gradnB[i] + gradB2[i]);
						grade2B[2] += dn_dtau[i] * (photon[n].n[2] * gradnB[i] + gradB3[i]);
					}
				}

				rotate_vector(photon[n].n, tmp, coordsys);

				ngradphi = tmp[0] * gradphi[0] + tmp[1] * gradphi[1] + tmp[2] * gradphi[2];

				nntidalphi = tmp[0]*tmp[0]*tidalphi[0] + 2*tmp[0]*tmp[1]*tidalphi[1] + 2*tmp[0]*tmp[2]*tidalphi[2] + tmp[1]*tmp[1]*tidalphi[3] + 2*tmp[1]*tmp[2]*tidalphi[4] + tmp[2]*tmp[2]*tidalphi[5];

				for (int i = 0; i < 3; i++)
				{
					grade1B[0] += photon[n].n[0] * tmp[i] * (photon[n].e1[0] * gradB1[i] + photon[n].e1[1] * gradB2[i] + photon[n].e1[2] * gradB3[i]);
					grade1B[1] += photon[n].n[1] * tmp[i] * (photon[n].e1[0] * gradB1[i] + photon[n].e1[1] * gradB2[i] + photon[n].e1[2] * gradB3[i]);
					grade1B[2] += photon[n].n[2] * tmp[i] * (photon[n].e1[0] * gradB1[i] + photon[n].e1[1] * gradB2[i] + photon[n].e1[2] * gradB3[i]);
					grade2B[0] += photon[n].n[0] * tmp[i] * (photon[n].e2[0] * gradB1[i] + photon[n].e2[1] * gradB2[i] + photon[n].e2[2] * gradB3[i]);
					grade2B[1] += photon[n].n[1] * tmp[i] * (photon[n].e2[0] * gradB1[i] + photon[n].e2[1] * gradB2[i] + photon[n].e2[2] * gradB3[i]);
					grade2B[2] += photon[n].n[2] * tmp[i] * (photon[n].e2[0] * gradB1[i] + photon[n].e2[1] * gradB2[i] + photon[n].e2[2] * gradB3[i]);
				}

				nngradB = tmp[0] * gradnB[0] + tmp[1] * gradnB[1] + tmp[2] * gradnB[2];

				for (int i = 0; i < 3; i++)
					tmp[i] = (tmp[i] * ngradphi - gradphi[i]) * 2 * exp2phi + tmp[i] * nngradB - gradnB[i];

				counterrotate_vector(tmp, dn_dtau, coordsys);

				dlnk0_dtau = -4. * ngradphi * exp2phi - 2. * phiprime - nngradB;

				if ((n+Nsource_skip) % 500000 == 0)
#pragma omp critical(evala)
				{
					rotate_vector(photon[n].e1, tmp, coordsys);
					photon[n].e1[0] = tmp[0]; photon[n].e1[1] = tmp[1]; photon[n].e1[2] = tmp[2];
					rotate_vector(photon[n].e2, tmp, coordsys);
					photon[n].e2[0] = tmp[0]; photon[n].e2[1] = tmp[1]; photon[n].e2[2] = tmp[2];
					counterrotate_vector(gradphi, tmp, coordsys);
					sprintf(filename, "ray%03d-pass%d.dat", (n+Nsource_skip)/500000, pass);
					outfile = fopen(filename, "a");
					fprintf(outfile, " %f  %f  %f  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e\n", photon[n].pos[0], photon[n].pos[1], photon[n].pos[2], exp(photon[n].lnk0) / gsl_spline_eval (aspline, tau, acc) / gsl_spline_eval (aspline, tau, acc), phi, photon[n].DA * gsl_spline_eval (aspline, tau, acc), photon[n].sigma[0], photon[n].sigma[1], photon[n].ellip[0], photon[n].ellip[1], tmp[0], tmp[1], tmp[2], nntidalphi, tidalphi[0] + tidalphi[3] + tidalphi[5], photon[n].e1[0]*photon[n].e1[0]*tidalphi[0] + 2*photon[n].e1[0]*photon[n].e1[1]*tidalphi[1] + 2*photon[n].e1[0]*photon[n].e1[2]*tidalphi[2] + photon[n].e1[1]*photon[n].e1[1]*tidalphi[3] + 2*photon[n].e1[1]*photon[n].e1[2]*tidalphi[4] + photon[n].e1[2]*photon[n].e1[2]*tidalphi[5], photon[n].e2[0]*photon[n].e2[0]*tidalphi[0] + 2*photon[n].e2[0]*photon[n].e2[1]*tidalphi[1] + 2*photon[n].e2[0]*photon[n].e2[2]*tidalphi[2] + photon[n].e2[1]*photon[n].e2[1]*tidalphi[3] + 2*photon[n].e2[1]*photon[n].e2[2]*tidalphi[4] + photon[n].e2[2]*photon[n].e2[2]*tidalphi[5], photon[n].e1[0]*photon[n].e2[0]*tidalphi[0] + photon[n].e1[0]*photon[n].e2[1]*tidalphi[1] + photon[n].e2[0]*photon[n].e1[1]*tidalphi[1] + photon[n].e1[0]*photon[n].e2[2]*tidalphi[2] + photon[n].e2[0]*photon[n].e1[2]*tidalphi[2] + photon[n].e1[1]*photon[n].e2[1]*tidalphi[3] + photon[n].e1[1]*photon[n].e2[2]*tidalphi[4] + photon[n].e2[1]*photon[n].e1[2]*tidalphi[4] + photon[n].e1[2]*photon[n].e2[2]*tidalphi[5], photon[n].n[0]*photon[n].n[0] + photon[n].n[1]*photon[n].n[1] + photon[n].n[2]*photon[n].n[2] - 1., photon[n].e1[0]*photon[n].e1[0] + photon[n].e1[1]*photon[n].e1[1] + photon[n].e1[2]*photon[n].e1[2] - 1., photon[n].e2[0]*photon[n].e2[0] + photon[n].e2[1]*photon[n].e2[1] + photon[n].e2[2]*photon[n].e2[2] - 1., photon[n].e1[0]*photon[n].e2[0] + photon[n].e1[1]*photon[n].e2[1] + photon[n].e1[2]*photon[n].e2[2], B[0], B[1], B[2]);
					fclose(outfile);
					counterrotate_vector(photon[n].e1, tmp, coordsys);
					photon[n].e1[0] = tmp[0]; photon[n].e1[1] = tmp[1]; photon[n].e1[2] = tmp[2];
					counterrotate_vector(photon[n].e2, tmp, coordsys);
					photon[n].e2[0] = tmp[0]; photon[n].e2[1] = tmp[1]; photon[n].e2[2] = tmp[2];
				}

				photon[n].DAprime *= (1. - 0.5 * integration_step * dlnk0_dtau);

				if (photon[n].DA) // treat first time step separately
				{
					photon[n].sigma[0] = (photon[n].sigma[0] - 0.5 * integration_step * photon[n].DA * photon[n].DA * (eetidalphi[0] * exp2phi * exp2phi + 0.5 * eetidalBterm[0])) / (1. + 0.5 * integration_step * dlnk0_dtau);
					photon[n].sigma[1] = (photon[n].sigma[1] - 0.5 * integration_step * photon[n].DA * photon[n].DA * (eetidalphi[1] * exp2phi * exp2phi + 0.5 * eetidalBterm[1])) / (1. + 0.5 * integration_step * dlnk0_dtau);

					photon[n].DAprime = (photon[n].DAprime + integration_step * photon[n].DA * ((nntidalphi - tidalphi[0] - tidalphi[3] - tidalphi[5]) * exp2phi * exp2phi - 0.5 * tidalBterm) - integration_step * (photon[n].sigma[0]*photon[n].sigma[0] + photon[n].sigma[1]*photon[n].sigma[1]) / photon[n].DA / photon[n].DA / photon[n].DA) / (1. + 0.5 * integration_step * dlnk0_dtau);
				}

				photon[n].sigma[0] = photon[n].sigma[0] * (1. - 0.5 * integration_step * dlnk0_dtau) - 0.5 * integration_step * photon[n].DA * photon[n].DA * (eetidalphi[0] * exp2phi * exp2phi + 0.5 * eetidalBterm[0]);
				photon[n].sigma[1] = photon[n].sigma[1] * (1. - 0.5 * integration_step * dlnk0_dtau) - 0.5 * integration_step * photon[n].DA * photon[n].DA * (eetidalphi[1] * exp2phi * exp2phi + 0.5 * eetidalBterm[1]);

				photon[n].DA += 0.5 * integration_step * photon[n].DAprime;

				tmp[2] = photon[n].ellip[0]*photon[n].ellip[0] + photon[n].ellip[1]*photon[n].ellip[1];
				tmp[0] = photon[n].ellip[0] + 2. * integration_step * photon[n].sigma[0] * sqrt(tmp[2] + 4.) / photon[n].DA / photon[n].DA;
				tmp[1] = photon[n].ellip[1] + 2. * integration_step * photon[n].sigma[1] * sqrt(tmp[2] + 4.) / photon[n].DA / photon[n].DA;
				tmp[2] = 0.5 * (tmp[2] + tmp[0]*tmp[0] + tmp[1]*tmp[1]);

				photon[n].ellip[0] += integration_step * photon[n].sigma[0] * sqrt(tmp[2] + 4.) / photon[n].DA / photon[n].DA;
				photon[n].ellip[1] += integration_step * photon[n].sigma[1] * sqrt(tmp[2] + 4.) / photon[n].DA / photon[n].DA;
				
				photon[n].omega += integration_step * (photon[n].sigma[1] * photon[n].ellip[0] - photon[n].sigma[0] * photon[n].ellip[1]) / photon[n].DA / photon[n].DA / (2. + sqrt(tmp[2] + 4.));
				
				photon[n].ellip[0] += integration_step * photon[n].sigma[0] * sqrt(tmp[2] + 4.) / photon[n].DA / photon[n].DA;
				photon[n].ellip[1] += integration_step * photon[n].sigma[1] * sqrt(tmp[2] + 4.) / photon[n].DA / photon[n].DA;

				photon[n].DA += 0.5 * integration_step * photon[n].DAprime;

				for (int i = 0; i < 3; i++)
				{
					photon[n].e1[i] += integration_step * (photon[n].n[i] * 2 * exp2phi * e1gradphi + 0.5 * grade1B[i]);
					photon[n].e2[i] += integration_step * (photon[n].n[i] * 2 * exp2phi * e2gradphi + 0.5 * grade2B[i]);
					photon[n].n[i] += integration_step * dn_dtau[i];
					photon[n].pos[i] += integration_step * (photon[n].n[i] * exp2phi + B[i]);
				}

				photon[n].lnk0 += integration_step * dlnk0_dtau;

				dist = sqrt(photon[n].pos[0]*photon[n].pos[0] + photon[n].pos[1]*photon[n].pos[1] + photon[n].pos[2]*photon[n].pos[2]);
				source_dist = sqrt(source_pos[3*n]*source_pos[3*n] + source_pos[1+3*n]*source_pos[1+3*n] + source_pos[2+3*n]*source_pos[2+3*n]);

				if (dist > source_dist)
				{
					photon[n].stop = 1;

					source_dist = tau_obs - tau - integration_step - dist - photon[n].delay;
					photon[n].delay = tau_obs - tau - integration_step - dist;

					for (int i = 0; i < 3; i++)
						source_pos[i+3*n] -= source_dist * source_vel[i+3*n];
			
					source_dist = sqrt(source_pos[3*n]*source_pos[3*n] + source_pos[1+3*n]*source_pos[1+3*n] + source_pos[2+3*n]*source_pos[2+3*n]);

					source_dist = (source_dist - dist) * source_dist / (source_pos[3*n] * photon[n].n[0] + source_pos[1+3*n] * photon[n].n[1] + source_pos[2+3*n] * photon[n].n[2]) / exp2phi; // positive time step

					for (int i = 0; i < 3; i++)
						photon[n].pos[i] += source_dist * (photon[n].n[i] * exp2phi + B[i]);

					photon[n].lnk0 += source_dist * dlnk0_dtau;

					photon[n].DA += source_dist * photon[n].DAprime;

					photon[n].ellip[0] += source_dist * photon[n].sigma[0] * sqrt(tmp[2] + 4.) / photon[n].DA / photon[n].DA;
					photon[n].ellip[1] += source_dist * photon[n].sigma[1] * sqrt(tmp[2] + 4.) / photon[n].DA / photon[n].DA;
					
					photon[n].omega += source_dist * (photon[n].sigma[1] * photon[n].ellip[0] - photon[n].sigma[0] * photon[n].ellip[1]) / photon[n].DA / photon[n].DA / (2. + sqrt(tmp[2] + 4.));
				
					photon[n].ellip[0] += source_dist * photon[n].sigma[0] * sqrt(tmp[2] + 4.) / photon[n].DA / photon[n].DA;
					photon[n].ellip[1] += source_dist * photon[n].sigma[1] * sqrt(tmp[2] + 4.) / photon[n].DA / photon[n].DA;

					if ((n+Nsource_skip) % 500000 == 0)
					{
						cout << " photon " << n+Nsource_skip << " reaching source position at (" << photon[n].pos[0] << ", " << photon[n].pos[1] << ", " << photon[n].pos[2] << ") with n = (" << photon[n].n[0] << ", " << photon[n].n[1] << ", " << photon[n].n[2] << ") and delay of " << photon[n].delay * sim.boxsize << " Mpc/h";
						if (pass == 1) cout << "; the source position is adjusted by (" << source_vel[3*n] * photon[n].delay * sim.boxsize << ", " << source_vel[1+3*n] * photon[n].delay * sim.boxsize << ", " << source_vel[2+3*n] * photon[n].delay * sim.boxsize << ") Mpc/h";
						cout << endl;
					}

					dist = get_direction(photon[n].pos, v, coordsys);

					interpolation(&phi_prev, &phi_next, (tau + integration_step + source_dist - cycle_prev.tau) / (cycle_next.tau - cycle_prev.tau), dist, dr, v, &phi);

#pragma omp critical(evala)
{
					photon[n].lnk0 = exp(photon[n].lnk0 + 4. * phi) * (sqrt(exp(-2. * phi) + source_vel[3*n]*source_vel[3*n] + source_vel[1+3*n]*source_vel[1+3*n] + source_vel[2+3*n]*source_vel[2+3*n]) - (photon[n].n[0]*source_vel[3*n] + photon[n].n[1]*source_vel[1+3*n] + photon[n].n[2]*source_vel[2+3*n])) / gsl_spline_eval (aspline, tau + integration_step + source_dist, acc);
					photon[n].DA *= gsl_spline_eval (aspline, tau + integration_step + source_dist, acc) * exp(-phi);
}
				}
				else
				{
					if (dist < min_dist) min_dist = dist;
					if (dist > max_dist) max_dist = dist;
				}

				if (stepcount % 1000 == 0)
				{
					tmp[0] = sqrt(photon[n].n[0] * photon[n].n[0] + photon[n].n[1] * photon[n].n[1] + photon[n].n[2] * photon[n].n[2]);
					tmp[1] = sqrt(photon[n].e1[0] * photon[n].e1[0] + photon[n].e1[1] * photon[n].e1[1] + photon[n].e1[2] * photon[n].e1[2]);
					tmp[2] = sqrt(photon[n].e2[0] * photon[n].e2[0] + photon[n].e2[1] * photon[n].e2[1] + photon[n].e2[2] * photon[n].e2[2]);
					for (int i = 0; i < 3; i++)
					{
						photon[n].n[i] /= tmp[0];
						photon[n].e1[i] /= tmp[1];
						photon[n].e2[i] /= tmp[2];
					}
				}

				num_live++;
			}
			
			tau += integration_step;
			stepcount++;

			if (stepcount % 1000 == 0)
			{
				if (numtasks > 1)
				{
					sprintf(filename, "%s.%d", restart_filename, taskindex);
					cout << " task " << taskindex << " writing restart point to " << filename << endl << endl;
					outfile = fopen(filename, "wb");
				}
				else
				{
					cout << " writing restart point to " << restart_filename << endl << endl;
					outfile = fopen(restart_filename, "wb");
				}
				
				if(outfile == NULL)
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": could not open restart file for output!" << endl;
					return -1;
				}

				if (fwrite((void *) &tau, sizeof(double), 1, outfile) != 1 || fwrite((void *) &stepcount, sizeof(long), 1, outfile) != 1 || fwrite((void *) &buffer_dist, sizeof(double), 1, outfile) != 1 || fwrite((void *) &min_dist, sizeof(double), 1, outfile) != 1 || fwrite((void *) &max_dist, sizeof(double), 1, outfile) != 1)
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": could not write header to restart file!" << endl;
					fclose(outfile);
					return -1;
				}
				if (fwrite((void *) photon, sizeof(photon_data), Nsource, outfile) != Nsource)
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": could not write photon data to restart file!" << endl;
					fclose(outfile);
					return -1;
				}
				if (fwrite((void *) source_pos, sizeof(double), 3l*Nsource, outfile) != 3l*Nsource)
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": could not write source data to restart file!" << endl;
					fclose(outfile);
					return -1;
				}
				fclose(outfile);
			}
		}
		while (num_live > 0);

		phi_prev.clear();
		phi_next.clear();

		if (lightcone_fields & MASK_B)
		{
			B1_prev.clear();
			B1_next.clear();
			B2_prev.clear();
			B2_next.clear();
			B3_prev.clear();
			B3_next.clear();
		}

		cout << " task " << taskindex << " pass " << pass << " done, correcting shooting angles..." << endl << endl;

		maxcorr = 0;
		meancorr = 0;
		count = 0;

		for (int n = 0; n < Nsource; n++)
		{
			if (photon[n].stop < 0) continue;

			source_dist = sqrt(source_pos[3*n]*source_pos[3*n] + source_pos[1+3*n]*source_pos[1+3*n] + source_pos[2+3*n]*source_pos[2+3*n]);

			for (int i = 0; i < 3; i++)
				v[i] = (source_pos[i+3*n] - photon[n].pos[i]) / source_dist;

			v[2] = (v[0] * cos(source_direction[1+2*n]) + v[1] * sin(source_direction[1+2*n])) * cos(source_direction[2*n]) - v[2] * sin(source_direction[2*n]);
			v[3] = v[1] * cos(source_direction[1+2*n]) - v[0] * sin(source_direction[1+2*n]);

			dist = v[3]*v[3] + v[2]*v[2];

			if (dist > maxcorr) maxcorr = dist;
			meancorr += dist;
			count++;

			v[0] = (sin(source_direction[2*n]) + cos(source_direction[2*n]) * v[2]) * cos(source_direction[1+2*n]) - v[3] * sin(source_direction[1+2*n]);
			v[1] = (sin(source_direction[2*n]) + cos(source_direction[2*n]) * v[2]) * sin(source_direction[1+2*n]) + v[3] * cos(source_direction[1+2*n]);
			v[2] = cos(source_direction[2*n]) - sin(source_direction[2*n]) * v[2];

			source_direction[2*n] = acos(v[2] / sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]));
			source_direction[1+2*n] = atan2(v[1], v[0]);
		}

		cout << " task " << taskindex << " rms / maximum shooting angle correction was " << sqrt(meancorr/count) << " / " << sqrt(maxcorr) << " radians." << endl << endl;

		if (numtasks > 1)
		{
			sprintf(filename, "%s.%d", output_filename, taskindex);
			outfile = fopen(filename, "w");	
		}
		else
			outfile = fopen(output_filename, "w");

		if (outfile == NULL)
		{
			cout << " error opening file for output!" << endl;
		}
		else
		{
			fprintf(outfile, "# LCARS analysis for source catalog %s\n", catalog_filename);
			fprintf(outfile, "# output for pass %d\n", pass);
			fprintf(outfile, "# theta        phi          1+z_obs      D_A            Re(e)          Im(e)           omega           residual       mu'             phi'          alpha       cos(beta)     gamma\n");
			for (int n = 0; n < Nsource; n++)
			{
				if (photon[n].stop > 0)
					fprintf(outfile, "  %.8lf   %.8lf   %.8lf   %.10lf   %e   %e   %e   %e   %.10lf   %.10lf  %.8lf  %.8lf  %e\n", source_direction[2*n], (source_direction[1+2*n] < 0.) ? source_direction[1+2*n] + 2.*M_PI : source_direction[1+2*n], photon[n].lnk0, photon[n].DA, photon[n].ellip[0], photon[n].ellip[1], photon[n].omega, sqrt((photon[n].pos[0]-source_pos[3*n])*(photon[n].pos[0]-source_pos[3*n]) + (photon[n].pos[1]-source_pos[1+3*n])*(photon[n].pos[1]-source_pos[1+3*n]) + (photon[n].pos[2]-source_pos[2+3*n])*(photon[n].pos[2]-source_pos[2+3*n]))*sim.numpts, sin(source_direction[2*n])*(cos(source_direction[1+2*n])*coordsys[3]*coordsys[0] + sin(source_direction[1+2*n])*coordsys[3]*coordsys[1]) + cos(source_direction[2*n])*coordsys[2], atan2(sin(source_direction[2*n])*(coordsys[0]*sin(source_direction[1+2*n])-coordsys[1]*cos(source_direction[1+2*n])), sin(source_direction[2*n])*coordsys[2]*(coordsys[1]*sin(source_direction[1+2*n])+coordsys[0]*cos(source_direction[1+2*n])) - cos(source_direction[2*n])*coordsys[3]), atan2(-photon[n].n[1],-photon[n].n[0]), -photon[n].n[2], atan2(photon[n].e1[1]*photon[n].n[0]-photon[n].e1[0]*photon[n].n[1],photon[n].e2[1]*photon[n].n[0]-photon[n].e2[0]*photon[n].n[1]));
				else
					fprintf(outfile, "  %.8lf   %.8lf   -1.0000000   -1.000000000   0.000000e+00   0.000000e+00   0.000000e+00   -1.000000e+00   %.10lf   %.10lf  %.8lf  %.8lf  0.000000e+00\n", source_direction[2*n], (source_direction[1+2*n] < 0.) ? source_direction[1+2*n] + 2.*M_PI : source_direction[1+2*n], sin(source_direction[2*n])*(cos(source_direction[1+2*n])*coordsys[3]*coordsys[0] + sin(source_direction[1+2*n])*coordsys[3]*coordsys[1]) + cos(source_direction[2*n])*coordsys[2], atan2(sin(source_direction[2*n])*(coordsys[0]*sin(source_direction[1+2*n])-coordsys[1]*cos(source_direction[1+2*n])), sin(source_direction[2*n])*coordsys[2]*(coordsys[1]*sin(source_direction[1+2*n])+coordsys[0]*cos(source_direction[1+2*n])) - cos(source_direction[2*n])*coordsys[3]), (source_direction[1+2*n] < 0.) ? source_direction[1+2*n] + 2.*M_PI : source_direction[1+2*n], cos(source_direction[2*n]));
			}
			fclose(outfile);
		}

		if (numtasks > 1)
			cout << " task " << taskindex << " source data written to " << filename << endl << endl;
		else
			cout << " source data written to " << output_filename << endl << endl;

		tau = tau_obs;
		stepcount = 1;
		buffer_dist = READAHEAD * dr - integration_step;
		min_dist = 0;
		max_dist = 0;
	}

	free(photon);
	free(source_pos);
	free(source_vel);
	free(source_direction);

	gsl_spline_free(aspline);
	gsl_interp_accel_free(acc);

	cout << " task " << taskindex << " normal completion." << endl << endl;

	return 0;
}

int loadHealpixData(metric_container * field, double min_dist, double max_dist)
{
	metric_data metric;
	char filename[1024];
	uint32_t blocksize[2];
	FILE * infile = NULL;
	int count;
	long backtrack;

	metric.pixel = NULL;

	if (!field->healpix_data.empty())
	{
		if (field->healpix_data.begin()->second.hdr.distance < min_dist)
		{
			if ((min_dist = field->healpix_data.rbegin()->second.hdr.distance) > max_dist)
				return field->healpix_data.rbegin()->first;
		}	
		else if (field->healpix_data.begin()->second.hdr.distance > max_dist)
			max_dist = field->healpix_data.begin()->second.hdr.distance;
	}

	sprintf(filename, "%s%s_%04d_%s.map", field->dir, field->basename, field->cinfo.num, field->name);

	infile = fopen(filename, "rb");

	if (infile == NULL)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": map file " << filename << " could not be opened for reading!" << endl;
		return -1;
	}

	if (fread(blocksize+1, sizeof(uint32_t), 1, infile) != 1)
	{
		fclose(infile);
		return -1;
	}

	for (count = 0; !feof(infile) && !ferror(infile); count++)
	{
		if (blocksize[1] != 256)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid header block size in map file " << filename << "!" << endl;
			return -1;
		}

		if (fread(&metric.hdr, sizeof(metric.hdr), 1, infile) != 1)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read header block in map file " << filename << "!" << endl;
			return -1;
		}

		/*if (fread(&blocksize, sizeof(uint32_t), 1, infile) != 1)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read header block in map file " << filename << "!" << endl;
			return -1;
		}*/

		if (metric.hdr.header_blocksize != 256)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid header block size in map file " << filename << "!" << endl;
			return -1;
		}

		/*if (fread(&blocksize, sizeof(uint32_t), 1, infile) != 1)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read data block in map file " << filename << "!" << endl;
			return -1;
		}*/

		if (metric.hdr.data_blocksize != metric.hdr.precision * metric.hdr.Npix)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid data block size in map file " << filename << "!" << endl;
			return -1;
		}

		if (metric.hdr.distance > min_dist) break;

		backtrack = ftell(infile) - (256 + 2 * sizeof(uint32_t));

		if (fseek(infile, metric.hdr.data_blocksize, SEEK_CUR))
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to skip data block in map file " << filename << "!" << endl;
			fclose(infile);
			return -1;
		}

		if (fread(blocksize, sizeof(uint32_t), 2, infile) != 2)
		{
			fclose(infile);
			return -1;
		}
	}

	if (feof(infile) || ferror(infile))
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": read error occured in map file " << filename << "!" << endl;
		fclose(infile);
		return -1;
	}

	if (count > 0)
	{
		if (fseek(infile, backtrack, SEEK_SET))
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to backtrack in map file " << filename << "!" << endl;
			fclose(infile);
			return -1;
		}

		if (fread(&metric.hdr, sizeof(metric.hdr), 1, infile) != 1)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read header block in map file " << filename << "!" << endl;
			return -1;
		}

		/*if (fread(&blocksize, sizeof(uint32_t), 1, infile) != 1)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read header block in map file " << filename << "!" << endl;
			return -1;
		}*/

		if (metric.hdr.header_blocksize != 256)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid header block size in map file " << filename << "!" << endl;
			return -1;
		}

		/*if (fread(&blocksize, sizeof(uint32_t), 1, infile) != 1)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read data block in map file " << filename << "!" << endl;
			return -1;
		}*/

		if (metric.hdr.data_blocksize != metric.hdr.precision * metric.hdr.Npix)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid data block size in map file " << filename << "!" << endl;
			return -1;
		}

		count--;
	}

	while (true)
	{
		if (field->healpix_data.find(count) == field->healpix_data.end()) // data not present
		{
			if (metric.hdr.Nside < 2)
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid Nside = " << metric.hdr.Nside << " in map file " << filename << "!" << endl;
				fclose(infile);
				return -1;
			}

			if (metric.hdr.distance < 0)
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid distance = " << metric.hdr.distance << " in map file " << filename << "!" << endl;
				fclose(infile);
				return -1;
			}

			metric.pixel = (float *) malloc(metric.hdr.Npix * sizeof(float));

			if (metric.hdr.precision == 4)
			{
				if (fread(metric.pixel, sizeof(float), metric.hdr.Npix, infile) != metric.hdr.Npix)
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read data block in map file " << filename << "!" << endl;
					fclose(infile);
					free(metric.pixel);
					return -1;
				}
			}
			else if (metric.hdr.precision == 8)
			{
				double * dpix = (double *) malloc (metric.hdr.Npix * sizeof(double));
				if (fread(dpix, sizeof(double), metric.hdr.Npix, infile) != metric.hdr.Npix)
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read data block in map file " << filename << "!" << endl;
					fclose(infile);
					free(metric.pixel);
					return -1;
				}

				for (int i = 0; i < metric.hdr.Npix; i++) metric.pixel[i] = dpix[i];

				free(dpix);
			}
			else
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": precision " << metric.hdr.precision << " bytes not supported for map files!" << endl;
				free(metric.pixel);
			}

			field->healpix_data.insert(std::pair<int,metric_data>(count, metric));
		}
		else
		{
			if (fseek(infile, metric.hdr.data_blocksize, SEEK_CUR))
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to skip data block in map file " << filename << "!" << endl;
				fclose(infile);
				return -1;
			}
		}

		if (metric.hdr.distance > max_dist) break;

		/*if (fread(&blocksize, sizeof(uint32_t), 1, infile) != 1)
		{
			fclose(infile);
			return -1;
		}*/

		if (fread(blocksize, sizeof(uint32_t), 2, infile) != 2)
		{
			fclose(infile);
			return -2;
		}

		if (blocksize[1] != 256)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid header block size in map file " << filename << "!" << endl;
			return -1;
		}

		if (fread(&metric.hdr, sizeof(metric.hdr), 1, infile) != 1)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read header block in map file " << filename << "!" << endl;
			return -1;
		}

		/*if (fread(&blocksize, sizeof(uint32_t), 1, infile) != 1)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read header block in map file " << filename << "!" << endl;
			return -1;
		}*/

		if (metric.hdr.header_blocksize != 256)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid header block size in map file " << filename << "!" << endl;
			return -1;
		}

		/*if (fread(&blocksize, sizeof(uint32_t), 1, infile) != 1)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read data block in map file " << filename << "!" << endl;
			return -1;
		}*/

		if (metric.hdr.data_blocksize != metric.hdr.precision * metric.hdr.Npix)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid data block size in map file " << filename << "!" << endl;
			return -1;
		}

		count++;
	}

	fclose(infile);

	return count;
}

std::map<int,metric_data>::iterator get_iterator(metric_container * field, const double dist, const int dir = 1)
{
	std::map<int,metric_data>::iterator it;	

	for (it = field->healpix_data.begin(); it != field->healpix_data.end(); it++)
	{
		if (it->second.hdr.distance > dist)
		{
			if (dir < 0)
			{
				if (it != field->healpix_data.begin()) return --it;
				else return field->healpix_data.end();
			}
			else if (it != field->healpix_data.begin()) return it;
			else return field->healpix_data.end();
		}
	}

	return it;
}

void rotate_vector(double * input, double * output, double * coordsys)
{
	output[0] = input[0] * coordsys[0] * coordsys[2] + input[1] * coordsys[1] * coordsys[2] - input[2] * coordsys[3];
	output[1] = input[1] * coordsys[0] - input[0] * coordsys[1];
	output[2] = input[0] * coordsys[0] * coordsys[3] + input[1] * coordsys[1] * coordsys[3] + input[2] * coordsys[2];
}

void counterrotate_vector(double * input, double * output, double * coordsys)
{
	output[0] = input[0] * coordsys[0] * coordsys[2] - input[1] * coordsys[1] + input[2] * coordsys[0] * coordsys[3];
	output[1] = input[0] * coordsys[1] * coordsys[2] + input[1] * coordsys[0] + input[2] * coordsys[1] * coordsys[3];
	output[2] = input[2] * coordsys[2] - input[0] * coordsys[3];
}

bool interpolation(metric_container * field_prev, metric_container * field_next, const double weight, const double dist, const double dr, double * v, double * result)
{
	std::map<int,metric_data>::iterator it_prev;
	std::map<int,metric_data>::iterator it_next;
	double d;
	int64_t i;

	if(modf(dist / dr, &d) < 0.5)
	{
		it_prev = get_iterator(field_prev, dist, -1);
		it_next = get_iterator(field_next, dist, -1);
	}
	else
	{
		it_prev = get_iterator(field_prev, dist);
		it_next = get_iterator(field_next, dist);
	}

	if (it_prev == field_prev->healpix_data.end() || it_next == field_next->healpix_data.end()) return false;

	if (it_prev->second.hdr.distance != it_next->second.hdr.distance || it_prev->second.hdr.Nside != it_next->second.hdr.Nside) return false;

	if (dist == 0)
		i = 0;
	else
	{
		vec2pix_ring64(it_prev->second.hdr.Nside, v, &i);

		if (i >= it_prev->second.hdr.Npix || it_prev->second.pixel[i] < -1e30) return false;
		if (i >= it_next->second.hdr.Npix || it_next->second.pixel[i] < -1e30) return false;
	}

	(*result) = linear_interpolation(it_prev->second.pixel[i], it_next->second.pixel[i], weight);

	return true;
}

bool time_derivative(metric_container * field_prev, metric_container * field_next, const double dtau, const double dist, const double dr, double * v, double * result)
{
	std::map<int,metric_data>::iterator it_prev;
	std::map<int,metric_data>::iterator it_next;
	double d;
	int64_t i;

	if(modf(dist / dr, &d) < 0.5)
	{
		it_prev = get_iterator(field_prev, dist, -1);
		it_next = get_iterator(field_next, dist, -1);
	}
	else
	{
		it_prev = get_iterator(field_prev, dist);
		it_next = get_iterator(field_next, dist);
	}

	if (it_prev == field_prev->healpix_data.end() || it_next == field_next->healpix_data.end()) return false;

	if (it_prev->second.hdr.distance != it_next->second.hdr.distance || it_prev->second.hdr.Nside != it_next->second.hdr.Nside) return false;

	if (dist == 0)
		i = 0;
	else
	{
		vec2pix_ring64(it_prev->second.hdr.Nside, v, &i);

		if (i >= it_prev->second.hdr.Npix || it_prev->second.pixel[i] < -1e30) return false;
		if (i >= it_next->second.hdr.Npix || it_next->second.pixel[i] < -1e30) return false;
	}

	(*result) = (it_next->second.pixel[i] - it_prev->second.pixel[i]) / dtau;

	return true;
}

bool gradient(metric_container * field_prev, metric_container * field_next, const double weight, const double dist, const double dr, double * v, double * result, int debug)
{
	std::map<int,metric_data>::iterator inner_prev = get_iterator(field_prev, dist, -1);
	std::map<int,metric_data>::iterator inner_next = get_iterator(field_next, dist, -1);
	std::map<int,metric_data>::iterator outer_prev = get_iterator(field_prev, dist);
	std::map<int,metric_data>::iterator outer_next = get_iterator(field_next, dist);
	double d;
	int64_t i, j, k, l, ring, Npix;
	double grad[3];

	if (inner_prev == field_prev->healpix_data.end() || inner_next == field_next->healpix_data.end() || outer_prev == field_prev->healpix_data.end() || outer_next == field_next->healpix_data.end()) return false;

	if (inner_prev->second.hdr.distance != inner_next->second.hdr.distance || inner_prev->second.hdr.Nside != inner_next->second.hdr.Nside) return false;

	if (outer_prev->second.hdr.distance != outer_next->second.hdr.distance || outer_prev->second.hdr.Nside != outer_next->second.hdr.Nside) return false;

	if (v[0] == 0 && v[1] == 0)
	{
		if (v[2] < 0)
		{
			i = 12l * (int64_t) inner_prev->second.hdr.Nside * (int64_t) inner_prev->second.hdr.Nside - 1;
			j = 12l * (int64_t) outer_prev->second.hdr.Nside * (int64_t) outer_prev->second.hdr.Nside - 1;
			
		}
		else
		{
			i = 0;
			j = 0;
		}
	}
	else
	{
		vec2pix_ring64(inner_prev->second.hdr.Nside, v, &i);
		vec2pix_ring64(outer_prev->second.hdr.Nside, v, &j);
	}

	if (i >= inner_prev->second.hdr.Npix || inner_prev->second.pixel[i] < -1e30) return false;
	if (i >= inner_next->second.hdr.Npix || inner_next->second.pixel[i] < -1e30) return false;	

	if (j >= outer_prev->second.hdr.Npix || outer_prev->second.pixel[j] < -1e30) return false;
	if (j >= outer_next->second.hdr.Npix || outer_next->second.pixel[j] < -1e30) return false;

	grad[2] = linear_interpolation(outer_prev->second.pixel[j] - inner_prev->second.pixel[i], outer_next->second.pixel[j] - inner_next->second.pixel[i], weight) / dr;

	if(modf(dist / dr, &d) >= 0.5)
	{
		i = j;
		inner_prev = outer_prev;
		inner_next = outer_next;
	}

	if (inner_prev->second.hdr.distance == 0)
	{
		if (outer_prev->second.hdr.Nside < 2) return false;

		Npix = nside2npix64((int64_t) outer_prev->second.hdr.Nside);

		if (outer_prev->second.hdr.Npix < Npix || outer_next->second.hdr.Npix < Npix) return false;

		d = linear_interpolation(inner_prev->second.pixel[i], inner_next->second.pixel[i], weight);

		if (v[2] >= 0)
		{
			if (outer_prev->second.pixel[0] < -1e30 || outer_prev->second.pixel[1] < -1e30 || outer_prev->second.pixel[2] < -1e30 || outer_prev->second.pixel[3] < -1e30) return false;
			if (outer_next->second.pixel[0] < -1e30 || outer_next->second.pixel[1] < -1e30 || outer_next->second.pixel[2] < -1e30 || outer_next->second.pixel[3] < -1e30) return false;

			result[2] = linear_interpolation(outer_prev->second.pixel[0] + outer_prev->second.pixel[1] + outer_prev->second.pixel[2] + outer_prev->second.pixel[3], outer_next->second.pixel[0] + outer_next->second.pixel[1] + outer_next->second.pixel[2] + outer_next->second.pixel[3], weight);
			result[2] = (0.25 * result[2] - d) / dr;
		}
		else
		{
			if (outer_prev->second.pixel[Npix-1] < -1e30 || outer_prev->second.pixel[Npix-2] < -1e30 || outer_prev->second.pixel[Npix-3] < -1e30 || outer_prev->second.pixel[Npix-4] < -1e30) return false;
			if (outer_next->second.pixel[Npix-1] < -1e30 || outer_next->second.pixel[Npix-2] < -1e30 || outer_next->second.pixel[Npix-3] < -1e30 || outer_next->second.pixel[Npix-4] < -1e30) return false;
			result[2] = linear_interpolation(outer_prev->second.pixel[Npix-1] + outer_prev->second.pixel[Npix-2] + outer_prev->second.pixel[Npix-3] + outer_prev->second.pixel[Npix-4], outer_next->second.pixel[Npix-1] + outer_next->second.pixel[Npix-2] + outer_next->second.pixel[Npix-3] + outer_next->second.pixel[Npix-4], weight);
			result[2] = (d - 0.25 * result[2]) / dr;
		}

		if (v[0] >= 0)
		{
			j = (Npix/2) - 2 * outer_prev->second.hdr.Nside;
			k = (Npix/2) + 2 * outer_prev->second.hdr.Nside - 1;
			l = j - 4 * outer_prev->second.hdr.Nside;

			if (outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[k+1] < -1e30 || outer_prev->second.pixel[l] < -1e30) return false;
			if (outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[k+1] < -1e30 || outer_next->second.pixel[l] < -1e30) return false;

			result[0] = linear_interpolation(outer_prev->second.pixel[j] + outer_prev->second.pixel[k] + outer_prev->second.pixel[k+1] + outer_prev->second.pixel[l], outer_next->second.pixel[j] + outer_next->second.pixel[k] + outer_next->second.pixel[k+1] + outer_next->second.pixel[l], weight);
			result[0] = (0.25 * result[0] - d) / dr;
		}
		else
		{
			j = Npix/2;
			k = j - 4 * outer_prev->second.hdr.Nside;
			l = j + 4 * outer_prev->second.hdr.Nside;

			if (outer_prev->second.pixel[j-1] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[l] < -1e30) return false;
			if (outer_next->second.pixel[j-1] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[l] < -1e30) return false;

			result[0] = linear_interpolation(outer_prev->second.pixel[j-1] + outer_prev->second.pixel[j] + outer_prev->second.pixel[k] + outer_prev->second.pixel[l], outer_next->second.pixel[j-1] + outer_next->second.pixel[j] + outer_next->second.pixel[k] + outer_next->second.pixel[l], weight);
			result[0] = (d - 0.25 * result[0]) / dr;
		}

		if (v[1] >= 0)
		{
			j = Npix/2 - outer_prev->second.hdr.Nside;
			k = j - 4 * outer_prev->second.hdr.Nside;
			l = j + 4 * outer_prev->second.hdr.Nside;

			if (outer_prev->second.pixel[j-1] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[l] < -1e30) return false;
			if (outer_next->second.pixel[j-1] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[l] < -1e30) return false;

			result[1] = linear_interpolation(outer_prev->second.pixel[j-1] + outer_prev->second.pixel[j] + outer_prev->second.pixel[k] + outer_prev->second.pixel[l], outer_next->second.pixel[j-1] + outer_next->second.pixel[j] + outer_next->second.pixel[k] + outer_next->second.pixel[l], weight);
			result[1] = (0.25 * result[1] - d) / dr;
		}
		else
		{
			j = Npix/2 + outer_prev->second.hdr.Nside;
			k = j - 4 * outer_prev->second.hdr.Nside;
			l = j + 4 * outer_prev->second.hdr.Nside;

			if (outer_prev->second.pixel[j-1] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[l] < -1e30) return false;
			if (outer_next->second.pixel[j-1] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[l] < -1e30) return false;

			result[1] = linear_interpolation(outer_prev->second.pixel[j-1] + outer_prev->second.pixel[j] + outer_prev->second.pixel[k] + outer_prev->second.pixel[l], outer_next->second.pixel[j-1] + outer_next->second.pixel[j] + outer_next->second.pixel[k] + outer_next->second.pixel[l], weight);
			result[1] = (d - 0.25 * result[1]) / dr;
		}

		//return true;
	}
	else if (i < 4) // at north pole
	{
#ifdef DEBUG
		if (debug) cout << " at north pole ";
#endif

		d = 1. - 1. / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside);
		if (v[2] > d)
		{
#ifdef DEBUG
			if (debug) cout << "(innermost part)" << endl;
#endif
			if (inner_prev->second.hdr.Npix < 4 || inner_prev->second.pixel[0] < -1e30 || inner_prev->second.pixel[1] < -1e30 || inner_prev->second.pixel[2] < -1e30 || inner_prev->second.pixel[3] < -1e30) return false;
			if (inner_next->second.hdr.Npix < 4 || inner_next->second.pixel[0] < -1e30 || inner_next->second.pixel[1] < -1e30 || inner_next->second.pixel[2] < -1e30 || inner_next->second.pixel[3] < -1e30) return false;

			//d = sqrt(8.) * acos(d) * inner_prev->second.hdr.distance;
			d = sqrt(8. - 8. * d * d) * inner_prev->second.hdr.distance;
			result[0] = linear_interpolation(inner_prev->second.pixel[0] - inner_prev->second.pixel[2], inner_next->second.pixel[0] - inner_next->second.pixel[2], weight) / d;
			result[1] = linear_interpolation(inner_prev->second.pixel[1] - inner_prev->second.pixel[3], inner_next->second.pixel[1] - inner_next->second.pixel[3], weight) / d;

			d = result[0] - result[1]; // rotate
			result[1] += result[0];
			result[0] = d;
			result[2] = grad[2];

			//return true;
		}
		else if (inner_prev->second.hdr.Nside == 2)
		{
#ifdef DEBUG
			if (debug) cout << "(outer part, Nside=2)" << endl;
#endif
			if (inner_prev->second.hdr.Npix < 20 || inner_prev->second.pixel[2*i+4] < -1e30 || inner_prev->second.pixel[2*i+5] < -1e30 || inner_prev->second.pixel[2*i+13] < -1e30) return false;
			if (inner_next->second.hdr.Npix < 20 || inner_next->second.pixel[2*i+4] < -1e30 || inner_next->second.pixel[2*i+5] < -1e30 || inner_next->second.pixel[2*i+13] < -1e30) return false;

			grad[0] = linear_interpolation(inner_prev->second.pixel[2*i+13] - inner_prev->second.pixel[i], inner_next->second.pixel[2*i+13] - inner_next->second.pixel[i], weight) /
 inner_prev->second.hdr.distance / 0.797054997187545; // 2. * sin(0.5 * (acos(1./3.) - acos(11./12.)))
// 0.819821555018427; // (acos(1./3.) - acos(11./12.))

			grad[1] = linear_interpolation(inner_prev->second.pixel[2*i+5] - inner_prev->second.pixel[2*i+4], inner_next->second.pixel[2*i+5] - inner_next->second.pixel[2*i+4], weight) / inner_prev->second.hdr.distance / 0.570470779087523; // (sin(acos(2./3.)) * 2. * sin(M_PI / 8.)
// 0.585401227586727; // (sin(acos(2./3.)) * MP_I / 4.)

			// rotate
			result[0] = ((i % 2) ? (-0.52704627669473 * grad[2] - 0.471404520791031 * grad[0]) : (0.52704627669473 * grad[2] + 0.471404520791031 * grad[0])) - 0.707106781186547 * grad[1];
			result[1] = 0.52704627669473 * grad[2] + 0.471404520791031 * grad[0] + ((i % 2) ? -0.707106781186547 : 0.707106781186547) * grad[1];
		
			if (i > 1)
			{
				result[0] *= -1;
				result[1] *= -1;
			}

			result[2] = (grad[2] / 1.5) - 0.74535599249993 * grad[0];

			//return true;
		}
		else // Nside > 2
		{
#ifdef DEBUG
			if (debug) cout << "(outer part, Nside>2)" << endl;
#endif
			if (inner_prev->second.hdr.Npix < 24 || inner_prev->second.pixel[2*i+4] < -1e30 || inner_prev->second.pixel[2*i+5] < -1e30 || inner_prev->second.pixel[3*i+13] < -1e30) return false;
			if (inner_next->second.hdr.Npix < 24 || inner_next->second.pixel[2*i+4] < -1e30 || inner_next->second.pixel[2*i+5] < -1e30 || inner_next->second.pixel[3*i+13] < -1e30) return false;

			grad[0] = linear_interpolation(inner_prev->second.pixel[3*i+13] - inner_prev->second.pixel[i], inner_next->second.pixel[3*i+13] - inner_next->second.pixel[i], weight) / inner_prev->second.hdr.distance / (sqrt(6.66666666666666 - sqrt(16. + (4. / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside - 26.6666666666666) / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside) - 2. / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside) / (double) inner_prev->second.hdr.Nside);
// (acos(1. - 3./(inner_prev->second.hdr.Nside * inner_prev->second.hdr.Nside)) - acos(1. - 1./(3. * inner_prev->second.hdr.Nside * inner_prev->second.hdr.Nside)));

			d = sqrt(6. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside - 4.); // sin(theta)

			grad[1] = linear_interpolation(inner_prev->second.pixel[2*i+5] - inner_prev->second.pixel[2*i+4], inner_next->second.pixel[2*i+5] - inner_next->second.pixel[2*i+4], weight) / inner_prev->second.hdr.distance / (d * 0.5102445764867864 / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside);

			//rotate

			result[0] = ((i % 2) ? -0.471404520791032 : 0.471404520791032) * (d * grad[2] + (1.5 * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside - 2.) * grad[0]) / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside - 0.707106781186547 * grad[1];
			result[1] = 0.471404520791032 * (d * grad[2] + (1.5 * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside - 2.) * grad[0]) / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside + ((i % 2) ? -0.707106781186547 : 0.707106781186547) * grad[1];

			if (i > 1)
			{
				result[0] *= -1;
				result[1] *= -1;
			}

			result[2] = ((1.5 * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside - 2.) * grad[2] - d * grad[0]) / 1.5 / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside;

			//return true;
		}
	}
	else if (i < (int64_t) inner_prev->second.hdr.Nside * ((int64_t) inner_prev->second.hdr.Nside - 1) * 2l) // in north polar cap
	{
#ifdef DEBUG
		if (debug) cout << " in north polar cap" << endl;
#endif

		ring = (1 + (int64_t) sqrt(1.5 + 2 * (double) i)) / 2;
			
		if ((int64_t) inner_prev->second.hdr.Npix < 2 * (ring+1) * (ring+2)) return false;

		j = i - 2 * ring * (ring-1);

		k = (j * (ring+1)) / ring;

		d = (1.5 + k) - ((double) (ring+1) * (0.5 + j) / (double) ring);

		k += 2 * ring * (ring+1);
		l = k + 1;

		if (inner_prev->second.pixel[k] < -1e30 || inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[k] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

		grad[0] = linear_interpolation(d * inner_prev->second.pixel[k] + (1.-d) * inner_prev->second.pixel[l], d * inner_next->second.pixel[k] + (1.-d) * inner_next->second.pixel[l], weight);

		if (j % ring == 0)
			k = ((4 + j / ring) * (ring-1) - 1) % (4 * (ring-1));
		else
			k = (j * (ring-1)) / ring;

		d = (1.5 + k) - ((double) (ring-1) * (0.5 + j) / (double) ring);

		if (j == 0) d -= 4 * (ring-1);

		l = ((k+1) % (4 * (ring-1))) + 2 * (ring-1) * (ring-2);
		k += 2 * (ring-1) * (ring-2);

		if (inner_prev->second.pixel[k] < -1e30 || inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[k] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

		grad[0] -= linear_interpolation(d * inner_prev->second.pixel[k] + (1.-d) * inner_prev->second.pixel[l], d * inner_next->second.pixel[k] + (1.-d) * inner_next->second.pixel[l], weight);

		grad[0] /= sqrt(12. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside * (ring*ring+1) - 2. * (ring*ring-1) * (ring*ring - 1 + sqrt(((2.44948974278318* (double) inner_prev->second.hdr.Nside + 1.)*(2.44948974278318* (double) inner_prev->second.hdr.Nside + 1.) - ring*ring)*((2.44948974278318* (double) inner_prev->second.hdr.Nside - 1.)*(2.44948974278318* (double) inner_prev->second.hdr.Nside - 1.) - ring*ring)))) * inner_prev->second.hdr.distance / 3. / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside;
// (acos(1. - (ring+1) * (ring+1) / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside)) - acos(1. - (ring-1) * (ring-1) / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside))) * inner_prev->second.hdr.distance;

		//d = sin(acos(1. - ring * ring / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside))) * inner_prev->second.hdr.distance * M_PI / (2 * ring);
		d = ring * sqrt(6. - (double) (ring*ring) / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside) / 3. / (double) inner_prev->second.hdr.Nside; // sin(theta)

		result[0] = d * grad[2] + (1. - (ring*ring) / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside)) * grad[0];

		if (4*ring * v[3] > (2*j + 1) * M_PI)
		{
			k = 2 * ring * (ring-1) + (j+1) % (4*ring);

			if (inner_prev->second.pixel[k] < -1e30 || inner_next->second.pixel[k] < -1e30) return false;

			grad[1] = tan(M_PI / (4 * ring)) * (linear_interpolation(inner_prev->second.pixel[k] - inner_prev->second.pixel[i], inner_next->second.pixel[k] - inner_next->second.pixel[i], weight) / d / inner_prev->second.hdr.distance / (1. - cos(M_PI / (2 * ring))) + result[0]);
		}
		else
		{
			k = 2 * ring * (ring-1) + (4*ring + j-1) % (4*ring);

			if (inner_prev->second.pixel[k] < -1e30 || inner_next->second.pixel[k] < -1e30) return false;

			grad[1] = tan(M_PI / (4 * ring)) * (linear_interpolation(inner_prev->second.pixel[i] - inner_prev->second.pixel[k], inner_next->second.pixel[i] - inner_next->second.pixel[k], weight) / d / inner_prev->second.hdr.distance / (1. - cos(M_PI / (2 * ring))) - result[0]);
		}

		// rotate
		result[1] = sin(M_PI * (j + 0.5) / (2 * ring));
		result[2] = cos(M_PI * (j + 0.5) / (2 * ring));

		result[0] = result[0] * result[2] - result[1] * grad[1];
		result[1] = (d * grad[2] + (1. - (ring*ring) / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside)) * grad[0]) * result[1] + result[2] * grad[1];
		result[2] = (1. - (ring*ring) / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside)) * grad[2] - d * grad[0];

		//return true;
	}
	else if (i < (int64_t) inner_prev->second.hdr.Nside * ((int64_t) inner_prev->second.hdr.Nside + 1) * 2) // on northern circle
	{
#ifdef DEBUG
		if (debug) cout << " on northern circle" << endl;
#endif

		ring = (int64_t) inner_prev->second.hdr.Nside;

		if ((int64_t) inner_prev->second.hdr.Npix < 2 * ring * (ring+3)) return false;

		j = i - 2 * ring * (ring-1);

		k = i + 4 * ring;
		l = (j == ring-1) ? i+1 : k+1;

		if (inner_prev->second.pixel[k] < -1e30 || inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[k] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

		grad[0] = 0.5 * linear_interpolation(inner_prev->second.pixel[k] + inner_prev->second.pixel[l], inner_next->second.pixel[k] + inner_next->second.pixel[l], weight);

		if (j % ring == 0)
			k = ((4 + j / ring) * (ring-1) - 1) % (4 * (ring-1));
		else
			k = (j * (ring-1)) / ring;

		d = (1.5 + k) - ((double) (ring-1) * (0.5 + j) / (double) ring);

		if (j == 0) d -= 4 * (ring-1);

		l = ((k+1) % (4 * (ring-1))) + 2 * (ring-1) * (ring-2);
		k += 2 * (ring-1) * (ring-2);

		if (inner_prev->second.pixel[k] < -1e30 || inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[k] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

		grad[0] -= linear_interpolation(d * inner_prev->second.pixel[k] + (1.-d) * inner_prev->second.pixel[l], d * inner_next->second.pixel[k] + (1.-d) * inner_next->second.pixel[l], weight);

		grad[0] /= sqrt((5*ring*ring + 6 - 2./ring - (ring-1) * sqrt((ring+2)*(5*ring-2) * (double) (ring*(5*ring-2)-1)) / ring) / ring) * inner_prev->second.hdr.distance / 3. / ring;
// (acos((ring-1) / (1.5 * ring)) - acos(1. - (ring-1) * (ring-1) / (3. * ring * ring))) * inner_prev->second.hdr.distance;

		d = 0.74535599249993 * inner_prev->second.hdr.distance * (1. - cos(0.5 * M_PI / ring));
// 1.17080245517345 * inner_prev->second.hdr.distance / ring; // sin(acos(2./3.)) * M_PI / 2.;

		if (4*ring * v[3] > (2*j + 1) * M_PI)
		{
			k = 2 * ring * (ring-1) + (j+1) % (4*ring);

			if (inner_prev->second.pixel[k] < -1e30 || inner_next->second.pixel[k] < -1e30) return false;

			grad[1] = tan(M_PI / (4 * ring)) * (linear_interpolation(inner_prev->second.pixel[k] - inner_prev->second.pixel[i], inner_next->second.pixel[k] - inner_next->second.pixel[i], weight) / d + 0.74535599249993 * grad[2] + grad[0] / 1.5);
		}
		else
		{
			k = 2 * ring * (ring-1) + (4*ring + j-1) % (4*ring);

			if (inner_prev->second.pixel[k] < -1e30 || inner_next->second.pixel[k] < -1e30) return false;

			grad[1] = tan(M_PI / (4 * ring)) * (linear_interpolation(inner_prev->second.pixel[i] - inner_prev->second.pixel[k], inner_next->second.pixel[i] - inner_next->second.pixel[k], weight) / d - 0.74535599249993 * grad[2] - grad[0] / 1.5);
		}

		// rotate
		result[1] = sin(M_PI * (j + 0.5) / (2 * ring));
		result[2] = cos(M_PI * (j + 0.5) / (2 * ring));

		result[0] = (0.74535599249993 * grad[2] + grad[0] / 1.5) * result[2] - result[1] * grad[1];
		result[1] = (0.74535599249993 * grad[2] + grad[0] / 1.5) * result[1] + result[2] * grad[1];
		result[2] = grad[2] / 1.5 - 0.74535599249993 * grad[0];

		//return true;
	}
	else if (i < (int64_t) inner_prev->second.hdr.Nside * (10 * (int64_t) inner_prev->second.hdr.Nside - 2)) // in equatorial region
	{
#ifdef DEBUG
		if (debug) cout << " in equatorial region" << endl;
#endif

		j = i - 2 * (int64_t) inner_prev->second.hdr.Nside * ((int64_t) inner_prev->second.hdr.Nside + 1);
		ring = (int64_t) inner_prev->second.hdr.Nside + 1 + j / (4 * (int64_t) inner_prev->second.hdr.Nside);		

		j %= 4 * (int64_t) inner_prev->second.hdr.Nside;

		if ((int64_t) inner_prev->second.hdr.Npix < i + 4 * (int64_t) inner_prev->second.hdr.Nside + 1) return false;

		d = (2*j + 1 - ring%2) * M_PI;

		if (4*ring * v[3] > d)
		{
			d += M_PI;
			if (j+1 < 4 * (int64_t) inner_prev->second.hdr.Nside)
			{
				k = i - 4 * (int64_t) inner_prev->second.hdr.Nside + 1 - ring%2;
				j = i + 1;
			}
			else
			{
				k = i - ((ring%2) ? 4 * (int64_t) inner_prev->second.hdr.Nside : 8 * (int64_t) inner_prev->second.hdr.Nside - 1);
				j = i + 1 - 4 * (int64_t) inner_prev->second.hdr.Nside;
			}
			if (inner_prev->second.pixel[j] < -1e30 || inner_next->second.pixel[j] < -1e30) return false;
		}
		else
		{
			d -= M_PI;
			if (j > 0)
			{
				k = i - 4 * (int64_t) inner_prev->second.hdr.Nside - ring%2;
				j = i--;
			}
			else
			{
				k = i - ((ring%2) ? 1 : 4 * (int64_t) inner_prev->second.hdr.Nside);
				j = i;
				i += 4 * (int64_t) inner_prev->second.hdr.Nside - 1;
			}
			if (inner_prev->second.pixel[i] < -1e30 || inner_next->second.pixel[i] < -1e30) return false;
		}
		if (inner_prev->second.pixel[k] < -1e30 || inner_next->second.pixel[k] < -1e30) return false;

		l = k + 8 * (int64_t) inner_prev->second.hdr.Nside;

		if (inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

		grad[0] = linear_interpolation(inner_prev->second.pixel[l] - inner_prev->second.pixel[k], inner_next->second.pixel[l] - inner_next->second.pixel[k], weight) / inner_prev->second.hdr.distance / (sqrt(8. + 2. * (2 * ring - inner_prev->second.hdr.Nside) * (7 * inner_prev->second.hdr.Nside - 2 * ring) - 2. * sqrt((double) (2 * ring - 2 - inner_prev->second.hdr.Nside) * (double) (2 * ring + 2 - inner_prev->second.hdr.Nside) * (double) (7 * inner_prev->second.hdr.Nside + 2 - 2 * ring) * (double) (7 * inner_prev->second.hdr.Nside - 2 - 2 * ring))) / 3. / (double) inner_prev->second.hdr.Nside);
// (acos((2 * (int64_t) inner_prev->second.hdr.Nside - ring - 1) / (1.5 * inner_prev->second.hdr.Nside)) - acos((2 * (int64_t) inner_prev->second.hdr.Nside - ring + 1) / (1.5 * inner_prev->second.hdr.Nside)));

		result[1] = sin(d / (4 * ring));
		result[2] = cos(d / (4 * ring));

		d = sqrt(1. - (double) (2 * inner_prev->second.hdr.Nside - ring) * (double) (2 * inner_prev->second.hdr.Nside - ring) / (2.25 * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside)); // sin(theta)

		grad[1] = linear_interpolation(inner_prev->second.pixel[j] - inner_prev->second.pixel[i], inner_next->second.pixel[j] - inner_next->second.pixel[i], weight) / inner_prev->second.hdr.distance / (2. * d * sin(0.25 * M_PI / (double) inner_prev->second.hdr.Nside));
// (sin(acos((2 * (int64_t) inner_prev->second.hdr.Nside - ring) / (1.5 * inner_prev->second.hdr.Nside))) * M_PI / (2 * (int64_t) inner_prev->second.hdr.Nside));

		// rotate

		result[0] = (d * grad[2] + ((2 * inner_prev->second.hdr.Nside - ring) / 1.5 / (double) inner_prev->second.hdr.Nside) * grad[0]) * result[2] - result[1] * grad[1];
		result[1] = (d * grad[2] + ((2 * inner_prev->second.hdr.Nside - ring) / 1.5 / (double) inner_prev->second.hdr.Nside) * grad[0]) * result[1] + result[2] * grad[1];
		result[2] = ((2 * inner_prev->second.hdr.Nside - ring) / 1.5 / (double) inner_prev->second.hdr.Nside) * grad[2] - d * grad[0];

		//return true;
	}
	else if (i < (int64_t) inner_prev->second.hdr.Nside * (10 * (int64_t) inner_prev->second.hdr.Nside + 2)) // on southern circle
	{
#ifdef DEBUG
		if (debug) cout << " on southern circle" << endl;
#endif

		ring = (int64_t) inner_prev->second.hdr.Nside; // counted from south pole

		if ((int64_t) inner_prev->second.hdr.Npix < ring * (10 * ring + 6) - 4) return false;

		j = i - ring * (10 * ring - 2);

		if (j % ring == 0)
			k = ((4 + j / ring) * (ring-1) - 1) % (4 * (ring-1));
		else
			k = (j * (ring-1)) / ring;

		d = (1.5 + k) - ((double) (ring-1) * (0.5 + j) / (double) ring);

		if (j == 0) d -= 4 * (ring-1);

		l = ((k+1) % (4 * (ring-1))) + ring * (10 * ring + 2);
		k += ring * (10 * ring + 2);

		if (inner_prev->second.pixel[k] < -1e30 || inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[k] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

		grad[0] = linear_interpolation(d * inner_prev->second.pixel[k] + (1.-d) * inner_prev->second.pixel[l], d * inner_next->second.pixel[k] + (1.-d) * inner_next->second.pixel[l], weight);

		k = i - 4 * ring;
		l = (j == ring-1) ? k+1 - 4 * ring : k+1;

		if (inner_prev->second.pixel[k] < -1e30 || inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[k] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

		grad[0] -= 0.5 * linear_interpolation(inner_prev->second.pixel[k] + inner_prev->second.pixel[l], inner_next->second.pixel[k] + inner_next->second.pixel[l], weight);

		grad[0] /= sqrt((5*ring*ring + 6 - 2./ring - (ring-1) * sqrt((ring+2)*(5*ring-2) * (double) (ring*(5*ring-2)-1)) / ring) / ring) * inner_prev->second.hdr.distance / 3. / ring;
// (acos(-1. + (ring-1) * (ring-1) / (3. * ring * ring)) - acos((1-ring) / (1.5 * ring))) * inner_prev->second.hdr.distance;

		d = 0.74535599249993 * inner_prev->second.hdr.distance * (1. - cos(0.5 * M_PI / ring));
// 1.17080245517345 * inner_prev->second.hdr.distance / ring; // sin(acos(-2./3.)) * M_PI / 2.;

		if (4*ring * v[3] > (2*j + 1) * M_PI)
		{
			k = ring * (10 * ring - 2) + (j+1) % (4*ring);

			if (inner_prev->second.pixel[k] < -1e30 || inner_next->second.pixel[k] < -1e30) return false;

			grad[1] = tan(M_PI / (4 * ring)) * (linear_interpolation(inner_prev->second.pixel[k] - inner_prev->second.pixel[i], inner_next->second.pixel[k] - inner_next->second.pixel[i], weight) / d + 0.74535599249993 * grad[2] - grad[0] / 1.5);
		}
		else
		{
			k = ring * (10 * ring - 2) + (4*ring + j-1) % (4*ring);

			if (inner_prev->second.pixel[k] < -1e30 || inner_next->second.pixel[k] < -1e30) return false;

			grad[1] = tan(M_PI / (4 * ring)) * (linear_interpolation(inner_prev->second.pixel[i] - inner_prev->second.pixel[k], inner_next->second.pixel[i] - inner_next->second.pixel[k], weight) / d - 0.74535599249993 * grad[2] + grad[0] / 1.5);
		}

		// rotate
		result[1] = sin(M_PI * (j + 0.5) / (2 * ring));
		result[2] = cos(M_PI * (j + 0.5) / (2 * ring));

		result[0] = (0.74535599249993 * grad[2] - grad[0] / 1.5) * result[2] - result[1] * grad[1];
		result[1] = (0.74535599249993 * grad[2] - grad[0] / 1.5) * result[1] + result[2] * grad[1];
		result[2] = -grad[2] / 1.5 - 0.74535599249993 * grad[0];

		//return true;
	}
	else if (i < 12 * (int64_t) inner_prev->second.hdr.Nside * (int64_t) inner_prev->second.hdr.Nside - 4) // in south polar cap
	{
#ifdef DEBUG
		if (debug) cout << " in south polar cap" << endl;
#endif

		Npix = nside2npix64((int64_t) inner_prev->second.hdr.Nside);
		ring = (1 + (int64_t) sqrt(2 * (Npix-i) - 0.5)) / 2; // counted from south pole

		if ((int64_t) inner_prev->second.hdr.Npix < Npix - 2 * ring * (ring-1)) return false;

		j = i + 2 * ring * (ring+1) - Npix;

		if (j % ring == 0)
			k = ((4 + j / ring) * (ring-1) - 1) % (4 * (ring-1));
		else
			k = (j * (ring-1)) / ring;

		d = (1.5 + k) - ((double) (ring-1) * (0.5 + j) / (double) ring);

		if (j == 0) d -= 4 * (ring-1);

		l = ((k+1) % (4 * (ring-1))) + Npix - 2 * ring * (ring-1);
		k += Npix - 2 * ring * (ring-1);

		if (inner_prev->second.pixel[k] < -1e30 || inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[k] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

		grad[0] = linear_interpolation(d * inner_prev->second.pixel[k] + (1.-d) * inner_prev->second.pixel[l], d * inner_next->second.pixel[k] + (1.-d) * inner_next->second.pixel[l], weight);

		k = (j * (ring+1)) / ring;

		d = (1.5 + k) - ((double) (ring+1) * (0.5 + j) / (double) ring);

		k += Npix - 2 * (ring+1) * (ring+2);
		l = k + 1;

		if (inner_prev->second.pixel[k] < -1e30 || inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[k] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

		grad[0] -= linear_interpolation(d * inner_prev->second.pixel[k] + (1.-d) * inner_prev->second.pixel[l], d * inner_next->second.pixel[k] + (1.-d) * inner_next->second.pixel[l], weight);

		grad[0] /= sqrt(12. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside * (ring*ring+1) - 2. * (ring*ring-1) * (ring*ring - 1 + sqrt(((2.44948974278318* (double) inner_prev->second.hdr.Nside + 1.)*(2.44948974278318* (double) inner_prev->second.hdr.Nside + 1.) - ring*ring)*((2.44948974278318* (double) inner_prev->second.hdr.Nside - 1.)*(2.44948974278318* (double) inner_prev->second.hdr.Nside - 1.) - ring*ring)))) * inner_prev->second.hdr.distance / 3. / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside;
// (acos(1. - (ring+1) * (ring+1) / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside)) - acos(1. - (ring-1) * (ring-1) / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside))) * inner_prev->second.hdr.distance;

		//d = sin(acos(1. - ring * ring / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside))) * inner_prev->second.hdr.distance * M_PI / (2 * ring);
		//d = ring * sqrt(6. - (double) (ring*ring) / (double) inner_prev->second.hdr.Nside / inner_prev->second.hdr.Nside) * inner_prev->second.hdr.distance * sin(M_PI / (4 * ring)) / 1.5 / inner_prev->second.hdr.Nside;
		d = ring * sqrt(6. - (double) (ring*ring) / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside) / 3. / (double) inner_prev->second.hdr.Nside; // sin(theta)

		result[0] = d * grad[2] - (1. - (ring*ring) / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside)) * grad[0];

		if (4*ring * v[3] > (2*j + 1) * M_PI)
		{
			k = Npix - 2 * ring * (ring+1) + (j+1) % (4*ring);

			if (inner_prev->second.pixel[k] < -1e30 || inner_next->second.pixel[k] < -1e30) return false;

			grad[1] = tan(M_PI / (4 * ring)) * (linear_interpolation(inner_prev->second.pixel[k] - inner_prev->second.pixel[i], inner_next->second.pixel[k] - inner_next->second.pixel[i], weight) / d / inner_prev->second.hdr.distance / (1. - cos(M_PI / (2 * ring))) + result[0]);
		}
		else
		{
			k = Npix - 2 * ring * (ring+1) + (4*ring + j-1) % (4*ring);

			if (inner_prev->second.pixel[k] < -1e30 || inner_next->second.pixel[k] < -1e30) return false;

			grad[1] = tan(M_PI / (4 * ring)) * (linear_interpolation(inner_prev->second.pixel[i] - inner_prev->second.pixel[k], inner_next->second.pixel[i] - inner_next->second.pixel[k], weight) / d / inner_prev->second.hdr.distance / (1. - cos(M_PI / (2 * ring))) - result[0]);
		}

		// rotate
		result[1] = sin(M_PI * (j + 0.5) / (2 * ring));
		result[2] = cos(M_PI * (j + 0.5) / (2 * ring));

		result[0] = result[0] * result[2] - result[1] * grad[1];
		result[1] = (d * grad[2] - (1. - (ring*ring) / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside)) * grad[0]) * result[1] + result[2] * grad[1];
		result[2] = (-1. + (ring*ring) / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside)) * grad[2] - d * grad[0];

		//return true;
	}
	else // at south pole
	{
#ifdef DEBUG
		if (debug) cout << " at south pole" << endl;
#endif

		Npix = nside2npix64((int64_t) inner_prev->second.hdr.Nside);
		if ((int64_t) inner_prev->second.hdr.Npix < Npix || (int64_t) inner_next->second.hdr.Npix < Npix) return false;
		d = -1. + 1. / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside);
		if (v[2] < d)
		{
			if (inner_prev->second.pixel[Npix-1] < -1e30 || inner_prev->second.pixel[Npix-2] < -1e30 || inner_prev->second.pixel[Npix-3] < -1e30 || inner_prev->second.pixel[Npix-4] < -1e30) return false;
			if (inner_next->second.pixel[Npix-1] < -1e30 || inner_next->second.pixel[Npix-2] < -1e30 || inner_next->second.pixel[Npix-3] < -1e30 || inner_next->second.pixel[Npix-4] < -1e30) return false;

			//d = sqrt(8.) * acos(d) * inner_prev->second.hdr.distance;
			d = sqrt(8. - 8. * d * d) * inner_prev->second.hdr.distance;
			result[0] = linear_interpolation(inner_prev->second.pixel[Npix-4] - inner_prev->second.pixel[Npix-2], inner_next->second.pixel[Npix-4] - inner_next->second.pixel[Npix-2], weight) / d;
			result[1] = linear_interpolation(inner_prev->second.pixel[Npix-3] - inner_prev->second.pixel[Npix-1], inner_next->second.pixel[Npix-3] - inner_next->second.pixel[Npix-1], weight) / d;

			d = result[0] - result[1]; // rotate
			result[1] += result[0];
			result[0] = d;
			result[2] = -grad[2];

			//return true;
		}
		else if (inner_prev->second.hdr.Nside == 2)
		{
			j = i-44;

			if (inner_prev->second.pixel[2*j+36] < -1e30 || inner_prev->second.pixel[2*j+37] < -1e30 || inner_prev->second.pixel[2*j+29] < -1e30) return false;
			if (inner_next->second.pixel[2*j+36] < -1e30 || inner_next->second.pixel[2*j+37] < -1e30 || inner_next->second.pixel[2*j+29] < -1e30) return false;

			grad[0] = linear_interpolation(inner_prev->second.pixel[i] - inner_prev->second.pixel[2*j+29], inner_next->second.pixel[i] - inner_next->second.pixel[2*j+29], weight) /
 inner_prev->second.hdr.distance / 0.797054997187545; // 2. * sin(0.5 * (acos(1./3.) - acos(11./12.))) / 0.819821555018427; // (acos(1./3.) - acos(11./12.))

			grad[1] = linear_interpolation(inner_prev->second.pixel[2*j+37] - inner_prev->second.pixel[2*j+36], inner_next->second.pixel[2*j+37] - inner_next->second.pixel[2*j+36], weight) / inner_prev->second.hdr.distance / 0.570470779087523; // (sin(acos(2./3.)) * 2. * sin(M_PI / 8.) / 0.585401227586727; // (sin(acos(2./3.)) * MP_I / 4.)

			// rotate
			result[0] = ((i % 2) ? (-0.52704627669473 * grad[2] + 0.471404520791031 * grad[0]) : (0.52704627669473 * grad[2] - 0.471404520791031 * grad[0])) - 0.707106781186547 * grad[1];
			result[1] = 0.52704627669473 * grad[2] - 0.471404520791031 * grad[0] + ((i % 2) ? -0.707106781186547 : 0.707106781186547) * grad[1];
		
			if (i > 1)
			{
				result[0] *= -1;
				result[1] *= -1;
			}

			result[2] = -(grad[2] / 1.5) - 0.74535599249993 * grad[0];

			//return true;
		}
		else // Nside > 2
		{
			j = i + 4 - Npix;

			if (inner_prev->second.pixel[2*j+Npix-12] < -1e30 || inner_prev->second.pixel[2*j+Npix-11] < -1e30 || inner_prev->second.pixel[3*j+Npix-23] < -1e30) return false;
			if (inner_next->second.pixel[2*j+Npix-12] < -1e30 || inner_next->second.pixel[2*j+Npix-11] < -1e30 || inner_next->second.pixel[3*j+Npix-23] < -1e30) return false;

			grad[0] = linear_interpolation(inner_prev->second.pixel[i] - inner_prev->second.pixel[3*j+Npix-23], inner_next->second.pixel[i] - inner_next->second.pixel[3*j+Npix-23], weight) / inner_prev->second.hdr.distance / (sqrt(6.66666666666666 - sqrt(16. + (4. / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside - 26.6666666666666) / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside) - 2. / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside) / (double) inner_prev->second.hdr.Nside);
// (acos(1. - 3./(inner_prev->second.hdr.Nside * inner_prev->second.hdr.Nside)) - acos(1. - 1./(3. * inner_prev->second.hdr.Nside * inner_prev->second.hdr.Nside)));

			d = sqrt(6. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside - 4.); // sin(theta)

			grad[1] = linear_interpolation(inner_prev->second.pixel[2*j+Npix-11] - inner_prev->second.pixel[2*j+Npix-12], inner_next->second.pixel[2*j+Npix-11] - inner_next->second.pixel[2*j+Npix-12], weight) / inner_prev->second.hdr.distance / (d * 0.5102445764867864 / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside);
// (sin(acos(1. - 4. / (3. * inner_prev->second.hdr.Nside * inner_prev->second.hdr.Nside))) * M_PI / 4.);

			//rotate

			result[0] = ((i % 2) ? -0.471404520791032 : 0.471404520791032) * (d * grad[2] - (1.5 * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside - 2.) * grad[0]) / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside - 0.707106781186547 * grad[1];
			result[1] = 0.471404520791032 * (d * grad[2] - (1.5 * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside - 2.) * grad[0]) / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside + ((i % 2) ? -0.707106781186547 : 0.707106781186547) * grad[1];

			if (i > 1)
			{
				result[0] *= -1;
				result[1] *= -1;
			}

			result[2] = ((2. - 1.5 * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside) * grad[2] - d * grad[0]) / 1.5 / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside;

			//return true;
		}
	}

/*	d = sqrt(1. - v[2] * v[2]); // sin(theta)

	result[0] = v[0] * grad[2] + (v[0] * v[2] * grad[0] - v[1] * grad[1]) / d;
	result[1] = v[1] * grad[2] + (v[1] * v[2] * grad[0] + v[0] * grad[1]) / d;
	result[2] = v[2] * grad[2] - d * grad[0];*/

	return true;
}

bool time_derivative_gradient(metric_container * field_prev, metric_container * field_next, const double dtau, const double dist, const double dr, double * v, double * result, int debug)
{
	std::map<int,metric_data>::iterator inner_prev = get_iterator(field_prev, dist, -1);
	std::map<int,metric_data>::iterator inner_next = get_iterator(field_next, dist, -1);
	std::map<int,metric_data>::iterator outer_prev = get_iterator(field_prev, dist);
	std::map<int,metric_data>::iterator outer_next = get_iterator(field_next, dist);
	double d;
	int64_t i, j, k, l, ring, Npix;
	double grad[3];

	if (inner_prev == field_prev->healpix_data.end() || inner_next == field_next->healpix_data.end() || outer_prev == field_prev->healpix_data.end() || outer_next == field_next->healpix_data.end()) return false;

	if (inner_prev->second.hdr.distance != inner_next->second.hdr.distance || inner_prev->second.hdr.Nside != inner_next->second.hdr.Nside) return false;

	if (outer_prev->second.hdr.distance != outer_next->second.hdr.distance || outer_prev->second.hdr.Nside != outer_next->second.hdr.Nside) return false;

	if (v[0] == 0 && v[1] == 0)
	{
		if (v[2] < 0)
		{
			i = 12l * (int64_t) inner_prev->second.hdr.Nside * (int64_t) inner_prev->second.hdr.Nside - 1;
			j = 12l * (int64_t) outer_prev->second.hdr.Nside * (int64_t) outer_prev->second.hdr.Nside - 1;
			
		}
		else
		{
			i = 0;
			j = 0;
		}
	}
	else
	{
		vec2pix_ring64(inner_prev->second.hdr.Nside, v, &i);
		vec2pix_ring64(outer_prev->second.hdr.Nside, v, &j);
	}

	if (i >= inner_prev->second.hdr.Npix || inner_prev->second.pixel[i] < -1e30) return false;
	if (i >= inner_next->second.hdr.Npix || inner_next->second.pixel[i] < -1e30) return false;	

	if (j >= outer_prev->second.hdr.Npix || outer_prev->second.pixel[j] < -1e30) return false;
	if (j >= outer_next->second.hdr.Npix || outer_next->second.pixel[j] < -1e30) return false;

	grad[2] = ((outer_next->second.pixel[j] - inner_next->second.pixel[i]) - (outer_prev->second.pixel[j] - inner_prev->second.pixel[i])) / dtau / dr;

	if(modf(dist / dr, &d) >= 0.5)
	{
		i = j;
		inner_prev = outer_prev;
		inner_next = outer_next;
	}

	if (inner_prev->second.hdr.distance == 0)
	{
		if (outer_prev->second.hdr.Nside < 2) return false;

		Npix = nside2npix64((int64_t) outer_prev->second.hdr.Nside);

		if (outer_prev->second.hdr.Npix < Npix || outer_next->second.hdr.Npix < Npix) return false;

		d = (inner_next->second.pixel[i] - inner_prev->second.pixel[i]) / dtau;

		if (v[2] >= 0)
		{
			if (outer_prev->second.pixel[0] < -1e30 || outer_prev->second.pixel[1] < -1e30 || outer_prev->second.pixel[2] < -1e30 || outer_prev->second.pixel[3] < -1e30) return false;
			if (outer_next->second.pixel[0] < -1e30 || outer_next->second.pixel[1] < -1e30 || outer_next->second.pixel[2] < -1e30 || outer_next->second.pixel[3] < -1e30) return false;

			result[2] = ((outer_next->second.pixel[0] + outer_next->second.pixel[1] + outer_next->second.pixel[2] + outer_next->second.pixel[3]) - (outer_prev->second.pixel[0] + outer_prev->second.pixel[1] + outer_prev->second.pixel[2] + outer_prev->second.pixel[3])) / dtau;
			result[2] = (0.25 * result[2] - d) / dr;
		}
		else
		{
			if (outer_prev->second.pixel[Npix-1] < -1e30 || outer_prev->second.pixel[Npix-2] < -1e30 || outer_prev->second.pixel[Npix-3] < -1e30 || outer_prev->second.pixel[Npix-4] < -1e30) return false;
			if (outer_next->second.pixel[Npix-1] < -1e30 || outer_next->second.pixel[Npix-2] < -1e30 || outer_next->second.pixel[Npix-3] < -1e30 || outer_next->second.pixel[Npix-4] < -1e30) return false;
			result[2] = ((outer_next->second.pixel[Npix-1] + outer_next->second.pixel[Npix-2] + outer_next->second.pixel[Npix-3] + outer_next->second.pixel[Npix-4]) - (outer_prev->second.pixel[Npix-1] + outer_prev->second.pixel[Npix-2] + outer_prev->second.pixel[Npix-3] + outer_prev->second.pixel[Npix-4])) /dtau;
			result[2] = (d - 0.25 * result[2]) / dr;
		}

		if (v[0] >= 0)
		{
			j = (Npix/2) - 2 * outer_prev->second.hdr.Nside;
			k = (Npix/2) + 2 * outer_prev->second.hdr.Nside - 1;
			l = j - 4 * outer_prev->second.hdr.Nside;

			if (outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[k+1] < -1e30 || outer_prev->second.pixel[l] < -1e30) return false;
			if (outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[k+1] < -1e30 || outer_next->second.pixel[l] < -1e30) return false;

			result[0] = ((outer_next->second.pixel[j] + outer_next->second.pixel[k] + outer_next->second.pixel[k+1] + outer_next->second.pixel[l]) - (outer_prev->second.pixel[j] + outer_prev->second.pixel[k] + outer_prev->second.pixel[k+1] + outer_prev->second.pixel[l])) / dtau;
			result[0] = (0.25 * result[0] - d) / dr;
		}
		else
		{
			j = Npix/2;
			k = j - 4 * outer_prev->second.hdr.Nside;
			l = j + 4 * outer_prev->second.hdr.Nside;

			if (outer_prev->second.pixel[j-1] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[l] < -1e30) return false;
			if (outer_next->second.pixel[j-1] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[l] < -1e30) return false;

			result[0] = ((outer_next->second.pixel[j-1] + outer_next->second.pixel[j] + outer_next->second.pixel[k] + outer_next->second.pixel[l]) - (outer_prev->second.pixel[j-1] + outer_prev->second.pixel[j] + outer_prev->second.pixel[k] + outer_prev->second.pixel[l])) / dtau;
			result[0] = (d - 0.25 * result[0]) / dr;
		}

		if (v[1] >= 0)
		{
			j = Npix/2 - outer_prev->second.hdr.Nside;
			k = j - 4 * outer_prev->second.hdr.Nside;
			l = j + 4 * outer_prev->second.hdr.Nside;

			if (outer_prev->second.pixel[j-1] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[l] < -1e30) return false;
			if (outer_next->second.pixel[j-1] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[l] < -1e30) return false;

			result[1] = ((outer_next->second.pixel[j-1] + outer_next->second.pixel[j] + outer_next->second.pixel[k] + outer_next->second.pixel[l]) - (outer_prev->second.pixel[j-1] + outer_prev->second.pixel[j] + outer_prev->second.pixel[k] + outer_prev->second.pixel[l])) / dtau;
			result[1] = (0.25 * result[1] - d) / dr;
		}
		else
		{
			j = Npix/2 + outer_prev->second.hdr.Nside;
			k = j - 4 * outer_prev->second.hdr.Nside;
			l = j + 4 * outer_prev->second.hdr.Nside;

			if (outer_prev->second.pixel[j-1] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[l] < -1e30) return false;
			if (outer_next->second.pixel[j-1] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[l] < -1e30) return false;

			result[1] = ((outer_next->second.pixel[j-1] + outer_next->second.pixel[j] + outer_next->second.pixel[k] + outer_next->second.pixel[l]) - (outer_prev->second.pixel[j-1] + outer_prev->second.pixel[j] + outer_prev->second.pixel[k] + outer_prev->second.pixel[l])) / dtau;
			result[1] = (d - 0.25 * result[1]) / dr;
		}

		//return true;
	}
	else if (i < 4) // at north pole
	{
#ifdef DEBUG
		if (debug) cout << " at north pole ";
#endif

		d = 1. - 1. / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside);
		if (v[2] > d)
		{
#ifdef DEBUG
			if (debug) cout << "(innermost part)" << endl;
#endif
			if (inner_prev->second.hdr.Npix < 4 || inner_prev->second.pixel[0] < -1e30 || inner_prev->second.pixel[1] < -1e30 || inner_prev->second.pixel[2] < -1e30 || inner_prev->second.pixel[3] < -1e30) return false;
			if (inner_next->second.hdr.Npix < 4 || inner_next->second.pixel[0] < -1e30 || inner_next->second.pixel[1] < -1e30 || inner_next->second.pixel[2] < -1e30 || inner_next->second.pixel[3] < -1e30) return false;

			//d = sqrt(8.) * acos(d) * inner_prev->second.hdr.distance;
			d = sqrt(8. - 8. * d * d) * inner_prev->second.hdr.distance;
			result[0] = ((inner_next->second.pixel[0] - inner_next->second.pixel[2]) - (inner_prev->second.pixel[0] - inner_prev->second.pixel[2])) / dtau / d;
			result[1] = ((inner_next->second.pixel[1] - inner_next->second.pixel[3]) - (inner_prev->second.pixel[1] - inner_prev->second.pixel[3])) / dtau / d;

			d = result[0] - result[1]; // rotate
			result[1] += result[0];
			result[0] = d;
			result[2] = grad[2];

			//return true;
		}
		else if (inner_prev->second.hdr.Nside == 2)
		{
#ifdef DEBUG
			if (debug) cout << "(outer part, Nside=2)" << endl;
#endif
			if (inner_prev->second.hdr.Npix < 20 || inner_prev->second.pixel[2*i+4] < -1e30 || inner_prev->second.pixel[2*i+5] < -1e30 || inner_prev->second.pixel[2*i+13] < -1e30) return false;
			if (inner_next->second.hdr.Npix < 20 || inner_next->second.pixel[2*i+4] < -1e30 || inner_next->second.pixel[2*i+5] < -1e30 || inner_next->second.pixel[2*i+13] < -1e30) return false;

			grad[0] = ((inner_next->second.pixel[2*i+13] - inner_next->second.pixel[i]) - (inner_prev->second.pixel[2*i+13] - inner_prev->second.pixel[i])) / dtau /
 inner_prev->second.hdr.distance / 0.797054997187545; // 2. * sin(0.5 * (acos(1./3.) - acos(11./12.)))
// 0.819821555018427; // (acos(1./3.) - acos(11./12.))

			grad[1] = ((inner_next->second.pixel[2*i+5] - inner_next->second.pixel[2*i+4]) - (inner_prev->second.pixel[2*i+5] - inner_prev->second.pixel[2*i+4])) / dtau / inner_prev->second.hdr.distance / 0.570470779087523; // (sin(acos(2./3.)) * 2. * sin(M_PI / 8.)
// 0.585401227586727; // (sin(acos(2./3.)) * MP_I / 4.)

			// rotate
			result[0] = ((i % 2) ? (-0.52704627669473 * grad[2] - 0.471404520791031 * grad[0]) : (0.52704627669473 * grad[2] + 0.471404520791031 * grad[0])) - 0.707106781186547 * grad[1];
			result[1] = 0.52704627669473 * grad[2] + 0.471404520791031 * grad[0] + ((i % 2) ? -0.707106781186547 : 0.707106781186547) * grad[1];
		
			if (i > 1)
			{
				result[0] *= -1;
				result[1] *= -1;
			}

			result[2] = (grad[2] / 1.5) - 0.74535599249993 * grad[0];

			//return true;
		}
		else // Nside > 2
		{
#ifdef DEBUG
			if (debug) cout << "(outer part, Nside>2)" << endl;
#endif
			if (inner_prev->second.hdr.Npix < 24 || inner_prev->second.pixel[2*i+4] < -1e30 || inner_prev->second.pixel[2*i+5] < -1e30 || inner_prev->second.pixel[3*i+13] < -1e30) return false;
			if (inner_next->second.hdr.Npix < 24 || inner_next->second.pixel[2*i+4] < -1e30 || inner_next->second.pixel[2*i+5] < -1e30 || inner_next->second.pixel[3*i+13] < -1e30) return false;

			grad[0] = ((inner_next->second.pixel[3*i+13] - inner_next->second.pixel[i]) - (inner_prev->second.pixel[3*i+13] - inner_prev->second.pixel[i])) / dtau / inner_prev->second.hdr.distance / (sqrt(6.66666666666666 - sqrt(16. + (4. / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside - 26.6666666666666) / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside) - 2. / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside) / (double) inner_prev->second.hdr.Nside);
// (acos(1. - 3./(inner_prev->second.hdr.Nside * inner_prev->second.hdr.Nside)) - acos(1. - 1./(3. * inner_prev->second.hdr.Nside * inner_prev->second.hdr.Nside)));

			d = sqrt(6. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside - 4.); // sin(theta)

			grad[1] = ((inner_next->second.pixel[2*i+5] - inner_next->second.pixel[2*i+4]) - (inner_prev->second.pixel[2*i+5] - inner_prev->second.pixel[2*i+4])) / dtau / inner_prev->second.hdr.distance / (d * 0.5102445764867864 / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside);

			//rotate

			result[0] = ((i % 2) ? -0.471404520791032 : 0.471404520791032) * (d * grad[2] + (1.5 * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside - 2.) * grad[0]) / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside - 0.707106781186547 * grad[1];
			result[1] = 0.471404520791032 * (d * grad[2] + (1.5 * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside - 2.) * grad[0]) / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside + ((i % 2) ? -0.707106781186547 : 0.707106781186547) * grad[1];

			if (i > 1)
			{
				result[0] *= -1;
				result[1] *= -1;
			}

			result[2] = ((1.5 * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside - 2.) * grad[2] - d * grad[0]) / 1.5 / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside;

			//return true;
		}
	}
	else if (i < (int64_t) inner_prev->second.hdr.Nside * ((int64_t) inner_prev->second.hdr.Nside - 1) * 2l) // in north polar cap
	{
#ifdef DEBUG
		if (debug) cout << " in north polar cap" << endl;
#endif

		ring = (1 + (int64_t) sqrt(1.5 + 2 * (double) i)) / 2;
			
		if ((int64_t) inner_prev->second.hdr.Npix < 2 * (ring+1) * (ring+2)) return false;

		j = i - 2 * ring * (ring-1);

		k = (j * (ring+1)) / ring;

		d = (1.5 + k) - ((double) (ring+1) * (0.5 + j) / (double) ring);

		k += 2 * ring * (ring+1);
		l = k + 1;

		if (inner_prev->second.pixel[k] < -1e30 || inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[k] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

		grad[0] = ((d * inner_next->second.pixel[k] + (1.-d) * inner_next->second.pixel[l]) - (d * inner_prev->second.pixel[k] + (1.-d) * inner_prev->second.pixel[l])) / dtau;

		if (j % ring == 0)
			k = ((4 + j / ring) * (ring-1) - 1) % (4 * (ring-1));
		else
			k = (j * (ring-1)) / ring;

		d = (1.5 + k) - ((double) (ring-1) * (0.5 + j) / (double) ring);

		if (j == 0) d -= 4 * (ring-1);

		l = ((k+1) % (4 * (ring-1))) + 2 * (ring-1) * (ring-2);
		k += 2 * (ring-1) * (ring-2);

		if (inner_prev->second.pixel[k] < -1e30 || inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[k] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

		grad[0] -= ((d * inner_next->second.pixel[k] + (1.-d) * inner_next->second.pixel[l]) - (d * inner_prev->second.pixel[k] + (1.-d) * inner_prev->second.pixel[l])) / dtau;

		grad[0] /= sqrt(12. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside * (ring*ring+1) - 2. * (ring*ring-1) * (ring*ring - 1 + sqrt(((2.44948974278318* (double) inner_prev->second.hdr.Nside + 1.)*(2.44948974278318* (double) inner_prev->second.hdr.Nside + 1.) - ring*ring)*((2.44948974278318* (double) inner_prev->second.hdr.Nside - 1.)*(2.44948974278318* (double) inner_prev->second.hdr.Nside - 1.) - ring*ring)))) * inner_prev->second.hdr.distance / 3. / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside;
// (acos(1. - (ring+1) * (ring+1) / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside)) - acos(1. - (ring-1) * (ring-1) / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside))) * inner_prev->second.hdr.distance;

		//d = sin(acos(1. - ring * ring / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside))) * inner_prev->second.hdr.distance * M_PI / (2 * ring);
		d = ring * sqrt(6. - (double) (ring*ring) / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside) / 3. / (double) inner_prev->second.hdr.Nside; // sin(theta)

		result[0] = d * grad[2] + (1. - (ring*ring) / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside)) * grad[0];

		if (4*ring * v[3] > (2*j + 1) * M_PI)
		{
			k = 2 * ring * (ring-1) + (j+1) % (4*ring);

			if (inner_prev->second.pixel[k] < -1e30 || inner_next->second.pixel[k] < -1e30) return false;

			grad[1] = tan(M_PI / (4 * ring)) * (((inner_next->second.pixel[k] - inner_next->second.pixel[i]) - (inner_prev->second.pixel[k] - inner_prev->second.pixel[i])) / dtau / d / inner_prev->second.hdr.distance / (1. - cos(M_PI / (2 * ring))) + result[0]);
		}
		else
		{
			k = 2 * ring * (ring-1) + (4*ring + j-1) % (4*ring);

			if (inner_prev->second.pixel[k] < -1e30 || inner_next->second.pixel[k] < -1e30) return false;

			grad[1] = tan(M_PI / (4 * ring)) * (((inner_next->second.pixel[i] - inner_next->second.pixel[k]) - (inner_prev->second.pixel[i] - inner_prev->second.pixel[k])) / dtau / d / inner_prev->second.hdr.distance / (1. - cos(M_PI / (2 * ring))) - result[0]);
		}

		// rotate
		result[1] = sin(M_PI * (j + 0.5) / (2 * ring));
		result[2] = cos(M_PI * (j + 0.5) / (2 * ring));

		result[0] = result[0] * result[2] - result[1] * grad[1];
		result[1] = (d * grad[2] + (1. - (ring*ring) / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside)) * grad[0]) * result[1] + result[2] * grad[1];
		result[2] = (1. - (ring*ring) / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside)) * grad[2] - d * grad[0];

		//return true;
	}
	else if (i < (int64_t) inner_prev->second.hdr.Nside * ((int64_t) inner_prev->second.hdr.Nside + 1) * 2) // on northern circle
	{
#ifdef DEBUG
		if (debug) cout << " on northern circle" << endl;
#endif

		ring = (int64_t) inner_prev->second.hdr.Nside;

		if ((int64_t) inner_prev->second.hdr.Npix < 2 * ring * (ring+3)) return false;

		j = i - 2 * ring * (ring-1);

		k = i + 4 * ring;
		l = (j == ring-1) ? i+1 : k+1;

		if (inner_prev->second.pixel[k] < -1e30 || inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[k] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

		grad[0] = 0.5 * ((inner_next->second.pixel[k] + inner_next->second.pixel[l]) - (inner_prev->second.pixel[k] + inner_prev->second.pixel[l])) / dtau;

		if (j % ring == 0)
			k = ((4 + j / ring) * (ring-1) - 1) % (4 * (ring-1));
		else
			k = (j * (ring-1)) / ring;

		d = (1.5 + k) - ((double) (ring-1) * (0.5 + j) / (double) ring);

		if (j == 0) d -= 4 * (ring-1);

		l = ((k+1) % (4 * (ring-1))) + 2 * (ring-1) * (ring-2);
		k += 2 * (ring-1) * (ring-2);

		if (inner_prev->second.pixel[k] < -1e30 || inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[k] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

		grad[0] -= ((d * inner_next->second.pixel[k] + (1.-d) * inner_next->second.pixel[l]) - (d * inner_prev->second.pixel[k] + (1.-d) * inner_prev->second.pixel[l])) / dtau;

		grad[0] /= sqrt((5*ring*ring + 6 - 2./ring - (ring-1) * sqrt((ring+2)*(5*ring-2) * (double) (ring*(5*ring-2)-1)) / ring) / ring) * inner_prev->second.hdr.distance / 3. / ring;
// (acos((ring-1) / (1.5 * ring)) - acos(1. - (ring-1) * (ring-1) / (3. * ring * ring))) * inner_prev->second.hdr.distance;

		d = 0.74535599249993 * inner_prev->second.hdr.distance * (1. - cos(0.5 * M_PI / ring));
// 1.17080245517345 * inner_prev->second.hdr.distance / ring; // sin(acos(2./3.)) * M_PI / 2.;

		if (4*ring * v[3] > (2*j + 1) * M_PI)
		{
			k = 2 * ring * (ring-1) + (j+1) % (4*ring);

			if (inner_prev->second.pixel[k] < -1e30 || inner_next->second.pixel[k] < -1e30) return false;

			grad[1] = tan(M_PI / (4 * ring)) * (((inner_next->second.pixel[k] - inner_next->second.pixel[i]) - (inner_prev->second.pixel[k] - inner_prev->second.pixel[i])) / dtau / d + 0.74535599249993 * grad[2] + grad[0] / 1.5);
		}
		else
		{
			k = 2 * ring * (ring-1) + (4*ring + j-1) % (4*ring);

			if (inner_prev->second.pixel[k] < -1e30 || inner_next->second.pixel[k] < -1e30) return false;

			grad[1] = tan(M_PI / (4 * ring)) * (((inner_next->second.pixel[i] - inner_next->second.pixel[k]) - (inner_prev->second.pixel[i] - inner_prev->second.pixel[k])) / dtau / d - 0.74535599249993 * grad[2] - grad[0] / 1.5);
		}

		// rotate
		result[1] = sin(M_PI * (j + 0.5) / (2 * ring));
		result[2] = cos(M_PI * (j + 0.5) / (2 * ring));

		result[0] = (0.74535599249993 * grad[2] + grad[0] / 1.5) * result[2] - result[1] * grad[1];
		result[1] = (0.74535599249993 * grad[2] + grad[0] / 1.5) * result[1] + result[2] * grad[1];
		result[2] = grad[2] / 1.5 - 0.74535599249993 * grad[0];

		//return true;
	}
	else if (i < (int64_t) inner_prev->second.hdr.Nside * (10 * (int64_t) inner_prev->second.hdr.Nside - 2)) // in equatorial region
	{
#ifdef DEBUG
		if (debug) cout << " in equatorial region" << endl;
#endif

		j = i - 2 * (int64_t) inner_prev->second.hdr.Nside * ((int64_t) inner_prev->second.hdr.Nside + 1);
		ring = (int64_t) inner_prev->second.hdr.Nside + 1 + j / (4 * (int64_t) inner_prev->second.hdr.Nside);		

		j %= 4 * (int64_t) inner_prev->second.hdr.Nside;

		if ((int64_t) inner_prev->second.hdr.Npix < i + 4 * (int64_t) inner_prev->second.hdr.Nside + 1) return false;

		d = (2*j + 1 - ring%2) * M_PI;

		if (4*ring * v[3] > d)
		{
			d += M_PI;
			if (j+1 < 4 * (int64_t) inner_prev->second.hdr.Nside)
			{
				k = i - 4 * (int64_t) inner_prev->second.hdr.Nside + 1 - ring%2;
				j = i + 1;
			}
			else
			{
				k = i - ((ring%2) ? 4 * (int64_t) inner_prev->second.hdr.Nside : 8 * (int64_t) inner_prev->second.hdr.Nside - 1);
				j = i + 1 - 4 * (int64_t) inner_prev->second.hdr.Nside;
			}
			if (inner_prev->second.pixel[j] < -1e30 || inner_next->second.pixel[j] < -1e30) return false;
		}
		else
		{
			d -= M_PI;
			if (j > 0)
			{
				k = i - 4 * (int64_t) inner_prev->second.hdr.Nside - ring%2;
				j = i--;
			}
			else
			{
				k = i - ((ring%2) ? 1 : 4 * (int64_t) inner_prev->second.hdr.Nside);
				j = i;
				i += 4 * (int64_t) inner_prev->second.hdr.Nside - 1;
			}
			if (inner_prev->second.pixel[i] < -1e30 || inner_next->second.pixel[i] < -1e30) return false;
		}
		if (inner_prev->second.pixel[k] < -1e30 || inner_next->second.pixel[k] < -1e30) return false;

		l = k + 8 * (int64_t) inner_prev->second.hdr.Nside;

		if (inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

		grad[0] = ((inner_next->second.pixel[l] - inner_next->second.pixel[k]) - (inner_prev->second.pixel[l] - inner_prev->second.pixel[k])) / dtau / inner_prev->second.hdr.distance / (sqrt(8. + 2. * (2 * ring - inner_prev->second.hdr.Nside) * (7 * inner_prev->second.hdr.Nside - 2 * ring) - 2. * sqrt((double) (2 * ring - 2 - inner_prev->second.hdr.Nside) * (double) (2 * ring + 2 - inner_prev->second.hdr.Nside) * (double) (7 * inner_prev->second.hdr.Nside + 2 - 2 * ring) * (double) (7 * inner_prev->second.hdr.Nside - 2 - 2 * ring))) / 3. / (double) inner_prev->second.hdr.Nside);
// (acos((2 * (int64_t) inner_prev->second.hdr.Nside - ring - 1) / (1.5 * inner_prev->second.hdr.Nside)) - acos((2 * (int64_t) inner_prev->second.hdr.Nside - ring + 1) / (1.5 * inner_prev->second.hdr.Nside)));

		result[1] = sin(d / (4 * ring));
		result[2] = cos(d / (4 * ring));

		d = sqrt(1. - (double) (2 * inner_prev->second.hdr.Nside - ring) * (double) (2 * inner_prev->second.hdr.Nside - ring) / (2.25 * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside)); // sin(theta)

		grad[1] = ((inner_next->second.pixel[j] - inner_next->second.pixel[i]) - (inner_prev->second.pixel[j] - inner_prev->second.pixel[i])) / dtau / inner_prev->second.hdr.distance / (2. * d * sin(0.25 * M_PI / (double) inner_prev->second.hdr.Nside));
// (sin(acos((2 * (int64_t) inner_prev->second.hdr.Nside - ring) / (1.5 * inner_prev->second.hdr.Nside))) * M_PI / (2 * (int64_t) inner_prev->second.hdr.Nside));

		// rotate

		result[0] = (d * grad[2] + ((2 * inner_prev->second.hdr.Nside - ring) / 1.5 / (double) inner_prev->second.hdr.Nside) * grad[0]) * result[2] - result[1] * grad[1];
		result[1] = (d * grad[2] + ((2 * inner_prev->second.hdr.Nside - ring) / 1.5 / (double) inner_prev->second.hdr.Nside) * grad[0]) * result[1] + result[2] * grad[1];
		result[2] = ((2 * inner_prev->second.hdr.Nside - ring) / 1.5 / (double) inner_prev->second.hdr.Nside) * grad[2] - d * grad[0];

		//return true;
	}
	else if (i < (int64_t) inner_prev->second.hdr.Nside * (10 * (int64_t) inner_prev->second.hdr.Nside + 2)) // on southern circle
	{
#ifdef DEBUG
		if (debug) cout << " on southern circle" << endl;
#endif

		ring = (int64_t) inner_prev->second.hdr.Nside; // counted from south pole

		if ((int64_t) inner_prev->second.hdr.Npix < ring * (10 * ring + 6) - 4) return false;

		j = i - ring * (10 * ring - 2);

		if (j % ring == 0)
			k = ((4 + j / ring) * (ring-1) - 1) % (4 * (ring-1));
		else
			k = (j * (ring-1)) / ring;

		d = (1.5 + k) - ((double) (ring-1) * (0.5 + j) / (double) ring);

		if (j == 0) d -= 4 * (ring-1);

		l = ((k+1) % (4 * (ring-1))) + ring * (10 * ring + 2);
		k += ring * (10 * ring + 2);

		if (inner_prev->second.pixel[k] < -1e30 || inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[k] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

		grad[0] = ((d * inner_next->second.pixel[k] + (1.-d) * inner_next->second.pixel[l]) - (d * inner_prev->second.pixel[k] + (1.-d) * inner_prev->second.pixel[l])) / dtau;

		k = i - 4 * ring;
		l = (j == ring-1) ? k+1 - 4 * ring : k+1;

		if (inner_prev->second.pixel[k] < -1e30 || inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[k] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

		grad[0] -= 0.5 * ((inner_next->second.pixel[k] + inner_next->second.pixel[l]) - (inner_prev->second.pixel[k] + inner_prev->second.pixel[l])) / dtau;

		grad[0] /= sqrt((5*ring*ring + 6 - 2./ring - (ring-1) * sqrt((ring+2)*(5*ring-2) * (double) (ring*(5*ring-2)-1)) / ring) / ring) * inner_prev->second.hdr.distance / 3. / ring;
// (acos(-1. + (ring-1) * (ring-1) / (3. * ring * ring)) - acos((1-ring) / (1.5 * ring))) * inner_prev->second.hdr.distance;

		d = 0.74535599249993 * inner_prev->second.hdr.distance * (1. - cos(0.5 * M_PI / ring));
// 1.17080245517345 * inner_prev->second.hdr.distance / ring; // sin(acos(-2./3.)) * M_PI / 2.;

		if (4*ring * v[3] > (2*j + 1) * M_PI)
		{
			k = ring * (10 * ring - 2) + (j+1) % (4*ring);

			if (inner_prev->second.pixel[k] < -1e30 || inner_next->second.pixel[k] < -1e30) return false;

			grad[1] = tan(M_PI / (4 * ring)) * (((inner_next->second.pixel[k] - inner_next->second.pixel[i]) -  (inner_prev->second.pixel[k] - inner_prev->second.pixel[i])) / dtau / d + 0.74535599249993 * grad[2] - grad[0] / 1.5);
		}
		else
		{
			k = ring * (10 * ring - 2) + (4*ring + j-1) % (4*ring);

			if (inner_prev->second.pixel[k] < -1e30 || inner_next->second.pixel[k] < -1e30) return false;

			grad[1] = tan(M_PI / (4 * ring)) * (((inner_next->second.pixel[i] - inner_next->second.pixel[k]) - (inner_prev->second.pixel[i] - inner_prev->second.pixel[k])) / dtau / d - 0.74535599249993 * grad[2] + grad[0] / 1.5);
		}

		// rotate
		result[1] = sin(M_PI * (j + 0.5) / (2 * ring));
		result[2] = cos(M_PI * (j + 0.5) / (2 * ring));

		result[0] = (0.74535599249993 * grad[2] - grad[0] / 1.5) * result[2] - result[1] * grad[1];
		result[1] = (0.74535599249993 * grad[2] - grad[0] / 1.5) * result[1] + result[2] * grad[1];
		result[2] = -grad[2] / 1.5 - 0.74535599249993 * grad[0];

		//return true;
	}
	else if (i < 12 * (int64_t) inner_prev->second.hdr.Nside * (int64_t) inner_prev->second.hdr.Nside - 4) // in south polar cap
	{
#ifdef DEBUG
		if (debug) cout << " in south polar cap" << endl;
#endif

		Npix = nside2npix64((int64_t) inner_prev->second.hdr.Nside);
		ring = (1 + (int64_t) sqrt(2 * (Npix-i) - 0.5)) / 2; // counted from south pole

		if ((int64_t) inner_prev->second.hdr.Npix < Npix - 2 * ring * (ring-1)) return false;

		j = i + 2 * ring * (ring+1) - Npix;

		if (j % ring == 0)
			k = ((4 + j / ring) * (ring-1) - 1) % (4 * (ring-1));
		else
			k = (j * (ring-1)) / ring;

		d = (1.5 + k) - ((double) (ring-1) * (0.5 + j) / (double) ring);

		if (j == 0) d -= 4 * (ring-1);

		l = ((k+1) % (4 * (ring-1))) + Npix - 2 * ring * (ring-1);
		k += Npix - 2 * ring * (ring-1);

		if (inner_prev->second.pixel[k] < -1e30 || inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[k] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

		grad[0] = ((d * inner_next->second.pixel[k] + (1.-d) * inner_next->second.pixel[l]) - (d * inner_prev->second.pixel[k] + (1.-d) * inner_prev->second.pixel[l])) / dtau;

		k = (j * (ring+1)) / ring;

		d = (1.5 + k) - ((double) (ring+1) * (0.5 + j) / (double) ring);

		k += Npix - 2 * (ring+1) * (ring+2);
		l = k + 1;

		if (inner_prev->second.pixel[k] < -1e30 || inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[k] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

		grad[0] -= ((d * inner_next->second.pixel[k] + (1.-d) * inner_next->second.pixel[l]) - (d * inner_prev->second.pixel[k] + (1.-d) * inner_prev->second.pixel[l])) / dtau;

		grad[0] /= sqrt(12. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside * (ring*ring+1) - 2. * (ring*ring-1) * (ring*ring - 1 + sqrt(((2.44948974278318* (double) inner_prev->second.hdr.Nside + 1.)*(2.44948974278318* (double) inner_prev->second.hdr.Nside + 1.) - ring*ring)*((2.44948974278318* (double) inner_prev->second.hdr.Nside - 1.)*(2.44948974278318* (double) inner_prev->second.hdr.Nside - 1.) - ring*ring)))) * inner_prev->second.hdr.distance / 3. / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside;
// (acos(1. - (ring+1) * (ring+1) / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside)) - acos(1. - (ring-1) * (ring-1) / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside))) * inner_prev->second.hdr.distance;

		//d = sin(acos(1. - ring * ring / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside))) * inner_prev->second.hdr.distance * M_PI / (2 * ring);
		//d = ring * sqrt(6. - (double) (ring*ring) / (double) inner_prev->second.hdr.Nside / inner_prev->second.hdr.Nside) * inner_prev->second.hdr.distance * sin(M_PI / (4 * ring)) / 1.5 / inner_prev->second.hdr.Nside;
		d = ring * sqrt(6. - (double) (ring*ring) / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside) / 3. / (double) inner_prev->second.hdr.Nside; // sin(theta)

		result[0] = d * grad[2] - (1. - (ring*ring) / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside)) * grad[0];

		if (4*ring * v[3] > (2*j + 1) * M_PI)
		{
			k = Npix - 2 * ring * (ring+1) + (j+1) % (4*ring);

			if (inner_prev->second.pixel[k] < -1e30 || inner_next->second.pixel[k] < -1e30) return false;

			grad[1] = tan(M_PI / (4 * ring)) * (((inner_next->second.pixel[k] - inner_next->second.pixel[i]) - (inner_prev->second.pixel[k] - inner_prev->second.pixel[i])) / dtau / d / inner_prev->second.hdr.distance / (1. - cos(M_PI / (2 * ring))) + result[0]);
		}
		else
		{
			k = Npix - 2 * ring * (ring+1) + (4*ring + j-1) % (4*ring);

			if (inner_prev->second.pixel[k] < -1e30 || inner_next->second.pixel[k] < -1e30) return false;

			grad[1] = tan(M_PI / (4 * ring)) * (((inner_next->second.pixel[i] - inner_next->second.pixel[k]) - (inner_prev->second.pixel[i] - inner_prev->second.pixel[k])) / dtau / d / inner_prev->second.hdr.distance / (1. - cos(M_PI / (2 * ring))) - result[0]);
		}

		// rotate
		result[1] = sin(M_PI * (j + 0.5) / (2 * ring));
		result[2] = cos(M_PI * (j + 0.5) / (2 * ring));

		result[0] = result[0] * result[2] - result[1] * grad[1];
		result[1] = (d * grad[2] - (1. - (ring*ring) / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside)) * grad[0]) * result[1] + result[2] * grad[1];
		result[2] = (-1. + (ring*ring) / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside)) * grad[2] - d * grad[0];

		//return true;
	}
	else // at south pole
	{
#ifdef DEBUG
		if (debug) cout << " at south pole" << endl;
#endif

		Npix = nside2npix64((int64_t) inner_prev->second.hdr.Nside);
		if ((int64_t) inner_prev->second.hdr.Npix < Npix || (int64_t) inner_next->second.hdr.Npix < Npix) return false;
		d = -1. + 1. / (3. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside);
		if (v[2] < d)
		{
			if (inner_prev->second.pixel[Npix-1] < -1e30 || inner_prev->second.pixel[Npix-2] < -1e30 || inner_prev->second.pixel[Npix-3] < -1e30 || inner_prev->second.pixel[Npix-4] < -1e30) return false;
			if (inner_next->second.pixel[Npix-1] < -1e30 || inner_next->second.pixel[Npix-2] < -1e30 || inner_next->second.pixel[Npix-3] < -1e30 || inner_next->second.pixel[Npix-4] < -1e30) return false;

			//d = sqrt(8.) * acos(d) * inner_prev->second.hdr.distance;
			d = sqrt(8. - 8. * d * d) * inner_prev->second.hdr.distance;
			result[0] = ((inner_next->second.pixel[Npix-4] - inner_next->second.pixel[Npix-2]) - (inner_prev->second.pixel[Npix-4] - inner_prev->second.pixel[Npix-2])) / dtau / d;
			result[1] = ((inner_next->second.pixel[Npix-3] - inner_next->second.pixel[Npix-1]) - (inner_prev->second.pixel[Npix-3] - inner_prev->second.pixel[Npix-1])) / dtau / d;

			d = result[0] - result[1]; // rotate
			result[1] += result[0];
			result[0] = d;
			result[2] = -grad[2];

			//return true;
		}
		else if (inner_prev->second.hdr.Nside == 2)
		{
			j = i-44;

			if (inner_prev->second.pixel[2*j+36] < -1e30 || inner_prev->second.pixel[2*j+37] < -1e30 || inner_prev->second.pixel[2*j+29] < -1e30) return false;
			if (inner_next->second.pixel[2*j+36] < -1e30 || inner_next->second.pixel[2*j+37] < -1e30 || inner_next->second.pixel[2*j+29] < -1e30) return false;

			grad[0] = ((inner_next->second.pixel[i] - inner_next->second.pixel[2*j+29]) - (inner_prev->second.pixel[i] - inner_prev->second.pixel[2*j+29])) / dtau /
 inner_prev->second.hdr.distance / 0.797054997187545; // 2. * sin(0.5 * (acos(1./3.) - acos(11./12.))) / 0.819821555018427; // (acos(1./3.) - acos(11./12.))

			grad[1] = ((inner_next->second.pixel[2*j+37] - inner_next->second.pixel[2*j+36]) - (inner_prev->second.pixel[2*j+37] - inner_prev->second.pixel[2*j+36])) / dtau / inner_prev->second.hdr.distance / 0.570470779087523; // (sin(acos(2./3.)) * 2. * sin(M_PI / 8.) / 0.585401227586727; // (sin(acos(2./3.)) * MP_I / 4.)

			// rotate
			result[0] = ((i % 2) ? (-0.52704627669473 * grad[2] + 0.471404520791031 * grad[0]) : (0.52704627669473 * grad[2] - 0.471404520791031 * grad[0])) - 0.707106781186547 * grad[1];
			result[1] = 0.52704627669473 * grad[2] - 0.471404520791031 * grad[0] + ((i % 2) ? -0.707106781186547 : 0.707106781186547) * grad[1];
		
			if (i > 1)
			{
				result[0] *= -1;
				result[1] *= -1;
			}

			result[2] = -(grad[2] / 1.5) - 0.74535599249993 * grad[0];

			//return true;
		}
		else // Nside > 2
		{
			j = i + 4 - Npix;

			if (inner_prev->second.pixel[2*j+Npix-12] < -1e30 || inner_prev->second.pixel[2*j+Npix-11] < -1e30 || inner_prev->second.pixel[3*j+Npix-23] < -1e30) return false;
			if (inner_next->second.pixel[2*j+Npix-12] < -1e30 || inner_next->second.pixel[2*j+Npix-11] < -1e30 || inner_next->second.pixel[3*j+Npix-23] < -1e30) return false;

			grad[0] = ((inner_next->second.pixel[i] - inner_next->second.pixel[3*j+Npix-23]) - (inner_prev->second.pixel[i] - inner_prev->second.pixel[3*j+Npix-23])) / dtau / inner_prev->second.hdr.distance / (sqrt(6.66666666666666 - sqrt(16. + (4. / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside - 26.6666666666666) / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside) - 2. / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside) / (double) inner_prev->second.hdr.Nside);
// (acos(1. - 3./(inner_prev->second.hdr.Nside * inner_prev->second.hdr.Nside)) - acos(1. - 1./(3. * inner_prev->second.hdr.Nside * inner_prev->second.hdr.Nside)));

			d = sqrt(6. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside - 4.); // sin(theta)

			grad[1] = ((inner_next->second.pixel[2*j+Npix-11] - inner_next->second.pixel[2*j+Npix-12]) - (inner_prev->second.pixel[2*j+Npix-11] - inner_prev->second.pixel[2*j+Npix-12])) / dtau / inner_prev->second.hdr.distance / (d * 0.5102445764867864 / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside);
// (sin(acos(1. - 4. / (3. * inner_prev->second.hdr.Nside * inner_prev->second.hdr.Nside))) * M_PI / 4.);

			//rotate

			result[0] = ((i % 2) ? -0.471404520791032 : 0.471404520791032) * (d * grad[2] - (1.5 * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside - 2.) * grad[0]) / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside - 0.707106781186547 * grad[1];
			result[1] = 0.471404520791032 * (d * grad[2] - (1.5 * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside - 2.) * grad[0]) / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside + ((i % 2) ? -0.707106781186547 : 0.707106781186547) * grad[1];

			if (i > 1)
			{
				result[0] *= -1;
				result[1] *= -1;
			}

			result[2] = ((2. - 1.5 * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside) * grad[2] - d * grad[0]) / 1.5 / (double) inner_prev->second.hdr.Nside / (double) inner_prev->second.hdr.Nside;

			//return true;
		}
	}

/*	d = sqrt(1. - v[2] * v[2]); // sin(theta)

	result[0] = v[0] * grad[2] + (v[0] * v[2] * grad[0] - v[1] * grad[1]) / d;
	result[1] = v[1] * grad[2] + (v[1] * v[2] * grad[0] + v[0] * grad[1]) / d;
	result[2] = v[2] * grad[2] - d * grad[0];*/

	return true;
}


bool tidal_matrix(metric_container * field_prev, metric_container * field_next, const double weight, const double dist, const double dr, double * v, double * result, int debug)
{
	std::map<int,metric_data>::iterator inner_prev = get_iterator(field_prev, dist, -1);
	std::map<int,metric_data>::iterator inner_next = get_iterator(field_next, dist, -1);
	std::map<int,metric_data>::iterator outer_prev = get_iterator(field_prev, dist);
	std::map<int,metric_data>::iterator outer_next = get_iterator(field_next, dist);
	std::map<int,metric_data>::iterator center_prev;
	std::map<int,metric_data>::iterator center_next;
	double d, cth, sth, cph, sph, tmp, ppp, mpp, mmp; //, wr, tmp2, ppp2, mpp2, mmp2, sph2, sth2, cph2;
	int64_t i, j, k, l, m, n, ring, Npix;
	double tidal[6];

	if (inner_prev == field_prev->healpix_data.end() || inner_next == field_next->healpix_data.end() || outer_prev == field_prev->healpix_data.end() || outer_next == field_next->healpix_data.end()) return false;

	if (inner_prev->second.hdr.distance != inner_next->second.hdr.distance || inner_prev->second.hdr.Nside != inner_next->second.hdr.Nside) return false;

	if (outer_prev->second.hdr.distance != outer_next->second.hdr.distance || outer_prev->second.hdr.Nside != outer_next->second.hdr.Nside) return false;

	//wr = modf(dist / dr, &d);

	if(modf(dist / dr, &d) >= 0.5)
	{
		center_prev = outer_prev;
		center_next = outer_next;
		
		outer_prev++; outer_next++;

		if (outer_prev == field_prev->healpix_data.end() || outer_next == field_next->healpix_data.end()) return false;

		if (outer_prev->second.hdr.distance != outer_next->second.hdr.distance || outer_prev->second.hdr.Nside != outer_next->second.hdr.Nside) return false;

		d = -1.;
	}
	else if (d > 0)
	{
		if (inner_prev->first == 0 || inner_next->first == 0) return false;

		center_prev = inner_prev;
		center_next = inner_next;

		if (inner_prev == field_prev->healpix_data.begin() || inner_next == field_next->healpix_data.begin()) return false;

		inner_prev--; inner_next--;

		if (inner_prev->second.hdr.distance != inner_next->second.hdr.distance || inner_prev->second.hdr.Nside != inner_next->second.hdr.Nside) return false;
	}
	else // at vertex
	{
		//if (debug) cout << " at vertex" << endl;
		if (inner_prev->second.hdr.distance != 0) return false;

		Npix = nside2npix64((int64_t) outer_prev->second.hdr.Nside);

		tidal[1] = 2. * linear_interpolation(inner_prev->second.pixel[0], inner_next->second.pixel[0], weight);

		i = Npix/2;
		j = i-1;
		k = i + 4*outer_prev->second.hdr.Nside;
		l = i - 4*outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[l] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[i] < -1e30 || outer_prev->second.pixel[k] < -1e30) return false;
		if (outer_next->second.pixel[l] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[i] < -1e30 || outer_next->second.pixel[k] < -1e30) return false;

		tidal[0] = linear_interpolation(outer_prev->second.pixel[l] + outer_prev->second.pixel[j] + outer_prev->second.pixel[i] + outer_prev->second.pixel[k], outer_next->second.pixel[l] + outer_next->second.pixel[j] + outer_next->second.pixel[i] + outer_next->second.pixel[k], weight);

		i -= 2*outer_prev->second.hdr.Nside;
		j += 2*outer_prev->second.hdr.Nside;
		k -= 2*outer_prev->second.hdr.Nside;
		l -= 2*outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[l] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[i] < -1e30 || outer_prev->second.pixel[k] < -1e30) return false;
		if (outer_next->second.pixel[l] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[i] < -1e30 || outer_next->second.pixel[k] < -1e30) return false;

		tidal[0] += linear_interpolation(outer_prev->second.pixel[l] + outer_prev->second.pixel[j] + outer_prev->second.pixel[i] + outer_prev->second.pixel[k], outer_next->second.pixel[l] + outer_next->second.pixel[j] + outer_next->second.pixel[i] + outer_next->second.pixel[k], weight);

		tidal[0] = (0.25 * tidal[0] - tidal[1]) / dr / dr;

		i += outer_prev->second.hdr.Nside;
		j -= outer_prev->second.hdr.Nside;
		k += outer_prev->second.hdr.Nside;
		l += outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[l] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[i] < -1e30 || outer_prev->second.pixel[k] < -1e30) return false;
		if (outer_next->second.pixel[l] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[i] < -1e30 || outer_next->second.pixel[k] < -1e30) return false;

		tidal[3] = linear_interpolation(outer_prev->second.pixel[l] + outer_prev->second.pixel[j] + outer_prev->second.pixel[i] + outer_prev->second.pixel[k], outer_next->second.pixel[l] + outer_next->second.pixel[j] + outer_next->second.pixel[i] + outer_next->second.pixel[k], weight);

		i += 2*outer_prev->second.hdr.Nside;
		j -= 2*outer_prev->second.hdr.Nside;
		k += 2*outer_prev->second.hdr.Nside;
		l += 2*outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[l] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[i] < -1e30 || outer_prev->second.pixel[k] < -1e30) return false;
		if (outer_next->second.pixel[l] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[i] < -1e30 || outer_next->second.pixel[k] < -1e30) return false;

		tidal[3] += linear_interpolation(outer_prev->second.pixel[l] + outer_prev->second.pixel[j] + outer_prev->second.pixel[i] + outer_prev->second.pixel[k], outer_next->second.pixel[l] + outer_next->second.pixel[j] + outer_next->second.pixel[i] + outer_next->second.pixel[k], weight);

		tidal[3] = (0.25 * tidal[3] - tidal[1]) / dr / dr;

		if (outer_prev->second.pixel[0] < -1e30 || outer_prev->second.pixel[1] < -1e30 || outer_prev->second.pixel[2] < -1e30 || outer_prev->second.pixel[3] < -1e30) return false;
		if (outer_next->second.pixel[0] < -1e30 || outer_next->second.pixel[1] < -1e30 || outer_next->second.pixel[2] < -1e30 || outer_next->second.pixel[3] < -1e30) return false;

		tidal[5] = linear_interpolation(outer_prev->second.pixel[0] + outer_prev->second.pixel[1] + outer_prev->second.pixel[2] + outer_prev->second.pixel[3], outer_next->second.pixel[0] + outer_next->second.pixel[1] + outer_next->second.pixel[2] + outer_next->second.pixel[3], weight);

		if (outer_prev->second.pixel[Npix-4] < -1e30 || outer_prev->second.pixel[Npix-3] < -1e30 || outer_prev->second.pixel[Npix-2] < -1e30 || outer_prev->second.pixel[Npix-1] < -1e30) return false;
		if (outer_next->second.pixel[Npix-4] < -1e30 || outer_next->second.pixel[Npix-3] < -1e30 || outer_next->second.pixel[Npix-2] < -1e30 || outer_next->second.pixel[Npix-1] < -1e30) return false;

		tidal[5] += linear_interpolation(outer_prev->second.pixel[Npix-4] + outer_prev->second.pixel[Npix-3] + outer_prev->second.pixel[Npix-2] + outer_prev->second.pixel[Npix-1], outer_next->second.pixel[Npix-4] + outer_next->second.pixel[Npix-3] + outer_next->second.pixel[Npix-2] + outer_next->second.pixel[Npix-1], weight);

		tidal[5] = (0.25 * tidal[5] - tidal[1]) / dr / dr;

		i += outer_prev->second.hdr.Nside/2;
		j = i - 1;
		k = i + 4*outer_prev->second.hdr.Nside;
		l = i - 4*outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[l] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[i] < -1e30 || outer_prev->second.pixel[k] < -1e30) return false;
		if (outer_next->second.pixel[l] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[i] < -1e30 || outer_next->second.pixel[k] < -1e30) return false;

		tidal[1] = linear_interpolation(outer_prev->second.pixel[l] + outer_prev->second.pixel[j] + outer_prev->second.pixel[i] + outer_prev->second.pixel[k], outer_next->second.pixel[l] + outer_next->second.pixel[j] + outer_next->second.pixel[i] + outer_next->second.pixel[k], weight);

		i -= outer_prev->second.hdr.Nside;
		j -= outer_prev->second.hdr.Nside;
		k -= outer_prev->second.hdr.Nside;
		l -= outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[l] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[i] < -1e30 || outer_prev->second.pixel[k] < -1e30) return false;
		if (outer_next->second.pixel[l] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[i] < -1e30 || outer_next->second.pixel[k] < -1e30) return false;

		tidal[1] -= linear_interpolation(outer_prev->second.pixel[l] + outer_prev->second.pixel[j] + outer_prev->second.pixel[i] + outer_prev->second.pixel[k], outer_next->second.pixel[l] + outer_next->second.pixel[j] + outer_next->second.pixel[i] + outer_next->second.pixel[k], weight);

		i -= outer_prev->second.hdr.Nside;
		j -= outer_prev->second.hdr.Nside;
		k -= outer_prev->second.hdr.Nside;
		l -= outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[l] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[i] < -1e30 || outer_prev->second.pixel[k] < -1e30) return false;
		if (outer_next->second.pixel[l] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[i] < -1e30 || outer_next->second.pixel[k] < -1e30) return false;

		tidal[1] += linear_interpolation(outer_prev->second.pixel[l] + outer_prev->second.pixel[j] + outer_prev->second.pixel[i] + outer_prev->second.pixel[k], outer_next->second.pixel[l] + outer_next->second.pixel[j] + outer_next->second.pixel[i] + outer_next->second.pixel[k], weight);

		i -= outer_prev->second.hdr.Nside;
		j -= outer_prev->second.hdr.Nside;
		k -= outer_prev->second.hdr.Nside;
		l -= outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[l] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[i] < -1e30 || outer_prev->second.pixel[k] < -1e30) return false;
		if (outer_next->second.pixel[l] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[i] < -1e30 || outer_next->second.pixel[k] < -1e30) return false;

		tidal[1] -= linear_interpolation(outer_prev->second.pixel[l] + outer_prev->second.pixel[j] + outer_prev->second.pixel[i] + outer_prev->second.pixel[k], outer_next->second.pixel[l] + outer_next->second.pixel[j] + outer_next->second.pixel[i] + outer_next->second.pixel[k], weight);

		tidal[1] /= -8. * dr * dr;

		i = 2 * (int64_t) outer_prev->second.hdr.Nside * (outer_prev->second.hdr.Nside + 1);
		j = i - 1;
		k = i - 4*outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[i] < -1e30) return false;
		if (outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[i] < -1e30) return false;

		tidal[2] = linear_interpolation(outer_prev->second.pixel[k] + outer_prev->second.pixel[j] + outer_prev->second.pixel[i], outer_next->second.pixel[k] + outer_next->second.pixel[j] + outer_next->second.pixel[i], weight);

		i += 2*outer_prev->second.hdr.Nside;
		j -= 2*outer_prev->second.hdr.Nside;
		k += 2*outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[i] < -1e30) return false;
		if (outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[i] < -1e30) return false;

		tidal[2] -= linear_interpolation(outer_prev->second.pixel[j] + outer_prev->second.pixel[k] + outer_prev->second.pixel[i], outer_next->second.pixel[j] + outer_next->second.pixel[k] + outer_next->second.pixel[i], weight);

		i += 8l * (int64_t) outer_prev->second.hdr.Nside * (outer_prev->second.hdr.Nside - 1);
		k = i + 4*outer_prev->second.hdr.Nside;
		j = k - 1;

		if (outer_prev->second.pixel[i] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[k] < -1e30) return false;
		if (outer_next->second.pixel[i] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[k] < -1e30) return false;

		tidal[2] += linear_interpolation(outer_prev->second.pixel[i] + outer_prev->second.pixel[j] + outer_prev->second.pixel[k], outer_next->second.pixel[i] + outer_next->second.pixel[j] + outer_next->second.pixel[k], weight);

		i -= 2*outer_prev->second.hdr.Nside;
		k += 2*outer_prev->second.hdr.Nside;
		j -= 2*outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[i] < -1e30 || outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[j] < -1e30) return false;
		if (outer_next->second.pixel[i] < -1e30 || outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[j] < -1e30) return false;

		tidal[2] -= linear_interpolation(outer_prev->second.pixel[i] + outer_prev->second.pixel[k] + outer_prev->second.pixel[j], outer_next->second.pixel[i] + outer_next->second.pixel[k] + outer_next->second.pixel[j], weight);

		i = (int64_t) outer_prev->second.hdr.Nside * (2 * outer_prev->second.hdr.Nside + 3);
		j = i - 4*outer_prev->second.hdr.Nside;
		k = j - 1;

		if (outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[i] < -1e30) return false;
		if (outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[i] < -1e30) return false;

		tidal[4] = linear_interpolation(outer_prev->second.pixel[k] + outer_prev->second.pixel[j] + outer_prev->second.pixel[i], outer_next->second.pixel[k] + outer_next->second.pixel[j] + outer_next->second.pixel[i], weight);

		i += 2*outer_prev->second.hdr.Nside;
		j += 2*outer_prev->second.hdr.Nside;
		k += 2*outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[i] < -1e30) return false;
		if (outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[i] < -1e30) return false;

		tidal[4] -= linear_interpolation(outer_prev->second.pixel[k] + outer_prev->second.pixel[j] + outer_prev->second.pixel[i], outer_next->second.pixel[k] + outer_next->second.pixel[j] + outer_next->second.pixel[i], weight);

		i += 8l * (int64_t) outer_prev->second.hdr.Nside * ((int64_t) outer_prev->second.hdr.Nside - 1);
		k = i + 4*outer_prev->second.hdr.Nside;
		j = k - 1;

		if (outer_prev->second.pixel[i] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[k] < -1e30) return false;
		if (outer_next->second.pixel[i] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[k] < -1e30) return false;

		tidal[4] += linear_interpolation(outer_prev->second.pixel[i] + outer_prev->second.pixel[j] + outer_prev->second.pixel[k], outer_next->second.pixel[i] + outer_next->second.pixel[j] + outer_next->second.pixel[k], weight);

		i -= 2*outer_prev->second.hdr.Nside;
		j -= 2*outer_prev->second.hdr.Nside;
		k -= 2*outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[i] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[k] < -1e30) return false;
		if (outer_next->second.pixel[i] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[k] < -1e30) return false;

		tidal[4] -= linear_interpolation(outer_prev->second.pixel[i] + outer_prev->second.pixel[j] + outer_prev->second.pixel[k], outer_next->second.pixel[i] + outer_next->second.pixel[j] + outer_next->second.pixel[k], weight);

		d = (6. * (double) outer_prev->second.hdr.Nside - 2.) / (9. * (double) outer_prev->second.hdr.Nside);
		d *= 12. * sqrt(1. - d*d) * dr * dr;

		tidal[2] /= d;
		tidal[4] /= d;

		for (i = 0; i < 6; i++)
			result[i] = tidal[i];

		return true;
	}

	vec2pix_ring64(inner_prev->second.hdr.Nside, v, &i);
	vec2pix_ring64(center_prev->second.hdr.Nside, v, &j);
	vec2pix_ring64(outer_prev->second.hdr.Nside, v, &k);

	if (inner_prev->second.hdr.Npix <= i || center_prev->second.hdr.Npix <= j || outer_prev->second.hdr.Npix <= k) return false;
	if (inner_next->second.hdr.Npix <= i || center_next->second.hdr.Npix <= j || outer_next->second.hdr.Npix <= k) return false;
	if (inner_prev->second.pixel[i] < -1e30 || center_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[k] < -1e30) return false;
	if (inner_next->second.pixel[i] < -1e30 || center_next->second.pixel[j] < -1e30 || outer_next->second.pixel[k] < -1e30) return false;

	tidal[5] = linear_interpolation(inner_prev->second.pixel[i], inner_next->second.pixel[i], weight);
	tidal[4] = linear_interpolation(center_prev->second.pixel[j], center_next->second.pixel[j], weight);
	tidal[0] = linear_interpolation(outer_prev->second.pixel[k], outer_next->second.pixel[k], weight) - 2. * tidal[4] + tidal[5];
	tidal[0] /= dr * dr;

	tidal[3] = tidal[4] - tidal[5];
	tidal[3] /= dr * center_prev->second.hdr.distance;

	tidal[5] = tidal[3];

	/*if (d > 0)
		result[4] = linear_interpolation(outer_prev->second.pixel[k], outer_next->second.pixel[k], weight);
	else
		result[4] = linear_interpolation(inner_prev->second.pixel[i], inner_next->second.pixel[i], weight);

	result[5] = tidal[5];
	result[3] = tidal[3];*/

	if (j < 4) // at north pole
	{
		//if (debug) cout << " at north pole" << endl;
		if (center_prev->second.hdr.Npix < 24 || center_next->second.hdr.Npix < 24) return false;
		if (center_prev->second.pixel[(j+1)%4] < -1e30 || center_prev->second.pixel[(j+3)%4] < -1e30 || center_next->second.pixel[(j+1)%4] < -1e30 || center_next->second.pixel[(j+3)%4] < -1e30) return false;
		if (center_prev->second.pixel[2*j+4] < -1e30 || center_prev->second.pixel[2*j+5] < -1e30 || center_next->second.pixel[2*j+4] < -1e30 || center_next->second.pixel[2*j+5] < -1e30) return false;
		if (center_prev->second.pixel[3*j+13] < -1e30 || center_next->second.pixel[3*j+13] < -1e30) return false;

		cth = 1. - 1. / (3. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside);
		sth = sqrt(1. - cth*cth);

		sph = (j > 1) ? -0.7071067811865475 : 0.7071067811865475;
		cph = (j % 2) ? -sph : sph;

		mmp = 1. - 3. / (double) center_prev->second.hdr.Nside / (double) center_prev->second.hdr.Nside;  // mu3
		ppp = 1. - 4. / 3. / (double) center_prev->second.hdr.Nside / (double) center_prev->second.hdr.Nside;
		ppp = sqrt(1. - ppp*ppp);
		mpp = (sqrt(1. - mmp*mmp) + sth) * sth; // d1(d1+d3)

		tidal[1] = linear_interpolation(center_prev->second.pixel[2*j+4] + center_prev->second.pixel[2*j+5], center_next->second.pixel[2*j+4] + center_next->second.pixel[2*j+5], weight);
		tidal[2] = linear_interpolation(center_prev->second.pixel[(j+1)%4] + center_prev->second.pixel[(j+3)%4], center_next->second.pixel[(j+1)%4] + center_next->second.pixel[(j+3)%4], weight);

		/*tidal[3] /= 3. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside - 0.5;
		tidal[3] += (0.25 * (tidal[1] + tidal[2]) - tidal[4]) / center_prev->second.hdr.distance / center_prev->second.hdr.distance;
		tidal[3] *= 3. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside;

		tidal[5] += (3. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside - 1.) * (tidal[4] - 0.5 * tidal[2]) / center_prev->second.hdr.distance / center_prev->second.hdr.distance;
		tidal[5] *= sqrt(36. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside - 6);
		tidal[5] += 9. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside * (tidal[2] - 2. * tidal[4]) / center_prev->second.hdr.distance / center_prev->second.hdr.distance;
		tidal[5] *= (double) center_prev->second.hdr.Nside / (6. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside - 1.);*/

		tidal[3] += (((4.*cth*cth - 2.) * (linear_interpolation(center_prev->second.pixel[3*j+13], center_next->second.pixel[3*j+13], weight) - tidal[4]) - (mpp + 2.*cth*cth - 2.) * (linear_interpolation(center_prev->second.pixel[(j+2)%4], center_next->second.pixel[(j+2)%4], weight) - tidal[4])) / (cth*cth - mmp*mmp) + (0.5 + cth*cth/mpp) * (linear_interpolation(center_prev->second.pixel[(j+2)%4], center_next->second.pixel[(j+2)%4], weight) - tidal[4])) / center_prev->second.hdr.distance / center_prev->second.hdr.distance;

		tidal[5] += (((2.*cth*cth * (linear_interpolation(center_prev->second.pixel[3*j+13], center_next->second.pixel[3*j+13], weight) - tidal[4]) - (mpp + 2.*cth*cth - 2.) * (linear_interpolation(center_prev->second.pixel[(j+2)%4], center_next->second.pixel[(j+2)%4], weight) - tidal[4])) / (cth*cth - mmp*mmp) + (0.5 + 1./mpp) * (linear_interpolation(center_prev->second.pixel[(j+2)%4], center_next->second.pixel[(j+2)%4], weight) - tidal[4])) - (tidal[4] - linear_interpolation(center_prev->second.pixel[(j+1)%4] - center_prev->second.pixel[(j+2)%4] + center_prev->second.pixel[(j+3)%4], center_next->second.pixel[(j+1)%4] - center_next->second.pixel[(j+2)%4] + center_next->second.pixel[(j+3)%4], weight)) / (1. - cth*cth)) / center_prev->second.hdr.distance / center_prev->second.hdr.distance;

		tidal[1] = 0.306186217847897 * (tidal[2] - tidal[1]) * (double) center_prev->second.hdr.Nside / center_prev->second.hdr.distance;

		tidal[2] = 1.95984444731456 * linear_interpolation(center_prev->second.pixel[2*j+5] - center_prev->second.pixel[2*j+4], center_next->second.pixel[2*j+5] - center_next->second.pixel[2*j+4], weight) * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside / sqrt(6. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside - 4.) / center_prev->second.hdr.distance;

		tidal[4] = tidal[2] - 1.5 * linear_interpolation(center_prev->second.pixel[(j+1)%4] - center_prev->second.pixel[(j+3)%4], center_next->second.pixel[(j+1)%4] - center_next->second.pixel[(j+3)%4], weight) * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside / sqrt(6. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside - 1.) / center_prev->second.hdr.distance;

		tidal[2] = tidal[2] - 0.5 * tidal[4];

		//tidal[4] *= 3. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside / sqrt(9. + 6.82842712474619 * (3. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside - 2.)) / center_prev->second.hdr.distance;

		tidal[4] = cth * (linear_interpolation(center_prev->second.pixel[2*j+5] - center_prev->second.pixel[2*j+4], center_next->second.pixel[2*j+5] - center_next->second.pixel[2*j+4], weight) / ppp / 0.7071067811865475244 - linear_interpolation(center_prev->second.pixel[(j+1)%4] - center_prev->second.pixel[(j+3)%4], center_next->second.pixel[(j+1)%4] - center_next->second.pixel[(j+3)%4], weight) / sth / 1.847759065022573512) / ppp / center_prev->second.hdr.distance / center_prev->second.hdr.distance;

		if (d > 0 && inner_prev->second.hdr.distance > 0 && inner_prev->second.hdr.Nside == center_prev->second.hdr.Nside && inner_prev->second.hdr.Npix >= 12 && inner_next->second.hdr.Npix >= 12)
		{
			if (inner_prev->second.pixel[(j+1)%4] < -1e30 || inner_prev->second.pixel[(j+3)%4] < -1e30 || inner_next->second.pixel[(j+1)%4] < -1e30 || inner_next->second.pixel[(j+3)%4] < -1e30) return false;
			if (inner_prev->second.pixel[2*j+4] < -1e30 || inner_prev->second.pixel[2*j+5] < -1e30 || inner_next->second.pixel[2*j+4] < -1e30 || inner_next->second.pixel[2*j+5] < -1e30) return false;

			tidal[1] += 0.306186217847897 * linear_interpolation(inner_prev->second.pixel[2*j+4] + inner_prev->second.pixel[2*j+5] - inner_prev->second.pixel[(j+1)%4] - inner_prev->second.pixel[(j+3)%4], inner_next->second.pixel[2*j+4] + inner_next->second.pixel[2*j+5] - inner_next->second.pixel[(j+1)%4] - inner_next->second.pixel[(j+3)%4], weight) * (double) inner_prev->second.hdr.Nside / inner_prev->second.hdr.distance;
			tidal[1] /= -dr;

			d = 0.9799222236572824 * linear_interpolation(inner_prev->second.pixel[2*j+5] - inner_prev->second.pixel[2*j+4], inner_next->second.pixel[2*j+5] - inner_next->second.pixel[2*j+4], weight) / sqrt(6. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside - 4.) + 0.75 * linear_interpolation(inner_prev->second.pixel[(j+1)%4] - inner_prev->second.pixel[(j+3)%4], inner_next->second.pixel[(j+1)%4] - inner_next->second.pixel[(j+3)%4], weight) / sqrt(6. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside - 1.);
			d *= (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside / inner_prev->second.hdr.distance;

			tidal[2] = (tidal[2] - d) / dr;
		}
		else if (outer_prev->second.hdr.Nside == center_prev->second.hdr.Nside)
		{
			if (outer_prev->second.hdr.Npix < 12 || outer_next->second.hdr.Npix < 12) return false;
			if (outer_prev->second.pixel[(j+1)%4] < -1e30 || outer_prev->second.pixel[(j+3)%4] < -1e30 || outer_next->second.pixel[(j+1)%4] < -1e30 || outer_next->second.pixel[(j+3)%4] < -1e30) return false;
			if (outer_prev->second.pixel[2*j+4] < -1e30 || outer_prev->second.pixel[2*j+5] < -1e30 || outer_next->second.pixel[2*j+4] < -1e30 || outer_next->second.pixel[2*j+5] < -1e30) return false;

			tidal[1] += 0.306186217847897 * linear_interpolation(outer_prev->second.pixel[2*j+4] + outer_prev->second.pixel[2*j+5] - outer_prev->second.pixel[(j+1)%4] - outer_prev->second.pixel[(j+3)%4], outer_next->second.pixel[2*j+4] + outer_next->second.pixel[2*j+5] - outer_next->second.pixel[(j+1)%4] - outer_next->second.pixel[(j+3)%4], weight) * (double) outer_prev->second.hdr.Nside / outer_prev->second.hdr.distance;
			tidal[1] /= dr;

			d = 0.9799222236572824 * linear_interpolation(outer_prev->second.pixel[2*j+5] - outer_prev->second.pixel[2*j+4], outer_next->second.pixel[2*j+5] - outer_next->second.pixel[2*j+4], weight) / sqrt(6. * (double) outer_prev->second.hdr.Nside * (double) outer_prev->second.hdr.Nside - 4.) + 0.75 * linear_interpolation(outer_prev->second.pixel[(j+1)%4] - outer_prev->second.pixel[(j+3)%4], outer_next->second.pixel[(j+1)%4] - outer_next->second.pixel[(j+3)%4], weight) / sqrt(6. * (double) outer_prev->second.hdr.Nside * (double) outer_prev->second.hdr.Nside - 1.);
			d *= (double) outer_prev->second.hdr.Nside * (double) outer_prev->second.hdr.Nside / outer_prev->second.hdr.distance;

			tidal[2] = (d - tidal[2]) / dr;
		}
		else
		{
			if (outer_prev->second.hdr.Nside != 2 * center_prev->second.hdr.Nside || outer_prev->second.hdr.Npix < 24 || outer_next->second.hdr.Npix < 24) return false;
			if (outer_prev->second.pixel[2*j+4] < -1e30 || outer_prev->second.pixel[2*j+5] < -1e30 || outer_next->second.pixel[2*j+4] < -1e30 || outer_next->second.pixel[2*j+5] < -1e30) return false;
			if (outer_prev->second.pixel[3*j+13] < -1e30 || outer_next->second.pixel[3*j+13] < -1e30) return false;

			tidal[1] += 1.73205080756888 * linear_interpolation(outer_prev->second.pixel[3*j+13] - outer_prev->second.pixel[j], outer_next->second.pixel[3*j+13] - outer_next->second.pixel[j], weight) * (double) outer_prev->second.hdr.Nside * (double) outer_prev->second.hdr.Nside / sqrt(20. * (double) outer_prev->second.hdr.Nside * (double) outer_prev->second.hdr.Nside - 6. - sqrt((144. * (double) outer_prev->second.hdr.Nside * (double) outer_prev->second.hdr.Nside - 240.) * (double) outer_prev->second.hdr.Nside * (double) outer_prev->second.hdr.Nside - 36.)) / outer_prev->second.hdr.distance;
			tidal[1] /= dr;

			d = 1.959844447314565 * linear_interpolation(outer_prev->second.pixel[2*j+5] - outer_prev->second.pixel[2*j+4], outer_next->second.pixel[2*j+5] - outer_next->second.pixel[2*j+4], weight) * (double) outer_prev->second.hdr.Nside * (double) outer_prev->second.hdr.Nside / sqrt(6. * (double) outer_prev->second.hdr.Nside * (double) outer_prev->second.hdr.Nside - 4.) / outer_prev->second.hdr.distance;

			tidal[2] = (d - tidal[2]) / dr;
		}
	}
	else if (j < (int64_t) center_prev->second.hdr.Nside * ((int64_t) center_prev->second.hdr.Nside - 1) * 2) // in north polar cap
	{
		//if (debug) cout << " in north polar cap: ";
		ring = (1 + (int64_t) sqrt(1.5 + 2 * (double) j)) / 2;
			
		if ((int64_t) center_prev->second.hdr.Npix < 2 * (ring+1) * (ring+2)) return false;

		i = j - 2 * ring * (ring-1);

		m = (i * (ring+1)) / ring;

		cth = (1.5 + m) - ((double) (ring+1) * (0.5 + i) / (double) ring);

		k = m + 2 * ring * (ring+1);
		l = k + 1;
		n = (m > 4*ring + 1) ? 2 * ring * (ring+1) : l + 1;
		m = (m > 0) ? k - 1 : 2 * (ring+1) * (ring+2) - 1;

		if (center_prev->second.pixel[k] < -1e30 || center_prev->second.pixel[l] < -1e30 || center_next->second.pixel[k] < -1e30 || center_next->second.pixel[l] < -1e30) return false;
		if (center_prev->second.pixel[m] < -1e30 || center_prev->second.pixel[n] < -1e30 || center_next->second.pixel[m] < -1e30 || center_next->second.pixel[n] < -1e30) return false;

		/*if (d > 0 && outer_prev->second.hdr.Nside == center_prev->second.hdr.Nside)
		{
			if (outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[l] < -1e30 || outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[l] < -1e30) return false;
			if (outer_prev->second.pixel[m] < -1e30 || outer_prev->second.pixel[n] < -1e30 || outer_next->second.pixel[m] < -1e30 || outer_next->second.pixel[n] < -1e30) return false;

			result[1] = linear_interpolation(cth * outer_prev->second.pixel[k] + (1.-cth) * outer_prev->second.pixel[l], cth * outer_next->second.pixel[k] + (1.-cth) * outer_next->second.pixel[l], weight);
			tmp2 = linear_interpolation(outer_prev->second.pixel[k] - outer_prev->second.pixel[l], outer_next->second.pixel[k] - outer_next->second.pixel[l], weight) * (ring+1);

			ppp2 = linear_interpolation(3 * outer_prev->second.pixel[k] - 3 * outer_prev->second.pixel[l] + outer_prev->second.pixel[n] - outer_prev->second.pixel[m], 3 * outer_next->second.pixel[k] - 3 * outer_next->second.pixel[l] + outer_next->second.pixel[n] - outer_next->second.pixel[m], weight);
			mpp2 = linear_interpolation(outer_prev->second.pixel[k] + outer_prev->second.pixel[l] - outer_prev->second.pixel[n] - outer_prev->second.pixel[m], outer_next->second.pixel[k] + outer_next->second.pixel[l] - outer_next->second.pixel[n] - outer_next->second.pixel[m], weight) - (1. - 2. * cth) * ppp2;
			mmp2 = -tmp2;
		}
		else if (d < 0 && inner_prev->second.hdr.Nside == center_prev->second.hdr.Nside)
		{
			if (inner_prev->second.pixel[k] < -1e30 || inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[k] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;
			if (inner_prev->second.pixel[m] < -1e30 || inner_prev->second.pixel[n] < -1e30 || inner_next->second.pixel[m] < -1e30 || inner_next->second.pixel[n] < -1e30) return false;

			result[1] = linear_interpolation(cth * inner_prev->second.pixel[k] + (1.-cth) * inner_prev->second.pixel[l], cth * inner_next->second.pixel[k] + (1.-cth) * inner_next->second.pixel[l], weight);
			tmp2 = linear_interpolation(inner_prev->second.pixel[k] - inner_prev->second.pixel[l], inner_next->second.pixel[k] - inner_next->second.pixel[l], weight) * (ring+1);

			ppp2 = linear_interpolation(3 * inner_prev->second.pixel[k] - 3 * inner_prev->second.pixel[l] + inner_prev->second.pixel[n] - inner_prev->second.pixel[m], 3 * inner_next->second.pixel[k] - 3 * inner_next->second.pixel[l] + inner_next->second.pixel[n] - inner_next->second.pixel[m], weight);
			mpp2 = linear_interpolation(inner_prev->second.pixel[k] + inner_prev->second.pixel[l] - inner_prev->second.pixel[n] - inner_prev->second.pixel[m], inner_next->second.pixel[k] + inner_next->second.pixel[l] - inner_next->second.pixel[n] - inner_next->second.pixel[m], weight) - (1. - 2. * cth) * ppp2;
			mmp2 = -tmp2;
		}*/

		//if (debug) cout << " in north polar cap: j = " << j << "; k+ = " << k << ", l+ = " << l << ", m+ = " << m << ", n+ = " << n << ", w+ = " << cth;

		tidal[1] = linear_interpolation(cth * center_prev->second.pixel[k] + (1.-cth) * center_prev->second.pixel[l], cth * center_next->second.pixel[k] + (1.-cth) * center_next->second.pixel[l], weight);
		tmp = linear_interpolation(center_prev->second.pixel[k] - center_prev->second.pixel[l], center_next->second.pixel[k] - center_next->second.pixel[l], weight) * (ring+1);

		ppp = linear_interpolation(3 * center_prev->second.pixel[k] - 3 * center_prev->second.pixel[l] + center_prev->second.pixel[n] - center_prev->second.pixel[m], 3 * center_next->second.pixel[k] - 3 * center_next->second.pixel[l] + center_next->second.pixel[n] - center_next->second.pixel[m], weight);
		mpp = linear_interpolation(center_prev->second.pixel[k] + center_prev->second.pixel[l] - center_prev->second.pixel[n] - center_prev->second.pixel[m], center_next->second.pixel[k] + center_next->second.pixel[l] - center_next->second.pixel[n] - center_next->second.pixel[m], weight) - (1. - 2. * cth) * ppp;
		mmp = -tmp;

		if (i % ring == 0)
			m = ((4 + i / ring) * (ring-1) - 1) % (4 * (ring-1));
		else
			m = (i * (ring-1)) / ring;

		sth = (1.5 + m) - ((double) (ring-1) * (0.5 + i) / (double) ring);

		if (i == 0) sth -= 4 * (ring-1);

		l = ((m+1) % (4 * (ring-1))) + 2 * (ring-1) * (ring-2);
		k = m + 2 * (ring-1) * (ring-2);

		if (center_prev->second.pixel[k] < -1e30 || center_prev->second.pixel[l] < -1e30 || center_next->second.pixel[k] < -1e30 || center_next->second.pixel[l] < -1e30) return false;

		/*if (d > 0 && outer_prev->second.hdr.Nside == center_prev->second.hdr.Nside)
		{
			if (outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[l] < -1e30 || outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[l] < -1e30) return false;

			result[2] = linear_interpolation(sth * outer_prev->second.pixel[k] + (1.-sth) * outer_prev->second.pixel[l], sth * outer_next->second.pixel[k] + (1.-sth) * outer_next->second.pixel[l], weight);
			tmp2 -= linear_interpolation(outer_prev->second.pixel[k] - outer_prev->second.pixel[l], outer_next->second.pixel[k] - outer_next->second.pixel[l], weight) * (ring-1);
			mmp2 -= linear_interpolation(outer_prev->second.pixel[k] - outer_prev->second.pixel[l], outer_next->second.pixel[k] - outer_next->second.pixel[l], weight) * (ring-1);
			tmp2 /= 0.5 * M_PI;
		}
		else if (d < 0 && inner_prev->second.hdr.Nside == center_prev->second.hdr.Nside)
		{
			if (inner_prev->second.pixel[k] < -1e30 || inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[k] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

			result[2] = linear_interpolation(sth * inner_prev->second.pixel[k] + (1.-sth) * inner_prev->second.pixel[l], sth * inner_next->second.pixel[k] + (1.-sth) * inner_next->second.pixel[l], weight);
			tmp2 -= linear_interpolation(inner_prev->second.pixel[k] - inner_prev->second.pixel[l], inner_next->second.pixel[k] - inner_next->second.pixel[l], weight) * (ring-1);
			mmp2 -= linear_interpolation(inner_prev->second.pixel[k] - inner_prev->second.pixel[l], inner_next->second.pixel[k] - inner_next->second.pixel[l], weight) * (ring-1);
			tmp2 /= 0.5 * M_PI;
		}*/

		//if (debug) cout << "; k- = " << k << ", l- = " << l << ", w- = " << sth << endl;

		tidal[2] = linear_interpolation(sth * center_prev->second.pixel[k] + (1.-sth) * center_prev->second.pixel[l], sth * center_next->second.pixel[k] + (1.-sth) * center_next->second.pixel[l], weight);
		tmp -= linear_interpolation(center_prev->second.pixel[k] - center_prev->second.pixel[l], center_next->second.pixel[k] - center_next->second.pixel[l], weight) * (ring-1);
		mmp -= linear_interpolation(center_prev->second.pixel[k] - center_prev->second.pixel[l], center_next->second.pixel[k] - center_next->second.pixel[l], weight) * (ring-1);
		tmp /= 0.5 * M_PI;

		k = (i == 0) ? j-1+4*ring : j-1;
		l = (i == 4*ring-1) ? j+1-4*ring : j+1;

		if (center_prev->second.pixel[k] < -1e30 || center_prev->second.pixel[l] < -1e30 || center_next->second.pixel[k] < -1e30 || center_next->second.pixel[l] < -1e30) return false;

		/*if (d > 0 && outer_prev->second.hdr.Nside == center_prev->second.hdr.Nside)
		{
			if (outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[l] < -1e30 || outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[l] < -1e30) return false;

			sph2 = (linear_interpolation(outer_prev->second.pixel[k] + outer_prev->second.pixel[l], outer_next->second.pixel[k] + outer_next->second.pixel[l], weight) - 2. * result[4]) * 4. * ring * ring;
			mmp2 += linear_interpolation(outer_prev->second.pixel[k] - outer_prev->second.pixel[l], outer_next->second.pixel[k] - outer_next->second.pixel[l], weight) * ring;
		}
		else if (d < 0 && inner_prev->second.hdr.Nside == center_prev->second.hdr.Nside)
		{
			if (inner_prev->second.pixel[k] < -1e30 || inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[k] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

			sph2 = (linear_interpolation(inner_prev->second.pixel[k] + inner_prev->second.pixel[l], inner_next->second.pixel[k] + inner_next->second.pixel[l], weight) - 2. * result[4]) * 4. * ring * ring;
			mmp2 += linear_interpolation(inner_prev->second.pixel[k] - inner_prev->second.pixel[l], inner_next->second.pixel[k] - inner_next->second.pixel[l], weight) * ring;
		}
		mmp2 -= 0.125 * ((2.*sth-1.) / (ring-1) + (2.*cth-1.) / (ring+1)) * sph2;

		mpp2 += 0.5 * sph2 / (ring+1) / (ring+1);

		mmp2 -= 0.25 * ((2.*sth-1.)*(2*ring-1)*(ring+1)/(2*ring+1)/(ring-1) - (2.*cth-1.)) * mpp2 * (ring+1);
		mmp2 += ppp2 * 0.5 * (ring+1) * ((cth-1.)*cth + (2. / 3. / ring / ring / (ring-1) / (ring-1)) + ((sth-1.)*sth*(ring+1)*(ring+1) - 2.) / (ring-1) / (ring-1));
		mmp2 /= 0.5 * M_PI;

		tmp2 += ((0.5-sth)/(ring-1) - (0.5-cth)/(ring+1)) * 0.5 * sph2 / M_PI;
		tmp2 += mmp2 * ring / (0.25 + ring * ring);
		tmp2 /= 4 * ring / 3. / center_prev->second.hdr.Nside / center_prev->second.hdr.Nside;
		tmp2 *= 1. + 0.25 / ring / ring;

		cph2 = (2*ring+1)*(result[2]-result[4]) + (2*ring-1)*(result[1]-result[4]) - 0.125 * ((2*ring+1)*sth*(1.-sth)/(ring-1)/(ring-1) + (2*ring-1)*cth*(1.-cth)/(ring+1)/(ring+1)) * sph2;
		cph2 += ppp2 * ((2*ring+1)*sth*(sth-1.)*(sth-0.5)*(ring+1)*(ring+1)*(ring+1)/(ring-1)/(ring-1)/(ring-1) + (2*ring-1)*cth*(cth-1.)*(cth-0.5)) / 3.;
		cph2 -= mpp2 * (ring-0.5) * (sth*(1.-sth)*(ring+1)*(ring+1)/(ring-1)/(ring-1) - cth*(1.-cth)) / 2.;
		cph2 /= 2. * (2*ring-1) * (2*ring+1) * ring / 9. / center_prev->second.hdr.Nside / center_prev->second.hdr.Nside / center_prev->second.hdr.Nside / center_prev->second.hdr.Nside;

		sth2 = result[2] - result[1] + 0.125 * (cth*(1.-cth)/(ring+1)/(ring+1) - sth*(1.-sth)/(ring-1)/(ring-1)) * sph2 + ppp2 * (sth*(sth-1.)*(sth-0.5)*(ring+1)*(ring+1)*(ring+1)/(ring-1)/(ring-1)/(ring-1) - cth*(cth-1.)*(cth-0.5)) / 3. - mpp2 * (sth*(1.-sth)*(ring+1)*(ring+1)*(ring-0.5)/(ring-1)/(ring-1)/(2*ring+1) + cth*(1.-cth)) / 2.;
		sth2 /= 4 * ring / 3. / center_prev->second.hdr.Nside / center_prev->second.hdr.Nside;
		sth2 += cph2 / 3. / center_prev->second.hdr.Nside / center_prev->second.hdr.Nside; ///////*/

		sph = (linear_interpolation(center_prev->second.pixel[k] + center_prev->second.pixel[l], center_next->second.pixel[k] + center_next->second.pixel[l], weight) - 2. * tidal[4]) * 4. * ring * ring;
		mmp += linear_interpolation(center_prev->second.pixel[k] - center_prev->second.pixel[l], center_next->second.pixel[k] - center_next->second.pixel[l], weight) * ring;
		mmp -= 0.125 * ((2.*sth-1.) / (ring-1) + (2.*cth-1.) / (ring+1)) * sph;

		mpp += 0.5 * sph / (ring+1) / (ring+1);

		mmp -= 0.25 * ((2.*sth-1.)*(2*ring-1)*(ring+1)/(2*ring+1)/(ring-1) - (2.*cth-1.)) * mpp * (ring+1);
		mmp += ppp * 0.5 * (ring+1) * ((cth-1.)*cth + (2. / 3. / ring / ring / (ring-1) / (ring-1)) + ((sth-1.)*sth*(ring+1)*(ring+1) - 2.) / (ring-1) / (ring-1));
		mmp /= 0.5 * M_PI;

		tmp += ((0.5-sth)/(ring-1) - (0.5-cth)/(ring+1)) * 0.5 * sph / M_PI;
		tmp += mmp * ring / (0.25 + ring * ring);
		tmp /= 4 * ring / 3. / (double) center_prev->second.hdr.Nside / (double) center_prev->second.hdr.Nside;
		tmp *= 1. + 0.25 /  ring / ring;

		cph = (2*ring+1)*(tidal[2]-tidal[4]) + (2*ring-1)*(tidal[1]-tidal[4]) - 0.125 * ((2*ring+1)*sth*(1.-sth)/(ring-1)/(ring-1) + (2*ring-1)*cth*(1.-cth)/(ring+1)/(ring+1)) * sph;
		cph += ppp * ((2*ring+1)*sth*(sth-1.)*(sth-0.5)*(ring+1)*(ring+1)*(ring+1)/(ring-1)/(ring-1)/(ring-1) + (2*ring-1)*cth*(cth-1.)*(cth-0.5)) / 3.;
		cph -= mpp * (ring-0.5) * (sth*(1.-sth)*(ring+1)*(ring+1)/(ring-1)/(ring-1) - cth*(1.-cth)) / 2.;
		cph /= 2. * (2*ring-1) * (2*ring+1) * ring / 9. / (double) center_prev->second.hdr.Nside / (double) center_prev->second.hdr.Nside / (double) center_prev->second.hdr.Nside / (double) center_prev->second.hdr.Nside;

		sth = tidal[2] - tidal[1] + 0.125 * (cth*(1.-cth)/(ring+1)/(ring+1) - sth*(1.-sth)/(ring-1)/(ring-1)) * sph + ppp * (sth*(sth-1.)*(sth-0.5)*(ring+1)*(ring+1)*(ring+1)/(ring-1)/(ring-1)/(ring-1) - cth*(cth-1.)*(cth-0.5)) / 3. - mpp * (sth*(1.-sth)*(ring+1)*(ring+1)*(ring-0.5)/(ring-1)/(ring-1)/(2*ring+1) + cth*(1.-cth)) / 2.;
		sth /= 4 * ring / 3. / (double) center_prev->second.hdr.Nside / (double) center_prev->second.hdr.Nside;
		sth += cph / 3. / (double) center_prev->second.hdr.Nside / (double) center_prev->second.hdr.Nside;

		cth = 1. - ring*ring / (3. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside);

		tidal[3] += ((1.-cth*cth)*cph - cth*sth) / center_prev->second.hdr.distance / center_prev->second.hdr.distance;

		tidal[5] += (sph / M_PI / M_PI / (1.-cth*cth) - cth*sth)  / center_prev->second.hdr.distance / center_prev->second.hdr.distance;

		tidal[1] -= tidal[2];

		tidal[2] = linear_interpolation(center_prev->second.pixel[l] - center_prev->second.pixel[k], center_next->second.pixel[l] - center_next->second.pixel[k], weight);

		tidal[4] = (cth*(tidal[2]*ring - ppp*(ring+1)*(ring+1)*(ring+1)/ring/ring/3.)/M_PI/(cth*cth-1.) - tmp) / center_prev->second.hdr.distance / center_prev->second.hdr.distance;

		tidal[1] /= center_prev->second.hdr.distance;

		/*if (d > 0 && outer_prev->second.hdr.Nside == center_prev->second.hdr.Nside)
		{
			result[3] += ((1.-cth*cth)*cph2 - cth*sth2) / outer_prev->second.hdr.distance / outer_prev->second.hdr.distance;
			result[5] += (sph2 / M_PI / M_PI / (1.-cth*cth) - cth*sth2) / outer_prev->second.hdr.distance / outer_prev->second.hdr.distance;
			result[4] = (cth*(linear_interpolation(outer_prev->second.pixel[l] - outer_prev->second.pixel[k], outer_next->second.pixel[l] - outer_next->second.pixel[k], weight)*ring - ppp2*(ring+1)*(ring+1)*(ring+1)/ring/ring/3.)/M_PI/(cth*cth-1.) - tmp2) / outer_prev->second.hdr.distance / outer_prev->second.hdr.distance;

			tidal[3] = tidal[3] * (1. - wr) + result[3] * wr;
			tidal[4] = tidal[4] * (1. - wr) + result[4] * wr;
			tidal[5] = tidal[5] * (1. - wr) + result[5] * wr;
		}
		else if (d < 0 && inner_prev->second.hdr.Nside == center_prev->second.hdr.Nside)
		{
			result[3] += ((1.-cth*cth)*cph2 - cth*sth2) / inner_prev->second.hdr.distance / inner_prev->second.hdr.distance;
			result[5] += (sph2 / M_PI / M_PI / (1.-cth*cth) - cth*sth2) / inner_prev->second.hdr.distance / inner_prev->second.hdr.distance;
			result[4] = (cth*(linear_interpolation(inner_prev->second.pixel[l] - inner_prev->second.pixel[k], inner_next->second.pixel[l] - inner_next->second.pixel[k], weight)*ring - ppp2*(ring+1)*(ring+1)*(ring+1)/ring/ring/3.)/M_PI/(cth*cth-1.) - tmp2) / inner_prev->second.hdr.distance / inner_prev->second.hdr.distance;

			tidal[3] = tidal[3] * wr + result[3] * (1. - wr);
			tidal[4] = tidal[4] * wr + result[4] * (1. - wr);
			tidal[5] = tidal[5] * wr + result[5] * (1. - wr);
		}*/

		sth = sqrt(1. - cth*cth);

		sph = sin(0.5 * M_PI / ring);

		tidal[2] *= 0.5 / sph / sth / center_prev->second.hdr.distance;

		cph = 3. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside / sqrt(12. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside * (ring*ring + 1) - 2. * (ring*ring - 1) * (ring*ring - 1 + sqrt((6. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside - 1.) * (6. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside - 1.) - ring*ring * (12. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside - ring*ring + 2))));

		if (d > 0 && inner_prev->second.hdr.distance > 0 && inner_prev->second.hdr.Nside == center_prev->second.hdr.Nside && inner_prev->second.hdr.Npix >= 2 * (ring+1) * (ring+2) && inner_next->second.hdr.Npix >= 2 * (ring+1) * (ring+2))
		{
			if (inner_prev->second.pixel[k] < -1e30 || inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[k] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

			tidal[2] = (tidal[2] - linear_interpolation(inner_prev->second.pixel[l] - inner_prev->second.pixel[k], inner_next->second.pixel[l] - inner_next->second.pixel[k], weight) * 0.5 / sph / sth / inner_prev->second.hdr.distance) / dr;

			k = (i * (ring+1)) / ring;

			d = (1.5 + k) - ((double) (ring+1) * (0.5 + i) / (double) ring);

			k += 2 * ring * (ring+1);
			l = k + 1;

			if (inner_prev->second.pixel[k] < -1e30 || inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[k] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

			tidal[1] -= linear_interpolation(d * inner_prev->second.pixel[k] + (1.-d) * inner_prev->second.pixel[l], d * inner_next->second.pixel[k] + (1.-d) * inner_next->second.pixel[l], weight) / inner_prev->second.hdr.distance;

			if (i % ring == 0)
				k = ((4 + i / ring) * (ring-1) - 1) % (4 * (ring-1));
			else
				k = (i * (ring-1)) / ring;

			d = (1.5 + k) - ((double) (ring-1) * (0.5 + i) / (double) ring);

			if (i == 0) d -= 4 * (ring-1);

			l = ((k+1) % (4 * (ring-1))) + 2 * (ring-1) * (ring-2);
			k += 2 * (ring-1) * (ring-2);

			if (inner_prev->second.pixel[k] < -1e30 || inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[k] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

			tidal[1] += linear_interpolation(d * inner_prev->second.pixel[k] + (1.-d) * inner_prev->second.pixel[l], d * inner_next->second.pixel[k] + (1.-d) * inner_next->second.pixel[l], weight) / inner_prev->second.hdr.distance;

			tidal[1] *= cph / dr;
		}
		else if (outer_prev->second.hdr.Nside == center_prev->second.hdr.Nside)
		{
			if (outer_prev->second.hdr.Npix < 2 * (ring+1) * (ring+2) || outer_next->second.hdr.Npix < 2 * (ring+1) * (ring+2)) return false;
			if (outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[l] < -1e30 || outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[l] < -1e30) return false;

			tidal[2] = (linear_interpolation(outer_prev->second.pixel[l] - outer_prev->second.pixel[k], outer_next->second.pixel[l] - outer_next->second.pixel[k], weight) * 0.5 / sph / sth / outer_prev->second.hdr.distance - tidal[2]) / dr;

			k = (i * (ring+1)) / ring;

			d = (1.5 + k) - ((double) (ring+1) * (0.5 + i) / (double) ring);

			k += 2 * ring * (ring+1);
			l = k + 1;

			if (outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[l] < -1e30 || outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[l] < -1e30) return false;

			tidal[1] -= linear_interpolation(d * outer_prev->second.pixel[k] + (1.-d) * outer_prev->second.pixel[l], d * outer_next->second.pixel[k] + (1.-d) * outer_next->second.pixel[l], weight) / outer_prev->second.hdr.distance;

			if (i % ring == 0)
				k = ((4 + i / ring) * (ring-1) - 1) % (4 * (ring-1));
			else
				k = (i * (ring-1)) / ring;

			d = (1.5 + k) - ((double) (ring-1) * (0.5 + i) / (double) ring);

			if (i == 0) d -= 4 * (ring-1);

			l = ((k+1) % (4 * (ring-1))) + 2 * (ring-1) * (ring-2);
			k += 2 * (ring-1) * (ring-2);

			if (outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[l] < -1e30 || outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[l] < -1e30) return false;

			tidal[1] += linear_interpolation(d * outer_prev->second.pixel[k] + (1.-d) * outer_prev->second.pixel[l], d * outer_next->second.pixel[k] + (1.-d) * outer_next->second.pixel[l], weight) / outer_prev->second.hdr.distance;

			tidal[1] *= -cph / dr;
		}
		else
		{
			if (outer_prev->second.hdr.Nside != 2 * center_prev->second.hdr.Nside || outer_prev->second.hdr.Npix < 2 * (2*ring+1) * (2*ring+2) || outer_next->second.hdr.Npix < 2 * (2*ring+1) * (2*ring+2)) return false;
			tidal[1] *= cph;

			k = 2*i + 4*ring*(2*ring - 1);
			l = k+1;

			if (outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[l] < -1e30 || outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[l] < -1e30) return false;

			cph = linear_interpolation(outer_prev->second.pixel[l] - outer_prev->second.pixel[k], outer_next->second.pixel[l] - outer_next->second.pixel[k], weight);

			tidal[2] = (cph * 0.5 / sin(0.125 * M_PI / ring)  / sth / outer_prev->second.hdr.distance - tidal[2]) / dr;

			k = (int64_t) floor((2*ring - 1)*(i + 0.5) / ring) + 4*(2*ring - 1)*(ring-1);
			l = (int64_t) floor((2*ring + 1)*(i + 0.5) / ring) + 4*ring*(2*ring + 1);

			if (outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[l] < -1e30 || outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[l] < -1e30) return false;

			d = (k + 0.5) / (2*ring-1) - (l + 0.5) / (2*ring + 1);

			tidal[1] = (linear_interpolation(outer_prev->second.pixel[l] - outer_prev->second.pixel[k], outer_next->second.pixel[l] - outer_next->second.pixel[k], weight) + d * cph) * 3. * (double) outer_prev->second.hdr.Nside * (double) outer_prev->second.hdr.Nside / sqrt(12. * (double) outer_prev->second.hdr.Nside * (double) outer_prev->second.hdr.Nside * (4*ring*ring + 1) - 2. * (4*ring*ring - 1) * (4*ring*ring - 1 + sqrt((6. * (double) outer_prev->second.hdr.Nside * (double) outer_prev->second.hdr.Nside - 1.) * (6. * (double) outer_prev->second.hdr.Nside * (double) outer_prev->second.hdr.Nside - 1.) - 4*ring*ring * (12. * (double) outer_prev->second.hdr.Nside * (double) outer_prev->second.hdr.Nside - 4*ring*ring + 2)))) / outer_prev->second.hdr.distance - tidal[1];

			tidal[1] /= dr;
		}

		cph = cos((i + 0.5) * M_PI / (2*ring));
		sph = (i < 2*ring) ? sqrt(1. - cph*cph) : -sqrt(1. - cph*cph);
	}
	else return false;

	// rotate

	result[0] = cph*cph * (tidal[0]*sth*sth + tidal[3]*cth*cth + 2*tidal[1]*cth*sth) + sph * (tidal[5]*sph - 2*cph * (tidal[2]*sth + tidal[4]*cth));
	result[1] = cph*sph * (tidal[0]*sth*sth + tidal[3]*cth*cth + 2*tidal[1]*cth*sth - tidal[5]) + (cph*cph - sph*sph) * (tidal[2]*sth + tidal[4]*cth);
	result[2] = cph * ((tidal[0] - tidal[3])*cth*sth + tidal[1]*(cth*cth - sth*sth)) + sph * (tidal[5]*sth - tidal[2]*cth);
	result[3] = sph*sph * (tidal[0]*sth*sth + tidal[3]*cth*cth + 2*tidal[1]*cth*sth) + cph * (tidal[5]*cph + 2*sph * (tidal[2]*sth + tidal[4]*cth));
	result[4] = sph * ((tidal[0] - tidal[3])*cth*sth + tidal[1]*(cth*cth - sth*sth)) - cph * (tidal[5]*sth - tidal[2]*cth);
	result[5] = tidal[0]*cth*cth - 2*tidal[1]*cth*sth + tidal[3]*sth*sth;

	return true;
}

bool tidal_matrix_old(metric_container * field_prev, metric_container * field_next, const double weight, const double dist, const double dr, double * v, double * result, int debug)
{
	std::map<int,metric_data>::iterator inner_prev = get_iterator(field_prev, dist, -1);
	std::map<int,metric_data>::iterator inner_next = get_iterator(field_next, dist, -1);
	std::map<int,metric_data>::iterator outer_prev = get_iterator(field_prev, dist);
	std::map<int,metric_data>::iterator outer_next = get_iterator(field_next, dist);
	std::map<int,metric_data>::iterator center_prev;
	std::map<int,metric_data>::iterator center_next;
	double d, cth, sth, cph, sph, p1, p2, h, corr1, corr2;
	int64_t i, j, k, l, m, ring, Npix;
	double tidal[6];

	if (inner_prev == field_prev->healpix_data.end() || inner_next == field_next->healpix_data.end() || outer_prev == field_prev->healpix_data.end() || outer_next == field_next->healpix_data.end()) return false;

	if (inner_prev->second.hdr.distance != inner_next->second.hdr.distance || inner_prev->second.hdr.Nside != inner_next->second.hdr.Nside) return false;

	if (outer_prev->second.hdr.distance != outer_next->second.hdr.distance || outer_prev->second.hdr.Nside != outer_next->second.hdr.Nside) return false;

	if(modf(dist / dr, &d) >= 0.5)
	{
		center_prev = outer_prev;
		center_next = outer_next;
		
		outer_prev++; outer_next++;

		if (outer_prev == field_prev->healpix_data.end() || outer_next == field_next->healpix_data.end()) return false;

		if (outer_prev->second.hdr.distance != outer_next->second.hdr.distance || outer_prev->second.hdr.Nside != outer_next->second.hdr.Nside) return false;

		d = -1.;
	}
	else if (d > 0)
	{
		if (inner_prev->first == 0 || inner_next->first == 0) return false;

		center_prev = inner_prev;
		center_next = inner_next;

		if (inner_prev == field_prev->healpix_data.begin() || inner_next == field_next->healpix_data.begin()) return false;

		inner_prev--; inner_next--;

		if (inner_prev->second.hdr.distance != inner_next->second.hdr.distance || inner_prev->second.hdr.Nside != inner_next->second.hdr.Nside) return false;
	}
	else // at vertex
	{
		if (debug) cout << " at vertex" << endl;
		if (inner_prev->second.hdr.distance != 0) return false;

		Npix = nside2npix64((int64_t) outer_prev->second.hdr.Nside);

		tidal[1] = 2. * linear_interpolation(inner_prev->second.pixel[0], inner_next->second.pixel[0], weight);

		i = Npix/2;
		j = i-1;
		k = i + 4*outer_prev->second.hdr.Nside;
		l = i - 4*outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[l] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[i] < -1e30 || outer_prev->second.pixel[k] < -1e30) return false;
		if (outer_next->second.pixel[l] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[i] < -1e30 || outer_next->second.pixel[k] < -1e30) return false;

		tidal[0] = linear_interpolation(outer_prev->second.pixel[l] + outer_prev->second.pixel[j] + outer_prev->second.pixel[i] + outer_prev->second.pixel[k], outer_next->second.pixel[l] + outer_next->second.pixel[j] + outer_next->second.pixel[i] + outer_next->second.pixel[k], weight);

		i -= 2*outer_prev->second.hdr.Nside;
		j += 2*outer_prev->second.hdr.Nside;
		k -= 2*outer_prev->second.hdr.Nside;
		l -= 2*outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[l] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[i] < -1e30 || outer_prev->second.pixel[k] < -1e30) return false;
		if (outer_next->second.pixel[l] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[i] < -1e30 || outer_next->second.pixel[k] < -1e30) return false;

		tidal[0] += linear_interpolation(outer_prev->second.pixel[l] + outer_prev->second.pixel[j] + outer_prev->second.pixel[i] + outer_prev->second.pixel[k], outer_next->second.pixel[l] + outer_next->second.pixel[j] + outer_next->second.pixel[i] + outer_next->second.pixel[k], weight);

		tidal[0] = (0.25 * tidal[0] - tidal[1]) / dr / dr;

		i += outer_prev->second.hdr.Nside;
		j -= outer_prev->second.hdr.Nside;
		k += outer_prev->second.hdr.Nside;
		l += outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[l] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[i] < -1e30 || outer_prev->second.pixel[k] < -1e30) return false;
		if (outer_next->second.pixel[l] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[i] < -1e30 || outer_next->second.pixel[k] < -1e30) return false;

		tidal[3] = linear_interpolation(outer_prev->second.pixel[l] + outer_prev->second.pixel[j] + outer_prev->second.pixel[i] + outer_prev->second.pixel[k], outer_next->second.pixel[l] + outer_next->second.pixel[j] + outer_next->second.pixel[i] + outer_next->second.pixel[k], weight);

		i += 2*outer_prev->second.hdr.Nside;
		j -= 2*outer_prev->second.hdr.Nside;
		k += 2*outer_prev->second.hdr.Nside;
		l += 2*outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[l] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[i] < -1e30 || outer_prev->second.pixel[k] < -1e30) return false;
		if (outer_next->second.pixel[l] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[i] < -1e30 || outer_next->second.pixel[k] < -1e30) return false;

		tidal[3] += linear_interpolation(outer_prev->second.pixel[l] + outer_prev->second.pixel[j] + outer_prev->second.pixel[i] + outer_prev->second.pixel[k], outer_next->second.pixel[l] + outer_next->second.pixel[j] + outer_next->second.pixel[i] + outer_next->second.pixel[k], weight);

		tidal[3] = (0.25 * tidal[3] - tidal[1]) / dr / dr;

		if (outer_prev->second.pixel[0] < -1e30 || outer_prev->second.pixel[1] < -1e30 || outer_prev->second.pixel[2] < -1e30 || outer_prev->second.pixel[3] < -1e30) return false;
		if (outer_next->second.pixel[0] < -1e30 || outer_next->second.pixel[1] < -1e30 || outer_next->second.pixel[2] < -1e30 || outer_next->second.pixel[3] < -1e30) return false;

		tidal[5] = linear_interpolation(outer_prev->second.pixel[0] + outer_prev->second.pixel[1] + outer_prev->second.pixel[2] + outer_prev->second.pixel[3], outer_next->second.pixel[0] + outer_next->second.pixel[1] + outer_next->second.pixel[2] + outer_next->second.pixel[3], weight);

		if (outer_prev->second.pixel[Npix-4] < -1e30 || outer_prev->second.pixel[Npix-3] < -1e30 || outer_prev->second.pixel[Npix-2] < -1e30 || outer_prev->second.pixel[Npix-1] < -1e30) return false;
		if (outer_next->second.pixel[Npix-4] < -1e30 || outer_next->second.pixel[Npix-3] < -1e30 || outer_next->second.pixel[Npix-2] < -1e30 || outer_next->second.pixel[Npix-1] < -1e30) return false;

		tidal[5] += linear_interpolation(outer_prev->second.pixel[Npix-4] + outer_prev->second.pixel[Npix-3] + outer_prev->second.pixel[Npix-2] + outer_prev->second.pixel[Npix-1], outer_next->second.pixel[Npix-4] + outer_next->second.pixel[Npix-3] + outer_next->second.pixel[Npix-2] + outer_next->second.pixel[Npix-1], weight);

		tidal[5] = (0.25 * tidal[5] - tidal[1]) / dr / dr;

		i += outer_prev->second.hdr.Nside/2;
		j = i - 1;
		k = i + 4*outer_prev->second.hdr.Nside;
		l = i - 4*outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[l] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[i] < -1e30 || outer_prev->second.pixel[k] < -1e30) return false;
		if (outer_next->second.pixel[l] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[i] < -1e30 || outer_next->second.pixel[k] < -1e30) return false;

		tidal[1] = linear_interpolation(outer_prev->second.pixel[l] + outer_prev->second.pixel[j] + outer_prev->second.pixel[i] + outer_prev->second.pixel[k], outer_next->second.pixel[l] + outer_next->second.pixel[j] + outer_next->second.pixel[i] + outer_next->second.pixel[k], weight);

		i -= outer_prev->second.hdr.Nside;
		j -= outer_prev->second.hdr.Nside;
		k -= outer_prev->second.hdr.Nside;
		l -= outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[l] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[i] < -1e30 || outer_prev->second.pixel[k] < -1e30) return false;
		if (outer_next->second.pixel[l] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[i] < -1e30 || outer_next->second.pixel[k] < -1e30) return false;

		tidal[1] -= linear_interpolation(outer_prev->second.pixel[l] + outer_prev->second.pixel[j] + outer_prev->second.pixel[i] + outer_prev->second.pixel[k], outer_next->second.pixel[l] + outer_next->second.pixel[j] + outer_next->second.pixel[i] + outer_next->second.pixel[k], weight);

		i -= outer_prev->second.hdr.Nside;
		j -= outer_prev->second.hdr.Nside;
		k -= outer_prev->second.hdr.Nside;
		l -= outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[l] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[i] < -1e30 || outer_prev->second.pixel[k] < -1e30) return false;
		if (outer_next->second.pixel[l] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[i] < -1e30 || outer_next->second.pixel[k] < -1e30) return false;

		tidal[1] += linear_interpolation(outer_prev->second.pixel[l] + outer_prev->second.pixel[j] + outer_prev->second.pixel[i] + outer_prev->second.pixel[k], outer_next->second.pixel[l] + outer_next->second.pixel[j] + outer_next->second.pixel[i] + outer_next->second.pixel[k], weight);

		i -= outer_prev->second.hdr.Nside;
		j -= outer_prev->second.hdr.Nside;
		k -= outer_prev->second.hdr.Nside;
		l -= outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[l] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[i] < -1e30 || outer_prev->second.pixel[k] < -1e30) return false;
		if (outer_next->second.pixel[l] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[i] < -1e30 || outer_next->second.pixel[k] < -1e30) return false;

		tidal[1] -= linear_interpolation(outer_prev->second.pixel[l] + outer_prev->second.pixel[j] + outer_prev->second.pixel[i] + outer_prev->second.pixel[k], outer_next->second.pixel[l] + outer_next->second.pixel[j] + outer_next->second.pixel[i] + outer_next->second.pixel[k], weight);

		tidal[1] /= -8. * dr * dr;

		i = 2 * (int64_t) outer_prev->second.hdr.Nside * (outer_prev->second.hdr.Nside + 1);
		j = i - 1;
		k = i - 4*outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[i] < -1e30) return false;
		if (outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[i] < -1e30) return false;

		tidal[2] = linear_interpolation(outer_prev->second.pixel[k] + outer_prev->second.pixel[j] + outer_prev->second.pixel[i], outer_next->second.pixel[k] + outer_next->second.pixel[j] + outer_next->second.pixel[i], weight);

		i += 2*outer_prev->second.hdr.Nside;
		j -= 2*outer_prev->second.hdr.Nside;
		k += 2*outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[i] < -1e30) return false;
		if (outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[i] < -1e30) return false;

		tidal[2] -= linear_interpolation(outer_prev->second.pixel[j] + outer_prev->second.pixel[k] + outer_prev->second.pixel[i], outer_next->second.pixel[j] + outer_next->second.pixel[k] + outer_next->second.pixel[i], weight);

		i += 8l * (int64_t) outer_prev->second.hdr.Nside * (outer_prev->second.hdr.Nside - 1);
		k = i + 4*outer_prev->second.hdr.Nside;
		j = k - 1;

		if (outer_prev->second.pixel[i] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[k] < -1e30) return false;
		if (outer_next->second.pixel[i] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[k] < -1e30) return false;

		tidal[2] += linear_interpolation(outer_prev->second.pixel[i] + outer_prev->second.pixel[j] + outer_prev->second.pixel[k], outer_next->second.pixel[i] + outer_next->second.pixel[j] + outer_next->second.pixel[k], weight);

		i -= 2*outer_prev->second.hdr.Nside;
		k += 2*outer_prev->second.hdr.Nside;
		j -= 2*outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[i] < -1e30 || outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[j] < -1e30) return false;
		if (outer_next->second.pixel[i] < -1e30 || outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[j] < -1e30) return false;

		tidal[2] -= linear_interpolation(outer_prev->second.pixel[i] + outer_prev->second.pixel[k] + outer_prev->second.pixel[j], outer_next->second.pixel[i] + outer_next->second.pixel[k] + outer_next->second.pixel[j], weight);

		i = (int64_t) outer_prev->second.hdr.Nside * (2 * outer_prev->second.hdr.Nside + 3);
		j = i - 4*outer_prev->second.hdr.Nside;
		k = j - 1;

		if (outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[i] < -1e30) return false;
		if (outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[i] < -1e30) return false;

		tidal[4] = linear_interpolation(outer_prev->second.pixel[k] + outer_prev->second.pixel[j] + outer_prev->second.pixel[i], outer_next->second.pixel[k] + outer_next->second.pixel[j] + outer_next->second.pixel[i], weight);

		i += 2*outer_prev->second.hdr.Nside;
		j += 2*outer_prev->second.hdr.Nside;
		k += 2*outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[i] < -1e30) return false;
		if (outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[i] < -1e30) return false;

		tidal[4] -= linear_interpolation(outer_prev->second.pixel[k] + outer_prev->second.pixel[j] + outer_prev->second.pixel[i], outer_next->second.pixel[k] + outer_next->second.pixel[j] + outer_next->second.pixel[i], weight);

		i += 8l * (int64_t) outer_prev->second.hdr.Nside * (outer_prev->second.hdr.Nside - 1);
		k = i + 4*outer_prev->second.hdr.Nside;
		j = k - 1;

		if (outer_prev->second.pixel[i] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[k] < -1e30) return false;
		if (outer_next->second.pixel[i] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[k] < -1e30) return false;

		tidal[4] += linear_interpolation(outer_prev->second.pixel[i] + outer_prev->second.pixel[j] + outer_prev->second.pixel[k], outer_next->second.pixel[i] + outer_next->second.pixel[j] + outer_next->second.pixel[k], weight);

		i -= 2*outer_prev->second.hdr.Nside;
		j -= 2*outer_prev->second.hdr.Nside;
		k -= 2*outer_prev->second.hdr.Nside;

		if (outer_prev->second.pixel[i] < -1e30 || outer_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[k] < -1e30) return false;
		if (outer_next->second.pixel[i] < -1e30 || outer_next->second.pixel[j] < -1e30 || outer_next->second.pixel[k] < -1e30) return false;

		tidal[4] -= linear_interpolation(outer_prev->second.pixel[i] + outer_prev->second.pixel[j] + outer_prev->second.pixel[k], outer_next->second.pixel[i] + outer_next->second.pixel[j] + outer_next->second.pixel[k], weight);

		d = (6. * outer_prev->second.hdr.Nside - 2.) / (9. * outer_prev->second.hdr.Nside);
		d *= 12. * sqrt(1. - d*d) * dr * dr;

		tidal[2] /= d;
		tidal[4] /= d;

		for (i = 0; i < 6; i++)
			result[i] = tidal[i];

		return true;
	}

	vec2pix_ring64(inner_prev->second.hdr.Nside, v, &i);
	vec2pix_ring64(center_prev->second.hdr.Nside, v, &j);
	vec2pix_ring64(outer_prev->second.hdr.Nside, v, &k);

	if (inner_prev->second.hdr.Npix <= i || center_prev->second.hdr.Npix <= j || outer_prev->second.hdr.Npix <= k) return false;
	if (inner_next->second.hdr.Npix <= i || center_next->second.hdr.Npix <= j || outer_next->second.hdr.Npix <= k) return false;
	if (inner_prev->second.pixel[i] < -1e30 || center_prev->second.pixel[j] < -1e30 || outer_prev->second.pixel[k] < -1e30) return false;
	if (inner_next->second.pixel[i] < -1e30 || center_next->second.pixel[j] < -1e30 || outer_next->second.pixel[k] < -1e30) return false;

	tidal[5] = linear_interpolation(inner_prev->second.pixel[i], inner_next->second.pixel[i], weight);
	tidal[4] = linear_interpolation(center_prev->second.pixel[j], center_next->second.pixel[j], weight);
	tidal[0] = linear_interpolation(outer_prev->second.pixel[k], outer_next->second.pixel[k], weight) - 2. * tidal[4] + tidal[5];
	tidal[0] /= dr * dr;

	tidal[3] = tidal[4] - tidal[5];
	tidal[3] /= dr * center_prev->second.hdr.distance;

	tidal[5] = tidal[3];

	if (j < 4) // at north pole
	{
		if (debug) cout << " at north pole" << endl;
		if (center_prev->second.hdr.Npix < 12 || center_next->second.hdr.Npix < 12) return false;
		if (center_prev->second.pixel[(j+1)%4] < -1e30 || center_prev->second.pixel[(j+3)%4] < -1e30 || center_next->second.pixel[(j+1)%4] < -1e30 || center_next->second.pixel[(j+3)%4] < -1e30) return false;
		if (center_prev->second.pixel[2*j+4] < -1e30 || center_prev->second.pixel[2*j+5] < -1e30 || center_next->second.pixel[2*j+4] < -1e30 || center_next->second.pixel[2*j+5] < -1e30) return false;

		cth = 1. - 1. / (3. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside);
		sth = sqrt(1. - cth*cth);

		sph = (j > 1) ? -0.7071067811865475 : 0.7071067811865475;
		cph = (j % 2) ? -sph : sph;

		tidal[1] = linear_interpolation(center_prev->second.pixel[2*j+4] + center_prev->second.pixel[2*j+5], center_next->second.pixel[2*j+4] + center_next->second.pixel[2*j+5], weight);
		tidal[2] = linear_interpolation(center_prev->second.pixel[(j+1)%4] + center_prev->second.pixel[(j+3)%4], center_next->second.pixel[(j+1)%4] + center_next->second.pixel[(j+3)%4], weight);

		tidal[3] /= 3. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside - 0.5;
		tidal[3] += (0.25 * (tidal[1] + tidal[2]) - tidal[4]) / center_prev->second.hdr.distance / center_prev->second.hdr.distance;
		tidal[3] *= 3. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside;

		tidal[5] += (3. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside - 1.) * (tidal[4] - 0.5 * tidal[2]) / center_prev->second.hdr.distance / center_prev->second.hdr.distance;
		tidal[5] *= sqrt(36. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside - 6);
		tidal[5] += 9. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside * (tidal[2] - 2. * tidal[4]) / center_prev->second.hdr.distance / center_prev->second.hdr.distance;
		tidal[5] *= (double) center_prev->second.hdr.Nside / (6. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside - 1.);

		tidal[1] = 0.306186217847897 * (tidal[2] - tidal[1]) * (double) center_prev->second.hdr.Nside / center_prev->second.hdr.distance;

		tidal[2] = 1.95984444731456 * linear_interpolation(center_prev->second.pixel[2*j+5] - center_prev->second.pixel[2*j+4], center_next->second.pixel[2*j+5] - center_next->second.pixel[2*j+4], weight) * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside / sqrt(6. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside - 4.) / center_prev->second.hdr.distance;

		tidal[4] = tidal[2] - 1.5 * linear_interpolation(center_prev->second.pixel[(j+1)%4] - center_prev->second.pixel[(j+3)%4], center_next->second.pixel[(j+1)%4] - center_next->second.pixel[(j+3)%4], weight) * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside / sqrt(6. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside - 1.) / center_prev->second.hdr.distance;

		tidal[2] = tidal[2] - 0.5 * tidal[4];

		tidal[4] *= 3. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside / sqrt(9. + 6.82842712474619 * (3. * (double) center_prev->second.hdr.Nside * (double) center_prev->second.hdr.Nside - 2.)) / center_prev->second.hdr.distance;

		if (d > 0 && inner_prev->second.hdr.distance > 0 && inner_prev->second.hdr.Nside == center_prev->second.hdr.Nside && inner_prev->second.hdr.Npix >= 12 && inner_next->second.hdr.Npix >= 12)
		{
			if (inner_prev->second.pixel[(j+1)%4] < -1e30 || inner_prev->second.pixel[(j+3)%4] < -1e30 || inner_next->second.pixel[(j+1)%4] < -1e30 || inner_next->second.pixel[(j+3)%4] < -1e30) return false;
			if (inner_prev->second.pixel[2*j+4] < -1e30 || inner_prev->second.pixel[2*j+5] < -1e30 || inner_next->second.pixel[2*j+4] < -1e30 || inner_next->second.pixel[2*j+5] < -1e30) return false;

			tidal[1] += 0.306186217847897 * linear_interpolation(inner_prev->second.pixel[2*j+4] + inner_prev->second.pixel[2*j+5] - inner_prev->second.pixel[(j+1)%4] - inner_prev->second.pixel[(j+3)%4], inner_next->second.pixel[2*j+4] + inner_next->second.pixel[2*j+5] - inner_next->second.pixel[(j+1)%4] - inner_next->second.pixel[(j+3)%4], weight) * (double) inner_prev->second.hdr.Nside / inner_prev->second.hdr.distance;
			tidal[1] /= -dr;

			d = 0.9799222236572824 * linear_interpolation(inner_prev->second.pixel[2*j+5] - inner_prev->second.pixel[2*j+4], inner_next->second.pixel[2*j+5] - inner_next->second.pixel[2*j+4], weight) / sqrt(6. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside - 4.) + 0.75 * linear_interpolation(inner_prev->second.pixel[(j+1)%4] - inner_prev->second.pixel[(j+3)%4], inner_next->second.pixel[(j+1)%4] - inner_next->second.pixel[(j+3)%4], weight) / sqrt(6. * (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside - 1.);
			d *= (double) inner_prev->second.hdr.Nside * (double) inner_prev->second.hdr.Nside / inner_prev->second.hdr.distance;

			tidal[2] = (tidal[2] - d) / dr;
		}
		else if (outer_prev->second.hdr.Nside == center_prev->second.hdr.Nside)
		{
			if (outer_prev->second.hdr.Npix < 12 || outer_next->second.hdr.Npix < 12) return false;
			if (outer_prev->second.pixel[(j+1)%4] < -1e30 || outer_prev->second.pixel[(j+3)%4] < -1e30 || outer_next->second.pixel[(j+1)%4] < -1e30 || outer_next->second.pixel[(j+3)%4] < -1e30) return false;
			if (outer_prev->second.pixel[2*j+4] < -1e30 || outer_prev->second.pixel[2*j+5] < -1e30 || outer_next->second.pixel[2*j+4] < -1e30 || outer_next->second.pixel[2*j+5] < -1e30) return false;

			tidal[1] += 0.306186217847897 * linear_interpolation(outer_prev->second.pixel[2*j+4] + outer_prev->second.pixel[2*j+5] - outer_prev->second.pixel[(j+1)%4] - outer_prev->second.pixel[(j+3)%4], outer_next->second.pixel[2*j+4] + outer_next->second.pixel[2*j+5] - outer_next->second.pixel[(j+1)%4] - outer_next->second.pixel[(j+3)%4], weight) * (double) outer_prev->second.hdr.Nside / outer_prev->second.hdr.distance;
			tidal[1] /= dr;

			d = 0.9799222236572824 * linear_interpolation(outer_prev->second.pixel[2*j+5] - outer_prev->second.pixel[2*j+4], outer_next->second.pixel[2*j+5] - outer_next->second.pixel[2*j+4], weight) / sqrt(6. * (double) outer_prev->second.hdr.Nside * (double) outer_prev->second.hdr.Nside - 4.) + 0.75 * linear_interpolation(outer_prev->second.pixel[(j+1)%4] - outer_prev->second.pixel[(j+3)%4], outer_next->second.pixel[(j+1)%4] - outer_next->second.pixel[(j+3)%4], weight) / sqrt(6. * (double) outer_prev->second.hdr.Nside * (double) outer_prev->second.hdr.Nside - 1.);
			d *= (double) outer_prev->second.hdr.Nside * (double) outer_prev->second.hdr.Nside / outer_prev->second.hdr.distance;

			tidal[2] = (d - tidal[2]) / dr;
		}
		else
		{
			if (outer_prev->second.hdr.Nside != 2 * center_prev->second.hdr.Nside || outer_prev->second.hdr.Npix < 24 || outer_next->second.hdr.Npix < 24) return false;
			if (outer_prev->second.pixel[2*j+4] < -1e30 || outer_prev->second.pixel[2*j+5] < -1e30 || outer_next->second.pixel[2*j+4] < -1e30 || outer_next->second.pixel[2*j+5] < -1e30) return false;
			if (outer_prev->second.pixel[3*j+13] < -1e30 || outer_next->second.pixel[3*j+13] < -1e30) return false;

			tidal[1] += 1.73205080756888 * linear_interpolation(outer_prev->second.pixel[3*j+13] - outer_prev->second.pixel[j], outer_next->second.pixel[3*j+13] - outer_next->second.pixel[j], weight) * (double) outer_prev->second.hdr.Nside * (double) outer_prev->second.hdr.Nside / sqrt(20. * (double) outer_prev->second.hdr.Nside * (double) outer_prev->second.hdr.Nside - 6. - sqrt((144. * (double) outer_prev->second.hdr.Nside * (double) outer_prev->second.hdr.Nside - 240.) * (double) outer_prev->second.hdr.Nside * (double) outer_prev->second.hdr.Nside - 36.)) / outer_prev->second.hdr.distance;
			tidal[1] /= dr;

			d = 1.959844447314565 * linear_interpolation(outer_prev->second.pixel[2*j+5] - outer_prev->second.pixel[2*j+4], outer_next->second.pixel[2*j+5] - outer_next->second.pixel[2*j+4], weight) * (double) outer_prev->second.hdr.Nside * (double) outer_prev->second.hdr.Nside / sqrt(6. * (double) outer_prev->second.hdr.Nside * (double) outer_prev->second.hdr.Nside - 4.) / outer_prev->second.hdr.distance;

			tidal[2] = (d - tidal[2]) / dr;
		}
	}
	else if (j < (int64_t) center_prev->second.hdr.Nside * ((int64_t) center_prev->second.hdr.Nside - 1) * 2) // in north polar cap
	{
		//if (debug) cout << " in north polar cap: ";
		ring = (1 + (int64_t) sqrt(1.5 + 2 * (double) j)) / 2;
			
		if ((int64_t) center_prev->second.hdr.Npix < 2 * (ring+1) * (ring+2)) return false;

		i = j - 2 * ring * (ring-1);

		m = (i * (ring+1)) / ring;

		tidal[2] = (1.5 + m) - ((double) (ring+1) * (0.5 + i) / (double) ring);
		result[5] = 1. / (sin(0.5 * tidal[2] * M_PI / (ring + 1)) + sin(0.5 * (1.-tidal[2]) * M_PI / (ring + 1)));
		tidal[2] = sin(0.5 * tidal[2] * M_PI / (ring + 1)) * result[5];
		corr2 = 2. * sin(0.25 * M_PI / (ring+1));
		corr2 *= corr2 * tidal[2] * (1. - tidal[2]);

		k = m + 2 * ring * (ring+1);
		l = k + 1;

		if (center_prev->second.pixel[k] < -1e30 || center_prev->second.pixel[l] < -1e30 || center_next->second.pixel[k] < -1e30 || center_next->second.pixel[l] < -1e30) return false;

		//if (debug) cout << " in north polar cap: j = " << j << "; k+ = " << k << ", l+ = " << l << ", w+ = " << tidal[2];

		if (tidal[2] > 0.5)
		{
			m = (m > 0) ? k - 1 : k + 4*ring + 3;

			if (center_prev->second.pixel[m] < -1e30 || center_next->second.pixel[m] < -1e30) return false;

			cph = linear_interpolation((0.5-tidal[2]) * center_prev->second.pixel[m] + (2.*tidal[2]-2.) * center_prev->second.pixel[k] + (1.5-tidal[2]) * center_prev->second.pixel[l], (0.5-tidal[2]) * center_next->second.pixel[m] + (2.*tidal[2]-2.) * center_next->second.pixel[k] + (1.5-tidal[2]) * center_next->second.pixel[l], weight);

			//tidal[1] = linear_interpolation(0.5*(tidal[2]-1.)*tidal[2] * center_prev->second.pixel[m] + (2.-tidal[2])*tidal[2] * center_prev->second.pixel[k] + (1. + 0.5*(tidal[2]-3.)*tidal[2]) * center_prev->second.pixel[l], 0.5*(tidal[2]-1.)*tidal[2] * center_next->second.pixel[m] + (2.-tidal[2])*tidal[2] * center_next->second.pixel[k] + (1. + 0.5*(tidal[2]-3.)*tidal[2]) * center_next->second.pixel[l], weight);
		}
		else
		{
			m = (m < 4*ring + 2) ? l + 1 : 2 * ring * (ring+1);

			if (center_prev->second.pixel[m] < -1e30 || center_next->second.pixel[m] < -1e30) return false;

			cph = linear_interpolation((0.5-tidal[2]) * center_prev->second.pixel[m] + 2.*tidal[2] * center_prev->second.pixel[l] - (0.5+tidal[2]) * center_prev->second.pixel[k], (0.5-tidal[2]) * center_next->second.pixel[m] + 2.*tidal[2] * center_next->second.pixel[l] - (0.5+tidal[2]) * center_next->second.pixel[k], weight);

			//tidal[1] = linear_interpolation(0.5*(tidal[2]+1.)*tidal[2] * center_prev->second.pixel[k] + (1.-tidal[2]*tidal[2]) * center_prev->second.pixel[l] + 0.5*(tidal[2]-1.)*tidal[2] * center_prev->second.pixel[m], 0.5*(tidal[2]+1.)*tidal[2] * center_next->second.pixel[k] + (1.-tidal[2]*tidal[2]) * center_next->second.pixel[l] + 0.5*(tidal[2]-1.)*tidal[2] * center_next->second.pixel[m], weight);
		}

/*		if (debug)
		{
			pix2ang_ring64(center_prev->second.hdr.Nside, j, result, result+1);
			pix2ang_ring64(center_prev->second.hdr.Nside, k, result+2, result+3);
			pix2ang_ring64(center_prev->second.hdr.Nside, l, result+4, result+5);
			cout << " ring+1 (j,k,l) = " << j << " (" << result[0] << "," << result[1] << "), " << k << " (" << result[2] << "," << result[3] << "), " << l << " (" << result[4] << "," << result[5] << "), wt = " << tidal[2];
		}*/

		result[5] *= sin(0.5 * M_PI / (ring + 1));
		tidal[1] = linear_interpolation(tidal[2] * center_prev->second.pixel[k] + (1.-tidal[2]) * center_prev->second.pixel[l], tidal[2] * center_next->second.pixel[k] + (1.-tidal[2]) * center_next->second.pixel[l], weight);

		if (i % ring == 0)
			m = ((4 + i / ring) * (ring-1) - 1) % (4 * (ring-1));
		else
			m = (i * (ring-1)) / ring;

		tidal[2] = (1.5 + m) - ((double) (ring-1) * (0.5 + i) / (double) ring);
		result[3] = 1. / (sin(0.5 * tidal[2] * M_PI / (ring - 1)) + sin(0.5 * (1.-tidal[2]) * M_PI / (ring - 1)));
		tidal[2] = sin(0.5 * tidal[2] * M_PI / (ring - 1)) * result[3];
		corr1 = 2. * sin(0.25 * M_PI / (ring-1));
		corr1 *= corr1 * tidal[2] * (1. - tidal[2]);

		if (i == 0) tidal[2] -= 4 * (ring-1);

		l = ((m+1) % (4 * (ring-1))) + 2 * (ring-1) * (ring-2);
		k = m + 2 * (ring-1) * (ring-2);

		if (center_prev->second.pixel[k] < -1e30 || center_prev->second.pixel[l] < -1e30 || center_next->second.pixel[k] < -1e30 || center_next->second.pixel[l] < -1e30) return false;

		//if (debug) cout << "; k- = " << k << ", l- = " << l << ", w- = " << tidal[2] << endl;

		if (tidal[2] > 0.5)
		{
			m = (m > 0) ? k - 1 : k + 4*ring - 5;

			if (center_prev->second.pixel[m] < -1e30 || center_next->second.pixel[m] < -1e30) return false;

			sph = linear_interpolation((0.5-tidal[2]) * center_prev->second.pixel[m] + (2.*tidal[2]-2.) * center_prev->second.pixel[k] + (1.5-tidal[2]) * center_prev->second.pixel[l], (0.5-tidal[2]) * center_next->second.pixel[m] + (2.*tidal[2]-2.) * center_next->second.pixel[k] + (1.5-tidal[2]) * center_next->second.pixel[l], weight);

			//tidal[2] = linear_interpolation(0.5*(tidal[2]-1.)*tidal[2] * center_prev->second.pixel[m] + (2.-tidal[2])*tidal[2] * center_prev->second.pixel[k] + (1. + 0.5*(tidal[2]-3.)*tidal[2]) * center_prev->second.pixel[l], 0.5*(tidal[2]-1.)*tidal[2] * center_next->second.pixel[m] + (2.-tidal[2])*tidal[2] * center_next->second.pixel[k] + (1. + 0.5*(tidal[2]-3.)*tidal[2]) * center_next->second.pixel[l], weight);
		}
		else
		{
			m = (m == 4*ring - 6) ? 2 * (ring-1) * (ring-2) : l + 1;

			if (center_prev->second.pixel[m] < -1e30 || center_next->second.pixel[m] < -1e30) return false;

			sph = linear_interpolation((0.5-tidal[2]) * center_prev->second.pixel[m] + 2.*tidal[2] * center_prev->second.pixel[l] - (0.5+tidal[2]) * center_prev->second.pixel[k], (0.5-tidal[2]) * center_next->second.pixel[m] + 2.*tidal[2] * center_next->second.pixel[l] - (0.5+tidal[2]) * center_next->second.pixel[k], weight);

			//tidal[2] = linear_interpolation(0.5*(tidal[2]+1.)*tidal[2] * center_prev->second.pixel[k] + (1.-tidal[2]*tidal[2]) * center_prev->second.pixel[l] + 0.5*(tidal[2]-1.)*tidal[2] * center_prev->second.pixel[m], 0.5*(tidal[2]+1.)*tidal[2] * center_next->second.pixel[k] + (1.-tidal[2]*tidal[2]) * center_next->second.pixel[l] + 0.5*(tidal[2]-1.)*tidal[2] * center_next->second.pixel[m], weight);
		}

/*		if (debug)
		{
			pix2ang_ring64(center_prev->second.hdr.Nside, k, result+2, result+3);
			pix2ang_ring64(center_prev->second.hdr.Nside, l, result+4, result+5);
			cout << ", ring-1 (k,l) = " << k << " (" << result[2] << "," << result[3] << "), " << l << " (" << result[4] << "," << result[5] << "), wt = " << tidal[2] << endl; // << ", factor = " << 3. / (ring*ring + 1 - (ring*ring - 1) * (ring*ring - 1 + sqrt((6. * center_prev->second.hdr.Nside * center_prev->second.hdr.Nside - 1.) * (6. * center_prev->second.hdr.Nside * center_prev->second.hdr.Nside - 1.) - (12. * center_prev->second.hdr.Nside * center_prev->second.hdr.Nside + 2) *ring*ring + (double) ring * ring * ring * ring)) / 6. / center_prev->second.hdr.Nside / center_prev->second.hdr.Nside) << endl;
		}*/

		result[3] *= sin(0.5 * M_PI / (ring - 1));
		tidal[2] = linear_interpolation(tidal[2] * center_prev->second.pixel[k] + (1.-tidal[2]) * center_prev->second.pixel[l], tidal[2] * center_next->second.pixel[k] + (1.-tidal[2]) * center_next->second.pixel[l], weight);

		//if (debug) cout << ", pixel values: " << tidal[1] << ", " << tidal[4] << ", " << tidal[2] << ", distance = " << center_prev->second.hdr.distance << endl;

		cph *= center_prev->second.hdr.Nside / sqrt(6. - (ring+1)*(ring+1) / (double) center_prev->second.hdr.Nside / (double) center_prev->second.hdr.Nside) / (ring+1) / sin(0.25 * M_PI / (ring+1));
		cph -= sph * center_prev->second.hdr.Nside / sqrt(6. - (ring-1)*(ring-1) / (double) center_prev->second.hdr.Nside / (double) center_prev->second.hdr.Nside) / (ring-1) / sin(0.25 * M_PI / (ring-1));

		result[0] = 1. - ring*ring / (3. * center_prev->second.hdr.Nside*center_prev->second.hdr.Nside);
		result[1] = sqrt(1. - result[0]*result[0]);
		/*result[2] = 1. - (ring-1)*(ring-1) / (3. * center_prev->second.hdr.Nside*center_prev->second.hdr.Nside);
		result[3] = sqrt(1. - result[2]*result[2]);
		result[4] = 1. - (ring+1)*(ring+1) / (3. * center_prev->second.hdr.Nside*center_prev->second.hdr.Nside);
		result[5] = sqrt(1. - result[4]*result[4]);*/

		p1 = 1. - (ring-1)*(ring-1) / (3. * center_prev->second.hdr.Nside*center_prev->second.hdr.Nside);
		p2 = 1. - (ring+1)*(ring+1) / (3. * center_prev->second.hdr.Nside*center_prev->second.hdr.Nside);

		result[3] *= sqrt(1. - p1*p1) / p1;
		result[5] *= sqrt(1. - p2*p2) / p2;

		result[2] = 1. / sqrt(1. + result[3]*result[3]); // cos(th1)
		result[4] = 1. / sqrt(1. + result[5]*result[5]); // cos(th2)

		result[3] *= result[2]; // sin(th1)
		result[5] *= result[4]; // sin(th2)

		p1 /= result[2];
		p2 /= result[4];

		h = sqrt(p1 * p2 * (result[5]*result[2] - result[4]*result[3]) / (p1 * (result[1]*result[2] - result[0]*result[3]) + p2 * (result[5]*result[0] - result[4]*result[1])));

		//sa = p1 * (result[1]*result[2] - result[0]*result[3]);

		p1 = sqrt(h*h + p1*p1 - 2.*h*p1*(result[0]*result[2] + result[1]*result[3]));
		p2 = sqrt(h*h + p2*p2 - 2.*h*p2*(result[0]*result[4] + result[1]*result[5]));

		//sa /= p1;

		corr1 *= result[3]*result[3] / (result[1]*result[2] - result[0]*result[3]);
		corr2 *= result[5]*result[5] / (result[5]*result[0] - result[4]*result[1]);

		corr1 = (corr1 + corr2) / (result[1]*result[2] - result[0]*result[3] + result[5]*result[0] - result[4]*result[1]);

		//tidal[3] /= 0.5 * (1. + result[0]*result[2] + result[1]*result[3]);

		tidal[3] *= (2. - 2. * h) / p1 * p2;

		//tidal[3] += (tidal[1] + tidal[2] - 2. * tidal[4]) * 2. / (1. - result[2]*result[4] - result[3]*result[5]) / center_prev->second.hdr.distance / center_prev->second.hdr.distance;

		tidal[3] += 2. * ((tidal[1] - tidal[4]) / p2 + (tidal[2] - tidal[4]) / p1) / (p1 + p2) / center_prev->second.hdr.distance / center_prev->second.hdr.distance;

		//tidal[3] += ((tidal[1] + tidal[2] - 2. * tidal[4]) / (result[1]*result[2] - result[0]*result[3])) / (result[1]*result[2] - result[0]*result[3]) / center_prev->second.hdr.distance / center_prev->second.hdr.distance;

		/*tidal[3] /= 1. - ((2 * ring * (ring-1) + 1) / 6. - ring * (ring-1) * (ring * (ring-1) + sqrt((6. * center_prev->second.hdr.Nside * center_prev->second.hdr.Nside - (ring-1) * (ring-1)) * (6. * center_prev->second.hdr.Nside * center_prev->second.hdr.Nside - ring * ring))) / 18. / center_prev->second.hdr.Nside / center_prev->second.hdr.Nside) / center_prev->second.hdr.Nside / center_prev->second.hdr.Nside;

		tidal[3] += (tidal[1] + tidal[2] - 2. * tidal[4]) * 3. * center_prev->second.hdr.Nside * center_prev->second.hdr.Nside / (ring*ring + 1 - (ring*ring - 1) * (ring*ring - 1 + sqrt((6. * center_prev->second.hdr.Nside * center_prev->second.hdr.Nside - 1.) * (6. * center_prev->second.hdr.Nside * center_prev->second.hdr.Nside - 1.) - (12. * center_prev->second.hdr.Nside * center_prev->second.hdr.Nside + 2) *ring*ring + (double) ring * ring * ring * ring)) / 6. / center_prev->second.hdr.Nside / center_prev->second.hdr.Nside) / center_prev->second.hdr.distance / center_prev->second.hdr.distance;*/

		//tidal[3] += ((tidal[2]-tidal[4])/(result[5]*result[0]-result[4]*result[1]) - (tidal[4]-tidal[1])/(result[1]*result[2]-result[0]*result[3]))*((2.*result[0]*result[0] - 1.)*(result[2]*result[4] - result[3]*result[5])+2.*result[0]*result[1]*(result[2]*result[5]+result[3]*result[4]))/sqrt(0.5*(1. - result[2]*result[4] - result[3]*result[5])) / center_prev->second.hdr.distance / center_prev->second.hdr.distance;

		//tidal[3] += ((tidal[1] - tidal[4])/((ring+1)*(3.*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside - ring*ring)*sqrt(6.*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside - (ring+1)*(ring+1)) - ring*(3.*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside - (ring+1)*(ring+1))*sqrt(6.*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside - ring*ring)) - (tidal[4]-tidal[2])/(ring*(3.*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside - (ring-1)*(ring-1))*sqrt(6.*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside - ring*ring) - (ring-1)*(3.*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside - ring*ring)*sqrt(6.*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside - (ring-1)*(ring-1))))*3.*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside*(ring*sqrt((6.*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside - ring*ring)*(6.*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside*(ring*ring + 1) + (ring*ring - 1)*(1 - ring*ring + sqrt((6.*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside - 1.)*(6.*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside - 1.) - 2. * (6.*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside + 1.)*ring*ring + (double)ring*ring*ring*ring)))) + (3.*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside - ring*ring)*sqrt(18.*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside - 6.*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside*(ring*ring + 1) - (ring*ring - 1)*(1 - ring*ring + sqrt((6.*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside - 1.)*(6.*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside - 1.) - 2. * (6.*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside + 1.)*ring*ring + (double)ring*ring*ring*ring))))/sqrt(6.*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside*(ring*ring + 1) - (ring*ring - 1)*(ring*ring - 1 + sqrt((6.*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside - 1.)*(6.*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside - 1.) - 2.*(6.*center_prev->second.hdr.Nside*center_prev->second.hdr.Nside + 1.)*ring*ring + (double)ring*ring*ring*ring))) / center_prev->second.hdr.distance / center_prev->second.hdr.distance;

		cth = 6*ring*(ring-1) + 3 - ring*(ring-1) * (ring*(ring-1) + sqrt((6. * center_prev->second.hdr.Nside * center_prev->second.hdr.Nside - (ring-1)*(ring-1)) * (6. * center_prev->second.hdr.Nside * center_prev->second.hdr.Nside - ring*ring))) / center_prev->second.hdr.Nside / center_prev->second.hdr.Nside;

		tidal[5] *= ((ring - 3. * center_prev->second.hdr.Nside * center_prev->second.hdr.Nside / (double) ring) * sqrt(cth / (12. * center_prev->second.hdr.Nside * center_prev->second.hdr.Nside - 2*ring*ring)) + sqrt(9. * center_prev->second.hdr.Nside * center_prev->second.hdr.Nside - 0.5 * cth)) / 3. / center_prev->second.hdr.Nside;
		tidal[5] += (tidal[4] - tidal[2]) * (9. * center_prev->second.hdr.Nside * center_prev->second.hdr.Nside - 3*ring*ring) / (double) ring / sqrt((12. - 2*ring*ring / (double) center_prev->second.hdr.Nside / (double) center_prev->second.hdr.Nside) * cth) / center_prev->second.hdr.distance / center_prev->second.hdr.distance;
		tidal[5] *= 2. / (1. + cos(0.5 * M_PI / ring));

		tidal[1] = (tidal[1] - tidal[2]) / center_prev->second.hdr.distance;

		k = (i == 0) ? j-1+4*ring : j-1;
		l = (i == 4*ring-1) ? j+1-4*ring : j+1;

		if (center_prev->second.pixel[k] < -1e30 || center_prev->second.pixel[l] < -1e30 || center_next->second.pixel[k] < -1e30 || center_next->second.pixel[l] < -1e30) return false;

		cth = 1. - ring*ring / (3. * center_prev->second.hdr.Nside * center_prev->second.hdr.Nside);
		sth = sqrt(1. - cth*cth);

		sph = sin(0.5 * M_PI / ring);

		tidal[5] += (linear_interpolation(center_prev->second.pixel[k] + center_prev->second.pixel[l], center_next->second.pixel[k] + center_next->second.pixel[l], weight) - 2. * tidal[4]) / sph / sph / sth / sth / center_prev->second.hdr.distance / center_prev->second.hdr.distance;

		tidal[3] -= corr1 * tidal[5];

		tidal[2] = linear_interpolation(center_prev->second.pixel[l] - center_prev->second.pixel[k], center_next->second.pixel[l] - center_next->second.pixel[k], weight) * 0.5 / sph / sth / center_prev->second.hdr.distance;

		tidal[4] = 1.5 * cph / center_prev->second.hdr.distance / center_prev->second.hdr.distance;

		cph = 3. * center_prev->second.hdr.Nside * center_prev->second.hdr.Nside / sqrt(12. * center_prev->second.hdr.Nside * center_prev->second.hdr.Nside * (ring*ring + 1) - 2. * (ring*ring - 1) * (ring*ring - 1 + sqrt((6. * center_prev->second.hdr.Nside * center_prev->second.hdr.Nside - 1.) * (6. * center_prev->second.hdr.Nside * center_prev->second.hdr.Nside - 1.) - ring*ring * (12. * center_prev->second.hdr.Nside * center_prev->second.hdr.Nside - ring*ring + 2))));

		tidal[4] *= cph;

		if (d > 0 && inner_prev->second.hdr.distance > 0 && inner_prev->second.hdr.Nside == center_prev->second.hdr.Nside && inner_prev->second.hdr.Npix >= 2 * (ring+1) * (ring+2) && inner_next->second.hdr.Npix >= 2 * (ring+1) * (ring+2))
		{
			if (inner_prev->second.pixel[k] < -1e30 || inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[k] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

			tidal[2] = (tidal[2] - linear_interpolation(inner_prev->second.pixel[l] - inner_prev->second.pixel[k], inner_next->second.pixel[l] - inner_next->second.pixel[k], weight) * 0.5 / sph / sth / inner_prev->second.hdr.distance) / dr;

			k = (i * (ring+1)) / ring;

			d = (1.5 + k) - ((double) (ring+1) * (0.5 + i) / (double) ring);

			k += 2 * ring * (ring+1);
			l = k + 1;

			if (inner_prev->second.pixel[k] < -1e30 || inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[k] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

			tidal[1] -= linear_interpolation(d * inner_prev->second.pixel[k] + (1.-d) * inner_prev->second.pixel[l], d * inner_next->second.pixel[k] + (1.-d) * inner_next->second.pixel[l], weight) / inner_prev->second.hdr.distance;

			if (i % ring == 0)
				k = ((4 + i / ring) * (ring-1) - 1) % (4 * (ring-1));
			else
				k = (i * (ring-1)) / ring;

			d = (1.5 + k) - ((double) (ring-1) * (0.5 + i) / (double) ring);

			if (i == 0) k -= 4 * (ring-1);

			l = ((k+1) % (4 * (ring-1))) + 2 * (ring-1) * (ring-2);
			k += 2 * (ring-1) * (ring-2);

			if (inner_prev->second.pixel[k] < -1e30 || inner_prev->second.pixel[l] < -1e30 || inner_next->second.pixel[k] < -1e30 || inner_next->second.pixel[l] < -1e30) return false;

			tidal[1] += linear_interpolation(d * inner_prev->second.pixel[k] + (1.-d) * inner_prev->second.pixel[l], d * inner_next->second.pixel[k] + (1.-d) * inner_next->second.pixel[l], weight) / inner_prev->second.hdr.distance;

			tidal[1] *= cph / dr;
		}
		else if (outer_prev->second.hdr.Nside == center_prev->second.hdr.Nside)
		{
			if (outer_prev->second.hdr.Npix < 2 * (ring+1) * (ring+2) || outer_next->second.hdr.Npix < 2 * (ring+1) * (ring+2)) return false;
			if (outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[l] < -1e30 || outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[l] < -1e30) return false;

			tidal[2] = (linear_interpolation(outer_prev->second.pixel[l] - outer_prev->second.pixel[k], outer_next->second.pixel[l] - outer_next->second.pixel[k], weight) * 0.5 / sph / sth / outer_prev->second.hdr.distance - tidal[2]) / dr;

			k = (i * (ring+1)) / ring;

			d = (1.5 + k) - ((double) (ring+1) * (0.5 + i) / (double) ring);

			k += 2 * ring * (ring+1);
			l = k + 1;

			if (outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[l] < -1e30 || outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[l] < -1e30) return false;

			tidal[1] -= linear_interpolation(d * outer_prev->second.pixel[k] + (1.-d) * outer_prev->second.pixel[l], d * outer_next->second.pixel[k] + (1.-d) * outer_next->second.pixel[l], weight) / outer_prev->second.hdr.distance;

			if (i % ring == 0)
				k = ((4 + i / ring) * (ring-1) - 1) % (4 * (ring-1));
			else
				k = (i * (ring-1)) / ring;

			d = (1.5 + k) - ((double) (ring-1) * (0.5 + i) / (double) ring);

			if (i == 0) k -= 4 * (ring-1);

			l = ((k+1) % (4 * (ring-1))) + 2 * (ring-1) * (ring-2);
			k += 2 * (ring-1) * (ring-2);

			if (outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[l] < -1e30 || outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[l] < -1e30) return false;

			tidal[1] += linear_interpolation(d * outer_prev->second.pixel[k] + (1.-d) * outer_prev->second.pixel[l], d * outer_next->second.pixel[k] + (1.-d) * outer_next->second.pixel[l], weight) / outer_prev->second.hdr.distance;

			tidal[1] *= -cph / dr;
		}
		else
		{
			if (outer_prev->second.hdr.Nside != 2 * center_prev->second.hdr.Nside || outer_prev->second.hdr.Npix < 2 * (2*ring+1) * (2*ring+2) || outer_next->second.hdr.Npix < 2 * (2*ring+1) * (2*ring+2)) return false;
			tidal[1] *= cph;

			k = 2*i + 4*ring*(2*ring - 1);
			l = k+1;

			if (outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[l] < -1e30 || outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[l] < -1e30) return false;

			cph = linear_interpolation(outer_prev->second.pixel[l] - outer_prev->second.pixel[k], outer_next->second.pixel[l] - outer_next->second.pixel[k], weight);

			tidal[2] = (cph * 0.5 / sin(0.125 * M_PI / ring)  / sth / outer_prev->second.hdr.distance - tidal[2]) / dr;

			k = (int64_t) floor((2*ring - 1)*(i + 0.5) / ring) + 4*(2*ring - 1)*(ring-1);
			l = (int64_t) floor((2*ring + 1)*(i + 0.5) / ring) + 4*ring*(2*ring + 1);

			if (outer_prev->second.pixel[k] < -1e30 || outer_prev->second.pixel[l] < -1e30 || outer_next->second.pixel[k] < -1e30 || outer_next->second.pixel[l] < -1e30) return false;

			d = (k + 0.5) / (2*ring-1) - (l + 0.5) / (2*ring + 1);

			tidal[1] = (linear_interpolation(outer_prev->second.pixel[l] - outer_prev->second.pixel[k], outer_next->second.pixel[l] - outer_next->second.pixel[k], weight) + d * cph) * 3. * outer_prev->second.hdr.Nside * outer_prev->second.hdr.Nside / sqrt(12. * outer_prev->second.hdr.Nside * outer_prev->second.hdr.Nside * (4*ring*ring + 1) - 2. * (4*ring*ring - 1) * (4*ring*ring - 1 + sqrt((6. * outer_prev->second.hdr.Nside * outer_prev->second.hdr.Nside - 1.) * (6. * outer_prev->second.hdr.Nside * outer_prev->second.hdr.Nside - 1.) - 4*ring*ring * (12. * outer_prev->second.hdr.Nside * outer_prev->second.hdr.Nside - 4*ring*ring + 2)))) / outer_prev->second.hdr.distance - tidal[1];

			tidal[1] /= dr;
		}

		cph = cos((i + 0.5) * M_PI / (2*ring));
		sph = (i < 2*ring) ? sqrt(1. - cph*cph) : -sqrt(1. - cph*cph);
	}
	else return false;

	// rotate

	result[0] = cph*cph * (tidal[0]*sth*sth + tidal[3]*cth*cth + 2*tidal[1]*cth*sth) + sph * (tidal[5]*sph - 2*cph * (tidal[2]*sth + tidal[4]*cth));
	result[1] = cph*sph * (tidal[0]*sth*sth + tidal[3]*cth*cth + 2*tidal[1]*cth*sth - tidal[5]) + (cph*cph - sph*sph) * (tidal[2]*sth + tidal[4]*cth);
	result[2] = cph * ((tidal[0] - tidal[3])*cth*sth + tidal[1]*(cth*cth - sth*sth)) + sph * (tidal[5]*sth - tidal[2]*cth);
	result[3] = sph*sph * (tidal[0]*sth*sth + tidal[3]*cth*cth + 2*tidal[1]*cth*sth) + cph * (tidal[5]*cph + 2*sph * (tidal[2]*sth + tidal[4]*cth));
	result[4] = sph * ((tidal[0] - tidal[3])*cth*sth + tidal[1]*(cth*cth - sth*sth)) - cph * (tidal[5]*sth - tidal[2]*cth);
	result[5] = tidal[0]*cth*cth - 2*tidal[1]*cth*sth + tidal[3]*sth*sth;

	return true;
}



