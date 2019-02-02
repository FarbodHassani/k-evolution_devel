#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <set>
#include "metadata.hpp"
#include "parser.hpp"
#include "background.hpp"


// g++ prepare_lc.cpp -o prepare_lc -std=c++11 -O3 -DLIGHTCONE_THICKNESS=3 -DCOLORTERMINAL -lgsl -lgslcblas

using namespace std;

#define BATCHSIZE 4194304

double pointInShell(double * pos, lightcone_geometry & lightcone, double & outer, double & inner, double * vertex)
{
	double d;

	d = sqrt((pos[0]-vertex[0])*(pos[0]-vertex[0]) + (pos[1]-vertex[1])*(pos[1]-vertex[1]) + (pos[2]-vertex[2])*(pos[2]-vertex[2]));

	if (d < inner || d >= outer) return -1.;

	if (lightcone.opening > -1.)
	{
#ifdef ACOS_HACK
		if (acos(((pos[0]-vertex[0])*lightcone.direction[0] + (pos[1]-vertex[1])*lightcone.direction[1] + (pos[2]-vertex[2])*lightcone.direction[2]) / d) < acos(lightcone.opening)) return d;
#else
		if (((pos[0]-vertex[0])*lightcone.direction[0] + (pos[1]-vertex[1])*lightcone.direction[1] + (pos[2]-vertex[2])*lightcone.direction[2]) / d > lightcone.opening) return d;
#endif
		else return -1.;
	}
	else return d;
}

int findIntersectingLightcones(lightcone_geometry & lightcone, double outer, double inner, double * domain, double vertex[MAX_INTERSECTS][3])
{
	int range = (int) ceil(outer) + 1;
	int u, v, w, n = 0;
	double corner[8][3];
	double rdom, dist;

	corner[0][0] = domain[0];
	corner[0][1] = domain[1];
	corner[0][2] = domain[2];

	corner[1][0] = domain[3];
	corner[1][1] = domain[1];
	corner[1][2] = domain[2];

	corner[2][0] = domain[0];
	corner[2][1] = domain[4];
	corner[2][2] = domain[2];

	corner[3][0] = domain[3];
	corner[3][1] = domain[4];
	corner[3][2] = domain[2];

	corner[4][0] = domain[0];
	corner[4][1] = domain[1];
	corner[4][2] = domain[5];

	corner[5][0] = domain[3];
	corner[5][1] = domain[1];
	corner[5][2] = domain[5];

	corner[6][0] = domain[0];
	corner[6][1] = domain[4];
	corner[6][2] = domain[5];

	corner[7][0] = domain[3];
	corner[7][1] = domain[4];
	corner[7][2] = domain[5];

	for (u = -range; u <= range; u++)
	{
		for (v = -range; v <= range; v++)
		{
			for (w = -range; w <= range; w++)
			{
				if (n >= MAX_INTERSECTS)
				{
					cout << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": maximum number of lightcone intersects exceeds MAX_INTERSECTS = " << MAX_INTERSECTS << " for domain (" << domain[0] << ", " << domain[1] << ", " << domain[2] << ") - (" << domain[3] << ", " << domain[4] << ", " << domain[5] << "); some data may be missing in output!" << endl;
					return MAX_INTERSECTS;
				}
				vertex[n][0] = lightcone.vertex[0] + u;
				vertex[n][1] = lightcone.vertex[1] + v;
				vertex[n][2] = lightcone.vertex[2] + w;

				// first, check if domain lies outside outer sphere
				if (vertex[n][0] < domain[0])
				{
					if (vertex[n][1] < domain[1])
					{
						if (vertex[n][2] < domain[2])
						{
							if (sqrt((vertex[n][0]-corner[0][0])*(vertex[n][0]-corner[0][0]) + (vertex[n][1]-corner[0][1])*(vertex[n][1]-corner[0][1]) + (vertex[n][2]-corner[0][2])*(vertex[n][2]-corner[0][2])) > outer) continue;
						}
						else if (vertex[n][2] > domain[5])
						{
							if (sqrt((vertex[n][0]-corner[4][0])*(vertex[n][0]-corner[4][0]) + (vertex[n][1]-corner[4][1])*(vertex[n][1]-corner[4][1]) + (vertex[n][2]-corner[4][2])*(vertex[n][2]-corner[4][2])) > outer) continue;
						}
						else if (sqrt((vertex[n][0]-domain[0])*(vertex[n][0]-domain[0]) + (vertex[n][1]-domain[1])*(vertex[n][1]-domain[1])) > outer) continue;
					}
					else if (vertex[n][1] > domain[4])
					{
						if (vertex[n][2] < domain[2])
						{
							if (sqrt((vertex[n][0]-corner[2][0])*(vertex[n][0]-corner[2][0]) + (vertex[n][1]-corner[2][1])*(vertex[n][1]-corner[2][1]) + (vertex[n][2]-corner[2][2])*(vertex[n][2]-corner[2][2])) > outer) continue;
						}
						else if (vertex[n][2] > domain[5])
						{
							if (sqrt((vertex[n][0]-corner[6][0])*(vertex[n][0]-corner[6][0]) + (vertex[n][1]-corner[6][1])*(vertex[n][1]-corner[6][1]) + (vertex[n][2]-corner[6][2])*(vertex[n][2]-corner[6][2])) > outer) continue;
						}
						else if (sqrt((vertex[n][0]-domain[0])*(vertex[n][0]-domain[0]) + (vertex[n][1]-domain[4])*(vertex[n][1]-domain[4])) > outer) continue;
					}
					else
					{
						if (vertex[n][2] < domain[2])
						{
							if (sqrt((vertex[n][0]-domain[0])*(vertex[n][0]-domain[0]) + (vertex[n][2]-domain[2])*(vertex[n][2]-domain[2])) > outer) continue;
						}
						else if (vertex[n][2] > domain[5])
						{
							if (sqrt((vertex[n][0]-domain[0])*(vertex[n][0]-domain[0]) + (vertex[n][2]-domain[5])*(vertex[n][2]-domain[5])) > outer) continue;
						}
						else if (domain[0]-vertex[n][0] > outer) continue;
					}
				}
				else if (vertex[n][0] > domain[3])
				{
					if (vertex[n][1] < domain[1])
					{
						if (vertex[n][2] < domain[2])
						{
							if (sqrt((vertex[n][0]-corner[1][0])*(vertex[n][0]-corner[1][0]) + (vertex[n][1]-corner[1][1])*(vertex[n][1]-corner[1][1]) + (vertex[n][2]-corner[1][2])*(vertex[n][2]-corner[1][2])) > outer) continue;
						}
						else if (vertex[n][2] > domain[5])
						{
							if (sqrt((vertex[n][0]-corner[5][0])*(vertex[n][0]-corner[5][0]) + (vertex[n][1]-corner[5][1])*(vertex[n][1]-corner[5][1]) + (vertex[n][2]-corner[5][2])*(vertex[n][2]-corner[5][2])) > outer) continue;
						}
						else if (sqrt((vertex[n][0]-domain[3])*(vertex[n][0]-domain[3]) + (vertex[n][1]-domain[1])*(vertex[n][1]-domain[1])) > outer) continue;
					}
					else if (vertex[n][1] > domain[4])
					{
						if (vertex[n][2] < domain[2])
						{
							if (sqrt((vertex[n][0]-corner[3][0])*(vertex[n][0]-corner[3][0]) + (vertex[n][1]-corner[3][1])*(vertex[n][1]-corner[3][1]) + (vertex[n][2]-corner[3][2])*(vertex[n][2]-corner[3][2])) > outer) continue;
						}
						else if (vertex[n][2] > domain[5])
						{
							if (sqrt((vertex[n][0]-corner[7][0])*(vertex[n][0]-corner[7][0]) + (vertex[n][1]-corner[7][1])*(vertex[n][1]-corner[7][1]) + (vertex[n][2]-corner[7][2])*(vertex[n][2]-corner[7][2])) > outer) continue;
						}
						else if (sqrt((vertex[n][0]-domain[3])*(vertex[n][0]-domain[3]) + (vertex[n][1]-domain[4])*(vertex[n][1]-domain[4])) > outer) continue;
					}
					else
					{
						if (vertex[n][2] < domain[2])
						{
							if (sqrt((vertex[n][0]-domain[3])*(vertex[n][0]-domain[3]) + (vertex[n][2]-domain[2])*(vertex[n][2]-domain[2])) > outer) continue;
						}
						else if (vertex[n][2] > domain[5])
						{
							if (sqrt((vertex[n][0]-domain[3])*(vertex[n][0]-domain[3]) + (vertex[n][2]-domain[5])*(vertex[n][2]-domain[5])) > outer) continue;
						}
						else if (vertex[n][0]-domain[3] > outer) continue;
					}
				}
				else
				{
					if (vertex[n][1] < domain[1])
					{
						if (vertex[n][2] < domain[2])
						{
							if (sqrt((vertex[n][1]-domain[1])*(vertex[n][1]-domain[1]) + (vertex[n][2]-domain[2])*(vertex[n][2]-domain[2])) > outer) continue;
						}
						else if (vertex[n][2] > domain[5])
						{
							if (sqrt((vertex[n][1]-domain[1])*(vertex[n][1]-domain[1]) + (vertex[n][2]-domain[5])*(vertex[n][2]-domain[5])) > outer) continue;
						}
						else if (domain[1]-vertex[n][1] > outer) continue;
					}
					else if (vertex[n][1] > domain[4])
					{
						if (vertex[n][2] < domain[2])
						{
							if (sqrt((vertex[n][1]-domain[4])*(vertex[n][1]-domain[4]) + (vertex[n][2]-domain[2])*(vertex[n][2]-domain[2])) > outer) continue;
						}
						else if (vertex[n][2] > domain[5])
						{
							if (sqrt((vertex[n][1]-domain[4])*(vertex[n][1]-domain[4]) + (vertex[n][2]-domain[5])*(vertex[n][2]-domain[5])) > outer) continue;
						}
						else if (vertex[n][1]-domain[4] > outer) continue;
					}
					else if (vertex[n][2]-domain[5] > outer || domain[2]-vertex[n][2] > outer) continue;
				}
				
				if (sqrt((corner[0][0]-vertex[n][0])*(corner[0][0]-vertex[n][0]) + (corner[0][1]-vertex[n][1])*(corner[0][1]-vertex[n][1]) + (corner[0][2]-vertex[n][2])*(corner[0][2]-vertex[n][2])) < inner && sqrt((corner[1][0]-vertex[n][0])*(corner[1][0]-vertex[n][0]) + (corner[1][1]-vertex[n][1])*(corner[1][1]-vertex[n][1]) + (corner[1][2]-vertex[n][2])*(corner[1][2]-vertex[n][2])) < inner && sqrt((corner[2][0]-vertex[n][0])*(corner[2][0]-vertex[n][0]) + (corner[2][1]-vertex[n][1])*(corner[2][1]-vertex[n][1]) + (corner[2][2]-vertex[n][2])*(corner[2][2]-vertex[n][2])) < inner && sqrt((corner[3][0]-vertex[n][0])*(corner[3][0]-vertex[n][0]) + (corner[3][1]-vertex[n][1])*(corner[3][1]-vertex[n][1]) + (corner[3][2]-vertex[n][2])*(corner[3][2]-vertex[n][2])) < inner && sqrt((corner[4][0]-vertex[n][0])*(corner[4][0]-vertex[n][0]) + (corner[4][1]-vertex[n][1])*(corner[4][1]-vertex[n][1]) + (corner[4][2]-vertex[n][2])*(corner[4][2]-vertex[n][2])) < inner && sqrt((corner[5][0]-vertex[n][0])*(corner[5][0]-vertex[n][0]) + (corner[5][1]-vertex[n][1])*(corner[5][1]-vertex[n][1]) + (corner[5][2]-vertex[n][2])*(corner[5][2]-vertex[n][2])) < inner && sqrt((corner[6][0]-vertex[n][0])*(corner[6][0]-vertex[n][0]) + (corner[6][1]-vertex[n][1])*(corner[6][1]-vertex[n][1]) + (corner[6][2]-vertex[n][2])*(corner[6][2]-vertex[n][2])) < inner && sqrt((corner[7][0]-vertex[n][0])*(corner[7][0]-vertex[n][0]) + (corner[7][1]-vertex[n][1])*(corner[7][1]-vertex[n][1]) + (corner[7][2]-vertex[n][2])*(corner[7][2]-vertex[n][2])) < inner) continue; // domain lies within inner sphere

				rdom = 0.5 * sqrt((domain[3]-domain[0])*(domain[3]-domain[0]) + (domain[4]-domain[1])*(domain[4]-domain[1]) + (domain[5]-domain[2])*(domain[5]-domain[2]));
				dist = sqrt((0.5*domain[0]+0.5*domain[3]-vertex[n][0])*(0.5*domain[0]+0.5*domain[3]-vertex[n][0]) + (0.5*domain[1]+0.5*domain[4]-vertex[n][1])*(0.5*domain[1]+0.5*domain[4]-vertex[n][1]) + (0.5*domain[2]+0.5*domain[5]-vertex[n][2])*(0.5*domain[2]+0.5*domain[5]-vertex[n][2]));

				if (dist + outer <= rdom) // outer sphere lies within domain enclosing sphere
				{
					n++;
					continue;
				}

				if (acos(((0.5*domain[0]+0.5*domain[3]-vertex[n][0])*lightcone.direction[0] + (0.5*domain[1]+0.5*domain[4]-vertex[n][1])*lightcone.direction[1] + (0.5*domain[2]+0.5*domain[5]-vertex[n][2])*lightcone.direction[2]) / dist) - acos(lightcone.opening) <= acos((outer*outer + dist*dist - rdom*rdom) / (2. * outer * dist))) // enclosing sphere within opening
				{
					n++;
				}
			}
		}
	}

	return n;
}

int main(int argc, char **argv)
{
	char * settingsfile = NULL;
	char * cycleparam = NULL;
	char * lightconeparam = NULL;
	parameter * params = NULL;
	metadata sim;
	cosmology cosmo;
	icsettings ic;
	int numparam = 0;
	int usedparams = 0;

	bool use_lightcone[MAX_OUTPUTS];
	int num_shells[MAX_OUTPUTS];
	int start_cycle[MAX_OUTPUTS];
	int end_cycle[MAX_OUTPUTS];
	int min_cycle = 2147483647;
	int max_cycle = 0;
	double * tau_info[MAX_OUTPUTS];
	double * a_info[MAX_OUTPUTS];
	double * otau_info[MAX_OUTPUTS];
	double * itau_info[MAX_OUTPUTS];

	double dummy_read[LIGHTCONE_THICKNESS+5];
	char filename[1024];
	char ofilename[1024];
	FILE * infile;
	FILE * outfile;
	size_t num_read1, num_read2;

	double * vertex = NULL;
	double vertexlist[MAX_INTERSECTS][3];
	double domain[6] = {0., 0., 0., 1., 1., 1.};
	int numvertex;
	int tryvertex = 0;
	double inner, outer, dtau;
	double maxvel = 0.1;
	double tau_obs = -1.;
	double z_obs;
	double fourpiG;

	uint64_t numpart_tot = 0;
	uint64_t numpart_write = 0;
	uint64_t numpart_reject = 0;
	long est_reject = 0;
	uint32_t blocksize = 0;
	int numfiles = 1;
	int numread = 0;
	int numwrite = 0;
	gadget2_header hdr;
	gadget2_header outhdr;

	vector<float> posbuffer;
	vector<float> velbuffer;
	long backtrack;
	long fastforward;
	float * posbatch = NULL;
	float * velbatch = NULL;
	uint32_t batch;
#if GADGET_ID_BYTES == 8
	vector<uint64_t> IDbuffer;
	set<uint64_t> IDbacklog;
	set<uint64_t> IDprelog;
	set<uint64_t> IDlookup;
	set<uint64_t> IDomit;
	uint64_t * IDbatch = NULL;
#else
	vector<uint32_t> IDbuffer;
	set<uint32_t> IDbacklog;
	set<uint32_t> IDprelog;
	set<uint32_t> IDlookup;
	set<uint32_t> IDomit;
	uint32_t * IDbatch = NULL;
#endif

	double pos[3];
	double dist;
	double rescale_vel;
	double maxpos = 1.;
	double offset = 0.;
	int restart = 0;

	for (int i = 1 ; i < argc ; i++ ){
		if ( argv[i][0] != '-' )
			continue;
		switch(argv[i][1]) {
			case 'c':
				cycleparam = argv[++i]; // cycle range selector
				break;
			case 'l':
				lightconeparam = argv[++i]; // light cone selector
				break;
			case 'n':
				numfiles = atoi(argv[++i]); // number of output files
				break;
			case 'R':
				restart = 1;
				break;
			case 'r':
				est_reject = atol(argv[++i]); // estimated number of rejected particles (to balance file sizes)
				break;
			case 'v':
				maxvel = 2. * atof(argv[++i]); // maximum velocity (to establish overlap size)
				break;
			case 'o':
				offset = atof(argv[++i]); // position offset (to make resulting contiguous data region fit into the final cube)
				break;
			case 's':
				settingsfile = argv[++i]; //settings file name
		}
	}

	if (settingsfile == NULL)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no settings file specified!" << endl;
		return -1;
	}
	else if (numfiles < 1 || !isfinite(numfiles))
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": number of output files not recognized!" << endl;
		return -1;
	}
	
	cout << COLORTEXT_WHITE << " LCARS tools: prepare light cone" << endl << COLORTEXT_RESET << " compiled with " << COLORTEXT_RED << "LIGHTCONE_THICKNESS=" << LIGHTCONE_THICKNESS << endl << COLORTEXT_RESET << endl << " opening settings file of simulation: " << settingsfile << endl << " parser output:" << endl << endl;

	numparam = loadParameterFile(settingsfile, params);
	usedparams = parseMetadata(params, numparam, sim, cosmo, ic);

	cout << endl << " file contains " << numparam << " parameters, " << usedparams << " of which could be parsed." << endl << endl;

	fourpiG = 1.5 * sim.boxsize * sim.boxsize / C_SPEED_OF_LIGHT / C_SPEED_OF_LIGHT;
	offset /= sim.boxsize;

	cout << " number of lightcones: " << sim.num_lightcone << endl << endl;

	for (int i = 0; i < sim.num_lightcone; i++)
	{
		cout << " lightcone " << i << " parameters:" << endl << "  vertex = " << sim.lightcone[i].vertex[0] << ", " << sim.lightcone[i].vertex[1] << ", " << sim.lightcone[i].vertex[2] << endl;
		cout << "  redshift of observation = " << sim.lightcone[i].z << endl;
		cout << "  direction = " << sim.lightcone[i].direction[0] << ", " << sim.lightcone[i].direction[1] << ", " << sim.lightcone[i].direction[2] << endl;
		cout << "  opening half-angle = " << ((sim.lightcone[i].opening > -1.) ? acos(sim.lightcone[i].opening) * 180. / M_PI : 180.) /* sim.lightcone[i].opening * 180. / M_PI */ << " degrees" << endl;
		cout << "  distance interval = " << sim.lightcone[i].distance[0] << ", " << sim.lightcone[i].distance[1] << endl << endl;
	}

	if (lightconeparam == NULL)
	{
		cout << " no light cones selected (parameter -l), using all light cones with particles: ";
		for (int i = 0; i < sim.num_lightcone; i++)
		{
			if (sim.out_lightcone[i] & MASK_GADGET)
			{
				cout << i << " ";
				use_lightcone[i] = true;
			}
			else use_lightcone[i] = false;
		}
		cout << endl << endl;
	}
	else
	{
		for (int i = 0; i < sim.num_lightcone; i++) use_lightcone[i] = false;
		
		char * token = NULL;

		sprintf(filename, "%s", lightconeparam);

		token = strtok(filename, ",");

		while (token != NULL)
		{
			int i = atoi(token);

			if (i < 0 || i >= sim.num_lightcone)
			{
				cout << " light cone " << i << " (parameter -l) does not exist!" << endl;
			}
			else if (!(sim.out_lightcone[i] & MASK_GADGET))
			{
				cout << " light cone " << i << " contains no particles and will be omitted!" << endl;
			}
			else
			{
				use_lightcone[i] = true;
			}

			token = strtok(NULL, ",");
		}

		cout << " using light cone(s) ";
		for (int i = 0; i < sim.num_lightcone; i++)
		{
			if (use_lightcone[i]) cout << i << " ";
		}
		cout << endl << endl;
	}

	cout << " reading light cone information files..." << endl;

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
				sprintf(filename, "%s%s%d_info.bin", sim.output_path, sim.basename_lightcone, i);
			else
				sprintf(filename, "%s%s_info.bin", sim.output_path, sim.basename_lightcone);

			infile = fopen(filename, "rb");

			if (infile == NULL)
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to open light cone information file " << filename << endl;
				return -1;
			}

			num_read1 = fread((void *) (start_cycle+i), sizeof(int), 1, infile);
			num_read2 = fread((void *) dummy_read, sizeof(double), LIGHTCONE_THICKNESS+5, infile);

			if (num_read1 != 1 || num_read2 != LIGHTCONE_THICKNESS+5)
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": read error in light cone information file " << filename << endl;
				fclose(infile);
				return -1;
			}

			for (num_shells[i] = 1; true; num_shells[i]++)
			{
				num_read1 = fread((void *) (end_cycle+i), sizeof(int), 1, infile);
				num_read2 = fread((void *) dummy_read, sizeof(double), LIGHTCONE_THICKNESS+5, infile);

				if (feof(infile) || ferror(infile)) break;
				else if (num_read1 != 1 || num_read2 != LIGHTCONE_THICKNESS+5)
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": read error in light cone information file " << filename << " (first pass)" << endl;
					fclose(infile);
					return -1;
				}
			}

			rewind(infile);

			tau_info[i] = (double *) malloc(num_shells[i] * sizeof(double));
			a_info[i] = (double *) malloc(num_shells[i] * sizeof(double));
			otau_info[i] = (double *) malloc(num_shells[i] * sizeof(double));
			itau_info[i] = (double *) malloc(num_shells[i] * sizeof(double));

			for (int j = 0; j < num_shells[i]; j++)
			{
				num_read1 = fread((void *) (end_cycle+i), sizeof(int), 1, infile);
				num_read2 = fread((void *) dummy_read, sizeof(double), LIGHTCONE_THICKNESS+5, infile);

				if (num_read1 != 1 || num_read2 != LIGHTCONE_THICKNESS+5)
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": read error in light cone information file " << filename << " (second pass)" << endl;
					fclose(infile);
					return -1;
				}

				tau_info[i][j] = dummy_read[0];
				a_info[i][j] = dummy_read[1];
				itau_info[i][j] = dummy_read[2];
				otau_info[i][j] = dummy_read[3];
			}

			fclose(infile);

			if (start_cycle[i] < min_cycle) min_cycle = start_cycle[i];
			if (end_cycle[i] > max_cycle) max_cycle = end_cycle[i];

			cout << " light cone " << i << " contains " << num_shells[i] << " shells (cycle " << start_cycle[i] << " to " << end_cycle[i] << ")" << endl;
		}
	}

	cout << endl;

	if (cycleparam == NULL)
	{
		cout << " no range of cycles selected (parameter -c), using entire range: " << min_cycle << "-" << max_cycle << endl;
	}
	else
	{
		sprintf(filename, "%s", cycleparam);
		char * dash = strchr(filename, '-');
		int tmp;

		if (dash == NULL)
		{
			cout << " range of cycle (parameter -c) could not be interpreted!" << endl;

			for (int i = 0; i < sim.num_lightcone; i++)
			{
				if (use_lightcone[i])
				{
					free(tau_info[i]);
					free(a_info[i]);
					free(otau_info[i]);
					free(itau_info[i]);
				}
			}

			return -1;
		}

		tmp = atoi(dash+1);

		if (tmp > max_cycle)
		{
			cout << " range of cycles (parameter -c) extends beyond available range; last cycle in range is " << max_cycle << endl;
		}
		else if (tmp < min_cycle)
		{
			cout << " range of cycles (parameter -c) outside available range!" << endl;

			for (int i = 0; i < sim.num_lightcone; i++)
			{
				if (use_lightcone[i])
				{
					free(tau_info[i]);
					free(a_info[i]);
					free(otau_info[i]);
					free(itau_info[i]);
				}
			}

			return -1;
		}
		else max_cycle = tmp;

		dash[0] = '\0';
		tmp = atoi(filename);

		if (tmp < min_cycle)
		{
			cout << " range of cycles (parameter -c) extends beyond available range; first cycle in range is " << min_cycle << endl;
		}
		else if (tmp > max_cycle)
		{
			cout << " range of cycles (parameter -c) outside available range!" << endl;

			for (int i = 0; i < sim.num_lightcone; i++)
			{
				if (use_lightcone[i])
				{
					free(tau_info[i]);
					free(a_info[i]);
					free(itau_info[i]);
					free(otau_info[i]);
				}
			}

			return -1;
		}
		else min_cycle = tmp;

		cout << " range of cycles set to: " << min_cycle << "-" << max_cycle << endl;
	}

	cout << endl;

	cout << " conformal age at observation = " << tau_obs << endl << endl;

	cout << " reading particle headers..." << endl << endl;

	for (int cycle = min_cycle; cycle <= max_cycle; cycle++)
	{
		for (int i = 0; i < sim.num_lightcone; i++)
		{
			if (!use_lightcone[i]) continue;
			else if (cycle < start_cycle[i] || cycle > end_cycle[i]) continue;

			if (sim.num_lightcone > 1)
				sprintf(filename, "%s%s%d_%04d_cdm", sim.output_path, sim.basename_lightcone, i, cycle);
			else
				sprintf(filename, "%s%s_%04d_cdm", sim.output_path, sim.basename_lightcone, cycle);

			infile = fopen(filename, "rb");

			if (infile != NULL)
			{
				fread(&blocksize, sizeof(uint32_t), 1, infile);

				if (blocksize != sizeof(hdr))
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unknown file format " << filename << "!" << endl;
					fclose(infile);
					continue;
				}

				if(fread(&hdr, sizeof(hdr), 1, infile) != 1)
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read header block from " << filename << "!" << endl;
					fclose(infile);
					continue;
				}

				numpart_tot += hdr.npart[1];
				numread++;
				fclose(infile);
			}
		}
	}

	cout << " " << numread << " particle headers read successfully. Total number of particles = " << numpart_tot << endl << endl;

	for (int i = 0; i < 6; i++)
	{
		outhdr.npart[i] = 0;
		outhdr.mass[i] = 0.;
		outhdr.npartTotal[i] = 0;
		outhdr.npartTotalHW[i] = 0;
	}
	for (int i = 0; i < 256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8 - 2 * 4 - 6 * 4; i++)
		outhdr.fill[i] = 0;

	outhdr.Omega0 = cosmo.Omega_m;
	outhdr.OmegaLambda = 1. - cosmo.Omega_m;
	outhdr.HubbleParam = cosmo.h;
	outhdr.BoxSize = sim.boxsize / GADGET_LENGTH_CONVERSION;
	outhdr.flag_sfr = 0;
	outhdr.flag_cooling = 0;
	outhdr.flag_feedback = 0;
	outhdr.flag_age = 0;
	outhdr.flag_metals = 0;
	outhdr.time = 1. / (z_obs + 1.);
	outhdr.redshift = z_obs;

	outhdr.npart[1] = (uint32_t) ((((long) numpart_tot - est_reject) / numfiles) % (1l << 32));
	outhdr.npartTotal[1] = (uint32_t) (((long) numpart_tot - est_reject) % (1l << 32));
	outhdr.npartTotalHW[1] = (uint32_t) (((long) numpart_tot - est_reject) >> 32);
	outhdr.mass[1] = hdr.mass[1];
	outhdr.num_files = numfiles;

	posbuffer.reserve(3*outhdr.npart[1]);
	velbuffer.reserve(3*outhdr.npart[1]);
	IDbuffer.reserve(outhdr.npart[1]);
	
	posbatch = (float *) malloc(3 * BATCHSIZE * sizeof(float));
	velbatch = (float *) malloc(3 * BATCHSIZE * sizeof(float));
#if GADGET_ID_BYTES == 8
	IDbatch = (uint64_t *) malloc(BATCHSIZE * sizeof(uint64_t));
#else
	IDbatch = (uint32_t *) malloc(BATCHSIZE * sizeof(uint32_t));
#endif

	if (posbatch == NULL || velbatch == NULL || IDbatch == NULL)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to allocate memory for data read!" << endl;
		return -1;
	}
	
	if (restart)
	{
		sprintf(filename, "%s%s_lcars_restart.bin", sim.output_path, sim.basename_lightcone);
		infile = fopen(filename, "rb");
		
		if (infile == NULL)
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to open restart file " << filename << "!" << endl << endl;
		}
		else
		{
			fread((void *) &min_cycle, sizeof(int), 1, infile);
			fread((void *) &numwrite, sizeof(int), 1, infile);
			fread((void *) &numpart_write, sizeof(uint64_t), 1, infile);
			fread((void *) &numpart_reject, sizeof(uint64_t), 1, infile);
			fread((void *) &maxpos, sizeof(double), 1, infile);
			uint32_t backlogsize = 0;
			uint32_t recordsize = 0;
			uint32_t omitsize = 0;
			fread((void *) &backlogsize, sizeof(uint32_t), 1, infile);
			fread((void *) &recordsize, sizeof(uint32_t), 1, infile);
			fread((void *) &omitsize, sizeof(uint32_t), 1, infile);
			for (uint32_t p = 0; p < recordsize; p+=batch)
			{
				batch = (recordsize - p >= BATCHSIZE) ? BATCHSIZE : (recordsize - p);
				fread((void *) posbatch, sizeof(float), 3*batch, infile);
				for (int it = 0; it < 3 * (int) batch; it++)
					posbuffer.push_back(posbatch[it]);
			}
			for (uint32_t p = 0; p < recordsize; p+=batch)
			{
				batch = (recordsize - p >= BATCHSIZE) ? BATCHSIZE : (recordsize - p);
				fread((void *) velbatch, sizeof(float), 3*batch, infile);
				for (int it = 0; it < 3 * (int) batch; it++)
					velbuffer.push_back(velbatch[it]);
			}
			for (uint32_t p = 0; p < recordsize; p+=batch)
			{
				batch = (recordsize - p >= BATCHSIZE) ? BATCHSIZE : (recordsize - p);
				fread((void *) IDbatch, GADGET_ID_BYTES, batch, infile);
				for (int it = 0; it < (int) batch; it++)
					IDbuffer.push_back(IDbatch[it]);
			}
			for (uint32_t p = 0; p < backlogsize; p+=batch)
			{
				batch = (backlogsize - p >= BATCHSIZE) ? BATCHSIZE : (backlogsize - p);
				fread((void *) IDbatch, GADGET_ID_BYTES, batch, infile);
				for (int it = 0; it < (int) batch; it++)
					IDbacklog.insert(IDbatch[it]);
			}
			for (uint32_t p = 0; p < omitsize; p+=batch)
			{
				batch = (omitsize - p >= BATCHSIZE) ? BATCHSIZE : (omitsize - p);
				fread((void *) IDbatch, GADGET_ID_BYTES, batch, infile);
				for (int it = 0; it < (int) batch; it++)
					IDomit.insert(IDbatch[it]);
			}
			fclose(infile);
			min_cycle++;
			
			cout << " resuming from cycle " << min_cycle << ". next output is " << numwrite << ", backlog retrieved = " << IDbacklog.size() << ", buffer contains " << recordsize << " particles." << endl << endl;
		}
	}

	cout << " building up particle light cone..." << endl << endl;

	for (int cycle = min_cycle; cycle <= max_cycle; cycle++)
	{
		cout << " cycle " << cycle << " ..." << endl;
		
		for (int i = 0; i < sim.num_lightcone; i++)
		{
			if (!use_lightcone[i]) continue;
			else if (cycle < start_cycle[i] || cycle > end_cycle[i]) continue;

			rescale_vel = GADGET_VELOCITY_CONVERSION * sqrt(a_info[i][cycle-start_cycle[i]]);

			dtau = otau_info[i][cycle-start_cycle[i]] - itau_info[i][cycle-start_cycle[i]];
			outer = max<double>(otau_info[i][cycle-start_cycle[i]] + maxvel * dtau, 0.);
			inner = max<double>(itau_info[i][cycle-start_cycle[i]], 0.);
			numvertex = findIntersectingLightcones(sim.lightcone[i], outer, inner, domain, vertexlist);

			if (numvertex < 1) continue;

//			cout << " on lightcone " << i << ", outer = " << outer << ", inner = " << inner << ", number of possible vertex locations = " << numvertex << ", first vertex = " << vertexlist[0][0] << ", " << vertexlist[0][1] << ", " << vertexlist[0][2] << endl;

			if (sim.num_lightcone > 1)
				sprintf(filename, "%s%s%d_%04d_cdm", sim.output_path, sim.basename_lightcone, i, cycle);
			else
				sprintf(filename, "%s%s_%04d_cdm", sim.output_path, sim.basename_lightcone, cycle);

			infile = fopen(filename, "rb");

			if (infile != NULL)
			{
				fread(&blocksize, sizeof(uint32_t), 1, infile);

				if (blocksize != sizeof(hdr))
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unknown file format " << filename << "!" << endl;
					fclose(infile);
					continue;
				}

				if(fread(&hdr, sizeof(hdr), 1, infile) != 1)
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read header block from " << filename << "!" << endl;
					fclose(infile);
					continue;
				}
				
				fread(&blocksize, sizeof(uint32_t), 1, infile);
				fread(&blocksize, sizeof(uint32_t), 1, infile);

				long blockoffset = 3l * sizeof(float) * (long) hdr.npart[1] + 2l * sizeof(uint32_t);

				backtrack = ftell(infile);

				if (fseek(infile, 2 * blockoffset, SEEK_CUR))
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to fast forward to ID block in " << filename << "!" << endl;
					fclose(infile);
					continue;
				}

				for (uint32_t p = 0; p < hdr.npart[1]; p+=batch)
				{
					int v;

					batch = (hdr.npart[1] - p >= BATCHSIZE) ? BATCHSIZE : (hdr.npart[1] - p);

					if (
#if GADGET_ID_BYTES == 8
					fread(IDbatch, sizeof(uint64_t), batch, infile)
#else
					fread(IDbatch, sizeof(uint32_t), batch, infile)
#endif
					!= batch)
					{
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read ID batch from " << filename << "!" << endl;
						fclose(infile);
						return -1;
					}

					fastforward = ftell(infile);

					if (fseek(infile, backtrack, SEEK_SET))
					{
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to rewind to positions block in " << filename << "!" << endl;
						fclose(infile);
						return -1;
					}

					if (fread(posbatch, sizeof(float), 3l*batch, infile) != 3l*batch)
					{
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read position data from " << filename << "!" << endl;
						fclose(infile);
						return -1;
					}

					backtrack = ftell(infile);

					if (fseek(infile, blockoffset - 3l * batch * sizeof(float), SEEK_CUR))
					{
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to fast forward to velocities block in " << filename << "!" << endl;
						fclose(infile);
						return -1;
					}

					if (fread(velbatch, sizeof(float), 3l*batch, infile) != 3l*batch)
					{
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read velocity data from " << filename << "!" << endl;
						fclose(infile);
						return -1;
					}

					if (fseek(infile, fastforward, SEEK_SET))
					{
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to fast forward to ID block in " << filename << "!" << endl;
						fclose(infile);
						return -1;
					}

					for (int it = 0; it < batch; it++)
					{
#ifdef CHECK_LATERAL_OVERLAP
						if (IDbacklog.find(IDbatch[it]) == IDbacklog.end() && IDlookup.find(IDbatch[it]) == IDlookup.end())
#else
						if (IDbacklog.find(IDbatch[it]) == IDbacklog.end())
#endif
						{
							pos[0] = posbatch[3*it] / outhdr.BoxSize;
							pos[1] = posbatch[3*it+1] / outhdr.BoxSize;
							pos[2] = posbatch[3*it+2] / outhdr.BoxSize;

							for (v = 0; v < numvertex; v++)
								if ((dist = pointInShell(pos, sim.lightcone[i], outer, inner, vertexlist[(tryvertex+v)%numvertex])) >= 0.) break;

							if (v < numvertex)
							{
								tryvertex = (tryvertex+v)%numvertex;

								if (dist - inner < maxvel * dtau) IDprelog.insert(IDbatch[it]);

#ifdef CHECK_LATERAL_OVERLAP
								if (i < sim.num_lightcone-1 && acos(sim.lightcone[i].opening) - acos(((pos[0] - vertexlist[tryvertex][0])*sim.lightcone[i].direction[0] + (pos[1] - vertexlist[tryvertex][1])*sim.lightcone[i].direction[1] + (pos[2] - vertexlist[tryvertex][2])*sim.lightcone[i].direction[2]) / dist) < CHECK_LATERAL_OVERLAP) IDlookup.insert(IDbatch[it]);
#endif

								if (IDomit.find(IDbatch[it]) != IDomit.end())
								{
									//cout << " particle with ID " << IDbatch[it] << " has been recovered." << endl;
									numpart_reject += IDomit.erase(IDbatch[it]);
								}

								for (int j = 0; j < 3; j++)
								{
									pos[j] += rescale_vel * velbatch[3*it+j] * (tau_obs - tau_info[i][cycle-start_cycle[i]] - dist) - vertexlist[tryvertex][j] + vertex[j] + offset;
									if (pos[j] > maxpos) maxpos = pos[j];
									posbuffer.push_back(pos[j] * outhdr.BoxSize);
									velbuffer.push_back(velbatch[3*it+j] * rescale_vel / GADGET_VELOCITY_CONVERSION);
								}

								IDbuffer.push_back(IDbatch[it]);

								if (IDbuffer.size() == outhdr.npart[1] && numwrite < numfiles-1)
								{
									if (numfiles > 1)
										sprintf(ofilename, "%s%s_cdm.%d", sim.output_path, sim.basename_lightcone, numwrite);
									else
										sprintf(ofilename, "%s%s_cdm", sim.output_path, sim.basename_lightcone);

									outfile = fopen(ofilename, "wb");
	
									if (outfile == NULL)
									{
										fclose(infile);
										cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to open file " << ofilename << " for output!" << endl;
										return -1;
									}
	
									blocksize = sizeof(outhdr);

									cout << COLORTEXT_CYAN << " writing" << COLORTEXT_RESET << " output file " << ofilename << " ..." << endl;
									cout << " current prelog size = " << IDprelog.size()
#ifdef CHECK_LATERAL_OVERLAP
									<< ", current lookup size = " << IDlookup.size()
#endif
									<< endl;

									fwrite(&blocksize, sizeof(uint32_t), 1, outfile);
									fwrite(&outhdr, sizeof(outhdr), 1, outfile);
									fwrite(&blocksize, sizeof(uint32_t), 1, outfile);

									blocksize = posbuffer.size() * sizeof(float);

									fwrite(&blocksize, sizeof(uint32_t), 1, outfile);

									if (fwrite(posbuffer.data(), sizeof(float), posbuffer.size(), outfile) != posbuffer.size())
									{
										fclose(infile);
										fclose(outfile);
										cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to write position block!" << endl;
										return -1;
									}

									fwrite(&blocksize, sizeof(uint32_t), 1, outfile);
									fwrite(&blocksize, sizeof(uint32_t), 1, outfile);

									if (fwrite(velbuffer.data(), sizeof(float), velbuffer.size(), outfile) != velbuffer.size())
									{
										fclose(infile);
										fclose(outfile);
										cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to write velocity block!" << endl;
										return -1;
									}

									fwrite(&blocksize, sizeof(uint32_t), 1, outfile);

									blocksize = IDbuffer.size() * GADGET_ID_BYTES;

									fwrite(&blocksize, sizeof(uint32_t), 1, outfile);

									if (fwrite(IDbuffer.data(), GADGET_ID_BYTES, IDbuffer.size(), outfile) != IDbuffer.size())
									{
										fclose(infile);
										fclose(outfile);
										cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to write ID block!" << endl;
										return -1;
									}

									fwrite(&blocksize, sizeof(uint32_t), 1, outfile);

									fclose(outfile);

									cout << " output file written, contains " << outhdr.npart[1] << " particles." << endl;
	
									numpart_write += outhdr.npart[1];
									numwrite++;
									IDbuffer.clear();
									posbuffer.clear();
									velbuffer.clear();
								}
							}
							else
							{
								//cout << " warning: particle at (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ") with ID " << IDbatch[it] << " does not seem to belong to lightcone " << i << " and will be omitted. I've tried " << v << " alternative vertex location(s)." << endl;
								IDomit.insert(IDbatch[it]);
							}
						}
						else numpart_reject++;
					}
				}
				
				fclose(infile);
			}
		}

		IDbacklog = IDprelog;
	
		cout << " backlog size = " << IDbacklog.size()
#ifdef CHECK_LATERAL_OVERLAP
		 << ", lookup size = " << IDlookup.size()
#endif
		 << endl;

		IDprelog.clear();
		IDlookup.clear();
		
		if (cycle % 10 == 0)
		{
			cout << " writing restart point ..." << endl;
			
			sprintf(ofilename, "%s%s_lcars_restart.bin", sim.output_path, sim.basename_lightcone);
			outfile = fopen(ofilename, "wb");
	
			if (outfile == NULL)
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to open file " << ofilename << " for output!" << endl;
				return -1;
			}
			
			fwrite((void *) &cycle, sizeof(int), 1, outfile);
			fwrite((void *) &numwrite, sizeof(int), 1, outfile);
			fwrite((void *) &numpart_write, sizeof(uint64_t), 1, outfile);
			fwrite((void *) &numpart_reject, sizeof(uint64_t), 1, outfile);
			fwrite((void *) &maxpos, sizeof(double), 1, outfile);
			uint32_t backlogsize = IDbacklog.size();
			uint32_t recordsize = IDbuffer.size();
			uint32_t omitsize = IDomit.size();
			fwrite((void *) &backlogsize, sizeof(uint32_t), 1, outfile);
			fwrite((void *) &recordsize, sizeof(uint32_t), 1, outfile);
			fwrite((void *) &omitsize, sizeof(uint32_t), 1, outfile);
			if (fwrite(posbuffer.data(), sizeof(float), posbuffer.size(), outfile) != posbuffer.size())
			{
				fclose(outfile);
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to write position block!" << endl;
				return -1;
			}
			if (fwrite(velbuffer.data(), sizeof(float), velbuffer.size(), outfile) != velbuffer.size())
			{
				fclose(outfile);
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to write velocity block!" << endl;
				return -1;
			}
			if (fwrite(IDbuffer.data(), GADGET_ID_BYTES, IDbuffer.size(), outfile) != IDbuffer.size())
			{
				fclose(outfile);
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to write ID block!" << endl;
				return -1;
			}
#if GADGET_ID_BYTES == 8
			set<uint64_t>::iterator it = IDbacklog.begin();
#else
			set<uint32_t>::iterator it = IDbacklog.begin();
#endif
			for (uint32_t p = 0; p < backlogsize; p += batch)
			{
				batch = (backlogsize - p >= BATCHSIZE) ? BATCHSIZE : (backlogsize - p);
				for (int it2 = 0; it2 < (int) batch; it2++, it++)
					IDbatch[it2] = (*it);
				fwrite((void *) IDbatch, GADGET_ID_BYTES, batch, outfile);
			}

			it = IDomit.begin();
			for (uint32_t p = 0; p < omitsize; p += batch)
			{
				batch = (omitsize - p >= BATCHSIZE) ? BATCHSIZE : (omitsize - p);
				for (int it2 = 0; it2 < (int) batch; it2++, it++)
					IDbatch[it2] = (*it);
				fwrite((void *) IDbatch, GADGET_ID_BYTES, batch, outfile);
			}
			
			cout << " buffer with " << recordsize << " particles saved." << endl << endl;
			
			fclose(outfile);
		}
	}

	cout << endl << COLORTEXT_GREEN << " particle light cone complete." << endl << COLORTEXT_RESET << endl;

	cout << " number of duplicate particles rejected: " << numpart_reject << " (" << 100. * (double) numpart_reject / (double) numpart_tot << "%)";
	if (est_reject != 0)
		cout << " -- user estimated number was " << est_reject;
	cout << endl << " number of particles omitted: " << IDomit.size() << endl << endl;

	for (int i = 0; i < sim.num_lightcone; i++)
	{
		if (use_lightcone[i])
		{
			free(tau_info[i]);
			free(a_info[i]);
			free(itau_info[i]);
			free(otau_info[i]);
		}
	}

	free(IDbatch);
	free(posbatch);
	free(velbatch);

	cout << " correcting header information ..." << endl;

	numpart_write += IDbuffer.size();
	outhdr.npartTotal[1] = (uint32_t) (numpart_write % (1l << 32));
	outhdr.npartTotalHW[1] = (uint32_t) (numpart_write >> 32);
	outhdr.BoxSize *= maxpos;
	blocksize = sizeof(outhdr);

	for (int i = 0; i < numwrite; i++)
	{
		sprintf(ofilename, "%s%s_cdm.%d", sim.output_path, sim.basename_lightcone, i);

		outfile = fopen(ofilename, "r+b");

		if (outfile == NULL)
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to open file " << ofilename << " for output!" << endl;
			return -1;
		}

		fwrite(&blocksize, sizeof(uint32_t), 1, outfile);
		fwrite(&outhdr, sizeof(outhdr), 1, outfile);

		fclose(outfile);
	}

	outhdr.npart[1] = IDbuffer.size();

	if (numfiles > 1)
		sprintf(ofilename, "%s%s_cdm.%d", sim.output_path, sim.basename_lightcone, numwrite);
	else
		sprintf(ofilename, "%s%s_cdm", sim.output_path, sim.basename_lightcone);

	outfile = fopen(ofilename, "wb");

	if (outfile == NULL)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to open file " << ofilename << " for output!" << endl;
		return -1;
	}
	
	blocksize = sizeof(outhdr);

	cout << COLORTEXT_CYAN << " writing" << COLORTEXT_RESET << " output file " << ofilename << " ..." << endl;

	fwrite(&blocksize, sizeof(uint32_t), 1, outfile);
	fwrite(&outhdr, sizeof(outhdr), 1, outfile);
	fwrite(&blocksize, sizeof(uint32_t), 1, outfile);

	blocksize = posbuffer.size() * sizeof(float);

	fwrite(&blocksize, sizeof(uint32_t), 1, outfile);

	if (fwrite(posbuffer.data(), sizeof(float), posbuffer.size(), outfile) != posbuffer.size())
	{
		fclose(outfile);
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to write position block!" << endl;
		return -1;
	}

	fwrite(&blocksize, sizeof(uint32_t), 1, outfile);
	fwrite(&blocksize, sizeof(uint32_t), 1, outfile);

	if (fwrite(velbuffer.data(), sizeof(float), velbuffer.size(), outfile) != velbuffer.size())
	{
		fclose(outfile);
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to write velocity block!" << endl;
		return -1;
	}

	fwrite(&blocksize, sizeof(uint32_t), 1, outfile);

	blocksize = IDbuffer.size() * GADGET_ID_BYTES;

	fwrite(&blocksize, sizeof(uint32_t), 1, outfile);

	if (fwrite(IDbuffer.data(), GADGET_ID_BYTES, IDbuffer.size(), outfile) != IDbuffer.size())
	{
		fclose(outfile);
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to write ID block!" << endl;
		return -1;
	}

	fwrite(&blocksize, sizeof(uint32_t), 1, outfile);

	fclose(outfile);

	cout << " output file written, contains " << outhdr.npart[1] << " particles." << endl;

	cout << endl << COLORTEXT_GREEN << " normal completion." << COLORTEXT_RESET << endl;

	return 0;
}

