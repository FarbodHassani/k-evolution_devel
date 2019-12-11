#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "LATfield2.hpp"
#include "metadata.hpp"
#include "parser.hpp"
#include "chealpix.h"

// mpic++ pixelise_lc.cpp -o pixelise_lc -std=c++11 -O3 -DLIGHTCONE_THICKNESS=3 -DCOLORTERMINAL -DSINGLE -DHDF5 -I../LATfield2 -I/astro/adamek/local/include -I/common/users/timis/anaconda3/include -L/astro/adamek/local/lib -L/common/users/timis/anaconda3/lib -lhdf5 -lgsl -lgslcblas -lchealpix -lcfitsio

#define PIXBUF 4194304

using namespace std;

using namespace LATfield2;

struct healpix_header
{
	uint32_t Nside;
	uint32_t Npix;
	uint32_t precision;
	uint32_t Ngrid;
	double direction[3];
	double distance;
	double boxsize;
	uint32_t Nside_ring;
	char fill[256 - 5 * 4 - 5 * 8]; /* fills to 256 Bytes */
};

bool pointInShell(double * pos, lightcone_geometry & lightcone, double & outer, double & inner, double * vertex, int * wrap)
{
	double d;

	d = sqrt((pos[0]-(vertex[0]-wrap[0]))*(pos[0]-(vertex[0]-wrap[0])) + (pos[1]-(vertex[1]-wrap[1]))*(pos[1]-(vertex[1]-wrap[1])) + (pos[2]-(vertex[2]-wrap[2]))*(pos[2]-(vertex[2]-wrap[2])));

	if (d < inner || d >= outer) return false;

	if (lightcone.opening > -1.)
	{
#ifdef ACOS_HACK
		if (acos(((pos[0]-(vertex[0]-wrap[0]))*lightcone.direction[0] + (pos[1]-(vertex[1]-wrap[1]))*lightcone.direction[1] + (pos[2]-(vertex[2]-wrap[2]))*lightcone.direction[2]) / d) < acos(lightcone.opening)) return true;
#else
		if (((pos[0]-(vertex[0]-wrap[0]))*lightcone.direction[0] + (pos[1]-(vertex[1]-wrap[1]))*lightcone.direction[1] + (pos[2]-(vertex[2]-wrap[2]))*lightcone.direction[2]) / d > lightcone.opening) return true;
#endif
		else return false;
	}
	else return true;
}

int main(int argc, char **argv)
{
	H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

	char * settingsfile = NULL;
	char * cycleparam = NULL;
	char * lightconeparam = NULL;
	char * fieldparam = NULL;
	char * Nsideparam = NULL;
	int fieldselector = 0;
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
	double * di_info[MAX_OUTPUTS];
	double * d01_info[MAX_OUTPUTS];
	double * d12_info[MAX_OUTPUTS];
	double * do_info[MAX_OUTPUTS];

	double dummy_read[LIGHTCONE_THICKNESS+5];
	char filename[1024];
	FILE * infile;
	FILE * outfile;
	size_t num_read1, num_read2;
	uint32_t blocksize;

	double * vertex = NULL;
	double vertexlist[MAX_INTERSECTS][3];
	double domain[6] = {0., 0., 0., 1., 1., 1.};
	int numvertex;
	int tryvertex = 0;
	double z_obs = -2;

	int Nside_min = 64;
	int Nside_max = 64;
	int Nside;
	long Npix;
	long batch;
	long p;
	Real * pix = NULL;
	healpix_header hdr;
	int min_dist;
	int max_dist;
	int dist;
	double temp;
	double * direction;
	double pos[3];
	double w[3];
	int base_pos[3];
	int wrap[3];
	double R[3][3];
	const Real Healpix_undef = -1.6375e30;

	int n = 0, m = 0;

	for (int i = 1; i < argc; i++ ){
		if ( argv[i][0] != '-' )
			continue;
		switch(argv[i][1]) {
			case 's':
				settingsfile = argv[++i]; //settings file name
				break;
			case 'l':
				lightconeparam = argv[++i]; // light cone selector
				break;
			case 'c':
				cycleparam = argv[++i]; // cycle range selector
				break;
			case 'f':
				fieldparam = argv[++i]; // field selector
				break;
			case 'N':
				Nsideparam = argv[++i]; // field selector
				break;
			case 'n':
				n = atoi(argv[++i]); //size of the dim 1 of the processor grid
				break;
			case 'm':
				m =  atoi(argv[++i]); //size of the dim 2 of the processor grid
		}
	}

	parallel.initialize(n,m);

	if (settingsfile == NULL)
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no settings file specified!" << endl;
		return -1;
	}
	else if (fieldparam == NULL)
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no field requested!" << endl;
		return -1;
	}

	if (strcmp(fieldparam, "phi") == 0)
	{
		fieldselector = MASK_PHI;
	}
	else if (strcmp(fieldparam, "chi") == 0)
	{
		fieldselector = MASK_CHI;
	}
	else if (strcmp(fieldparam, "B1") == 0 || strcmp(fieldparam, "B2") == 0 || strcmp(fieldparam, "B3") == 0)
	{
		fieldselector = MASK_B;
	}
	else
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": field request (parameter -f) not recognised!" << endl;
		return -1;
	}

	COUT << COLORTEXT_WHITE << " LCARS tools: pixelise light cone" << endl << COLORTEXT_RESET << " compiled with " << COLORTEXT_RED << "LIGHTCONE_THICKNESS=" << LIGHTCONE_THICKNESS << endl << COLORTEXT_RESET << endl << " opening settings file of simulation: " << settingsfile << endl << " parser output:" << endl << endl;

	numparam = loadParameterFile(settingsfile, params);
	usedparams = parseMetadata(params, numparam, sim, cosmo, ic);
	free(params);

	COUT << endl << " file contains " << numparam << " parameters, " << usedparams << " of which could be parsed." << endl << endl;

	COUT << " number of lightcones: " << sim.num_lightcone << endl << endl;

	for (int i = 0; i < sim.num_lightcone; i++)
	{
		COUT << " lightcone " << i << " parameters:" << endl << "  vertex = " << sim.lightcone[i].vertex[0] << ", " << sim.lightcone[i].vertex[1] << ", " << sim.lightcone[i].vertex[2] << endl;
		COUT << "  redshift of observation = " << sim.lightcone[i].z << endl;
		COUT << "  direction = " << sim.lightcone[i].direction[0] << ", " << sim.lightcone[i].direction[1] << ", " << sim.lightcone[i].direction[2] << endl;
		COUT << "  opening half-angle = " << /* sim.lightcone[i].opening * 180. / M_PI */ ((sim.lightcone[i].opening > -1.) ? acos(sim.lightcone[i].opening) * 180. / M_PI : 180.) << " degrees" << endl;
		COUT << "  distance interval = " << sim.lightcone[i].distance[0] << ", " << sim.lightcone[i].distance[1] << endl << endl;
	}

	if (Nsideparam == NULL)
	{
		COUT << " no Nside selected (parameter -N), default value (64)." << endl;
	}
	else
	{
		sprintf(filename, "%s", Nsideparam);
		char * dash = strchr(filename, '-');

		if (dash == NULL)
		{
			Nside_min = (Nside_max = atoi(Nsideparam));
		}
		else
		{
			Nside_max = atoi(dash+1);
			dash[0] = '\0';
			Nside_min = atoi(filename);
		}

		COUT << " Nside set to: " << Nside_min;
		if (Nside_max != Nside_min)
		{
			COUT << "-" << Nside_max;
		}
		COUT << endl;
	}

	if (lightconeparam == NULL)
	{
		COUT << " no light cones selected (parameter -l), using all light cones with requested field information: ";
		for (int i = 0; i < sim.num_lightcone; i++)
		{
			if (sim.out_lightcone[i] & fieldselector)
			{
				COUT << i << " ";
				use_lightcone[i] = true;
			}
			else use_lightcone[i] = false;
		}
		COUT << endl << endl;
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
				COUT << " light cone " << i << " (parameter -l) does not exist!" << endl;
			}
			else if (!(sim.out_lightcone[i] & fieldselector))
			{
				COUT << " light cone " << i << " does not contain required field information and will be omitted!" << endl;
			}
			else
			{
				use_lightcone[i] = true;
			}

			token = strtok(NULL, ",");
		}

		COUT << " using light cone(s) ";
		for (int i = 0; i < sim.num_lightcone; i++)
		{
			if (use_lightcone[i]) COUT << i << " ";
		}
		COUT << endl << endl;
	}

	COUT << " reading light cone information files..." << endl;

	for (int i = 0; i < sim.num_lightcone; i++)
	{
		if (use_lightcone[i])
		{
			if (z_obs < -1.)
			{
				vertex = sim.lightcone[i].vertex;
				z_obs = sim.lightcone[i].z;
				direction = sim.lightcone[i].direction;
			}
			else if (sim.lightcone[i].z != z_obs || sim.lightcone[i].vertex[0] != vertex[0] || sim.lightcone[i].vertex[1] != vertex[1] || sim.lightcone[i].vertex[2] != vertex[2])
			{
				COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": not all selected light cones refer to the same observation event!" << endl << " incompatible lightcone " << i << " will be omitted!" << endl << endl;
				use_lightcone[i] = false;
				continue;
			}

			if (sim.lightcone[i].opening > -0.999) direction = sim.lightcone[i].direction;

			if (parallel.isRoot())
			{
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
			}

			parallel.broadcast(num_shells[i], 0);

			di_info[i] = (double *) malloc(num_shells[i] * sizeof(double));
			d01_info[i] = (double *) malloc(num_shells[i] * sizeof(double));
			d12_info[i] = (double *) malloc(num_shells[i] * sizeof(double));
			do_info[i] = (double *) malloc(num_shells[i] * sizeof(double));

			if (parallel.isRoot())
			{
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
	
					di_info[i][j] = dummy_read[4];
					d01_info[i][j] = dummy_read[5];
					d12_info[i][j] = dummy_read[6];
					do_info[i][j] = dummy_read[7];
				}

				fclose(infile);
			}

			parallel.broadcast(di_info[i], num_shells[i], 0);
			parallel.broadcast(d01_info[i], num_shells[i], 0);
			parallel.broadcast(d12_info[i], num_shells[i], 0);
			parallel.broadcast(do_info[i], num_shells[i], 0);

			parallel.broadcast(start_cycle[i], 0);
			parallel.broadcast(end_cycle[i], 0);

			if (start_cycle[i] < min_cycle) min_cycle = start_cycle[i];
			if (end_cycle[i] > max_cycle) max_cycle = end_cycle[i];

			COUT << " light cone " << i << " contains " << num_shells[i] << " shells (cycle " << start_cycle[i] << " to " << end_cycle[i] << ")" << endl;
		}
	}

	COUT << endl;

	if (cycleparam == NULL)
	{
		COUT << " no range of cycles selected (parameter -c), using entire range: " << min_cycle << "-" << max_cycle << endl;
	}
	else
	{
		sprintf(filename, "%s", cycleparam);
		char * dash = strchr(filename, '-');
		int tmp;

		if (dash == NULL)
		{
			COUT << " range of cycle (parameter -c) could not be interpreted!" << endl;

			for (int i = 0; i < sim.num_lightcone; i++)
			{
				if (use_lightcone[i])
				{
					free(di_info[i]);
					free(d01_info[i]);
					free(d12_info[i]);
					free(do_info[i]);
				}
			}

			return -1;
		}

		tmp = atoi(dash+1);

		if (tmp > max_cycle)
		{
			COUT << " range of cycles (parameter -c) extends beyond available range; last cycle in range is " << max_cycle << endl;
		}
		else if (tmp < min_cycle)
		{
			COUT << " range of cycles (parameter -c) outside available range!" << endl;

			for (int i = 0; i < sim.num_lightcone; i++)
			{
				if (use_lightcone[i])
				{
					free(di_info[i]);
					free(d01_info[i]);
					free(d12_info[i]);
					free(do_info[i]);
				}
			}

			return -1;
		}
		else max_cycle = tmp;

		dash[0] = '\0';
		tmp = atoi(filename);

		if (tmp < min_cycle)
		{
			COUT << " range of cycles (parameter -c) extends beyond available range; first cycle in range is " << min_cycle << endl;
		}
		else if (tmp > max_cycle)
		{
			COUT << " range of cycles (parameter -c) outside available range!" << endl;

			for (int i = 0; i < sim.num_lightcone; i++)
			{
				if (use_lightcone[i])
				{
					free(di_info[i]);
					free(d01_info[i]);
					free(d12_info[i]);
					free(do_info[i]);
				}
			}

			return -1;
		}
		else min_cycle = tmp;

		COUT << " range of cycles set to: " << min_cycle << "-" << max_cycle << endl;
	}

	COUT << endl;

	if (direction[0] == 0 && direction[1] == 0)
	{
		R[0][0] = direction[2];
		R[0][1] = 0;
		R[0][2] = 0;
		R[1][0] = 0;
		R[1][1] = 1;
		R[1][2] = 0;
		R[2][0] = 0;
		R[2][1] = 0;
		R[2][2] = direction[2];
	}
	else
	{
		temp = atan2(direction[1], direction[0]);
		R[0][0] = cos(temp) * direction[2];
		R[0][1] = -sin(temp);
		R[0][2] = direction[0];
		R[1][0] = sin(temp) * direction[2];
		R[1][1] = cos(temp);
		R[1][2] = direction[1];
		R[2][0] = -sqrt(1. - direction[2]*direction[2]);
		R[2][1] = 0;
		R[2][2] = direction[2];
	}

	int dim = 3;
	int latSize[3] = {sim.numpts,sim.numpts,sim.numpts};
	int halo = 2;
	Lattice lat(dim,latSize,halo);

	Field<Real> *field[3];

	Site x(lat);

	field[0] = new Field<Real>;
	field[1] = new Field<Real>(lat, 1);
	field[2] = new Field<Real>(lat, 1);

	field[0]->initialize(lat, 1);
	field[0]->alloc();
	field[1]->initialize(lat, 1);
	field[1]->alloc();
	field[2]->initialize(lat, 1);
	field[2]->alloc();

	pix = (Real *) malloc(PIXBUF * sizeof(Real));

#ifdef SINGLE
	if (parallel.isRoot())
		MPI_Reduce(MPI_IN_PLACE, pix, PIXBUF, MPI_FLOAT, MPI_SUM, parallel.root(), parallel.lat_world_comm());
	else
		MPI_Reduce(pix, NULL, PIXBUF, MPI_FLOAT, MPI_SUM, parallel.root(), parallel.lat_world_comm());
#else
	if (parallel.isRoot())
		MPI_Reduce(MPI_IN_PLACE, pix, PIXBUF, MPI_DOUBLE, MPI_SUM, parallel.root(), parallel.lat_world_comm());
	else
		MPI_Reduce(pix, NULL, PIXBUF, MPI_DOUBLE, MPI_SUM, parallel.root(), parallel.lat_world_comm());
#endif

	string h5filename(sim.output_path);

	COUT << " loading HDF5 files..." << endl;

	COUT << " " << h5filename << sim.basename_lightcone << "_" << fieldparam << "_0.h5 ...";
	field[0]->loadHDF5(h5filename+sim.basename_lightcone+"_"+fieldparam+"_0.h5");
	field[0]->updateHalo();
	COUT << " OK" << endl << " " << h5filename << sim.basename_lightcone << "_" << fieldparam << "_1.h5 ...";
	field[1]->loadHDF5(h5filename+sim.basename_lightcone+"_"+fieldparam+"_1.h5");
	field[1]->updateHalo();
	COUT << " OK" << endl << " " << h5filename << sim.basename_lightcone << "_" << fieldparam << "_2.h5 ...";
	field[2]->loadHDF5(h5filename+sim.basename_lightcone+"_"+fieldparam+"_2.h5");
	field[2]->updateHalo();
	COUT << " OK" << endl << endl;
	
#ifdef LIGHTCONE_HOTFIX
	base_pos[0] = 2793;
	base_pos[1] = 1843;
	base_pos[2] = 3800;
	
	if (x.setCoord(base_pos))
	{
		cout << " proc#" << parallel.rank() << ": fixing broken light cone at (2793,1843,3800), [" << (*field[0])(x) << ", " << (*field[1])(x) << ", " << (*field[2])(x) << "] -> [";
		(*field[0])(x) = ((*field[0])(x+0) + (*field[0])(x-0) + (*field[0])(x+1) + (*field[0])(x-1) + (*field[0])(x+2) + (*field[0])(x-2)) / 6.;
		cout << (*field[0])(x) << ", ";
		(*field[1])(x) = ((*field[1])(x+0) + (*field[1])(x-0) + (*field[1])(x+1) + (*field[1])(x-1) + (*field[1])(x+2) + (*field[1])(x-2)) / 6.;
		cout << (*field[1])(x) << ", ";
		(*field[2])(x) = ((*field[2])(x+0) + (*field[2])(x-0) + (*field[2])(x+1) + (*field[2])(x-1) + (*field[2])(x+2) + (*field[2])(x-2)) / 6.;
		cout << (*field[2])(x) << "]" << endl;
	}
	
	parallel.barrier();
	
	base_pos[0] = 5586;
	base_pos[1] = 3686;
	base_pos[2] = 7600;
	
	if (x.setCoord(base_pos))
	{
		cout << " proc#" << parallel.rank() << ": fixing broken light cone at (5586,3686,7600), [" << (*field[0])(x) << ", " << (*field[1])(x) << ", " << (*field[2])(x) << "] -> [";
		(*field[0])(x) = ((*field[0])(x+0) + (*field[0])(x-0) + (*field[0])(x+1) + (*field[0])(x-1) + (*field[0])(x+2) + (*field[0])(x-2)) / 6.;
		cout << (*field[0])(x) << ", ";
		(*field[1])(x) = ((*field[1])(x+0) + (*field[1])(x-0) + (*field[1])(x+1) + (*field[1])(x-1) + (*field[1])(x+2) + (*field[1])(x-2)) / 6.;
		cout << (*field[1])(x) << ", ";
		(*field[2])(x) = ((*field[2])(x+0) + (*field[2])(x-0) + (*field[2])(x+1) + (*field[2])(x-1) + (*field[2])(x+2) + (*field[2])(x-2)) / 6.;
		cout << (*field[2])(x) << "]" << endl << endl;
	}
	
	field[0]->updateHalo();
	field[1]->updateHalo();
	field[2]->updateHalo();
#endif

	COUT << " pixelising the light cone..." << endl << endl;

	for (int c = min_cycle; c <= max_cycle; c++)
	{
		min_dist = -1;
		max_dist = -1;

		COUT << " cycle " << c << " is being processed..." << endl;

		for (int i = 0; i < sim.num_lightcone; i++)
		{
			if (!use_lightcone[i] || c < start_cycle[i] || c > end_cycle[i]) continue;

			if (d01_info[i][c-start_cycle[i]] >= sim.lightcone[i].distance[1])
				temp = (c < end_cycle[i]) ? d01_info[i][c+1-start_cycle[i]] : di_info[i][c-start_cycle[i]];
			else if (d12_info[i][c-start_cycle[i]] >= sim.lightcone[i].distance[1])
				temp = d01_info[i][c-start_cycle[i]];
			else if (do_info[i][c-start_cycle[i]] >= sim.lightcone[i].distance[1])
				temp = d12_info[i][c-start_cycle[i]];
			else continue;

			if (max_dist > 0)
				min_dist = max_dist+1;
			else if (temp > 0)
				min_dist = ceil((temp * sim.numpts) + 1.7321);
			else
				min_dist = 0;

			if (do_info[i][c-start_cycle[i]] < sim.lightcone[i].distance[0])
				temp = do_info[i][c-start_cycle[i]];
			else if (d12_info[i][c-start_cycle[i]] < sim.lightcone[i].distance[0])
				temp = d12_info[i][c-start_cycle[i]];
			else if (d01_info[i][c-start_cycle[i]] < sim.lightcone[i].distance[0])
				temp = d01_info[i][c-start_cycle[i]];
			else continue;
			
			max_dist = floor((temp * sim.numpts) - 1.7321);

			COUT << " stepping through distance " << min_dist << " to " << max_dist << " on light cone " << i << endl;

			for (dist = min_dist; dist <= max_dist; dist++)
			{
				for (Nside = Nside_min; Nside < Nside_max; Nside *= 2)
				{
					if (12. * (double) Nside * (double) Nside > 4. * M_PI * (double) dist * (double) dist) break;
				}

				Npix = 4;
				for (int ring = 2; Npix < (1. - sim.lightcone[i].opening) * 6 * Nside * Nside; ring++)
				{
					if (ring <= Nside) Npix += 4 * ring;
					else if (ring <= 3 * Nside - 1) Npix += 4 * Nside;
					else Npix += 4 * (4 * Nside - ring);
				}

				COUT << " shell at distance " << dist << " on light cone " << i << ": Nside=" << Nside << ", Npix=" << Npix << endl;
				
				if (parallel.isRoot())
				{
					sprintf(filename, "%s%s_%04d_%s.map", sim.output_path, sim.basename_lightcone, c, fieldparam);

					outfile = fopen(filename, "ab");

					if (outfile == NULL)
					{
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to open file " << filename << " for output!" << endl;
					}
					else
					{
						hdr.Nside = Nside;
						hdr.Npix = Npix;
						hdr.precision = sizeof(Real);
						hdr.Ngrid = sim.numpts;
						hdr.direction[0] = direction[0];
						hdr.direction[1] = direction[1];
						hdr.direction[2] = direction[2];
						hdr.distance = (double) dist / (double) sim.numpts;
						hdr.boxsize = sim.boxsize;
						hdr.Nside_ring = Nside;

						blocksize = 256;
						
						fwrite(&blocksize, sizeof(uint32_t), 1, outfile);
						fwrite(&hdr, sizeof(hdr), 1, outfile);
						fwrite(&blocksize, sizeof(uint32_t), 1, outfile);

						blocksize = Npix * sizeof(Real);

						fwrite(&blocksize, sizeof(uint32_t), 1, outfile);
							
						fclose(outfile);
					}
				}

				batch = (Npix > PIXBUF) ? PIXBUF : Npix;
				//pix = (Real *) malloc(batch * sizeof(Real));

				for (int64_t q = 0; q < Npix; q++)
				{
					p = q % batch;
					
					if (q > 0 && p == 0)
					{
						//parallel.sum(pix, batch);
						
						if (parallel.isRoot())
						{
#ifdef SINGLE
							MPI_Reduce(MPI_IN_PLACE, pix, batch, MPI_FLOAT, MPI_SUM, parallel.root(), parallel.lat_world_comm());
#else
							MPI_Reduce(MPI_IN_PLACE, pix, batch, MPI_DOUBLE, MPI_SUM, parallel.root(), parallel.lat_world_comm());
#endif

							sprintf(filename, "%s%s_%04d_%s.map", sim.output_path, sim.basename_lightcone, c, fieldparam);

							outfile = fopen(filename, "ab");

							if (outfile == NULL)
							{
								cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to open file " << filename << " for output!" << endl;
							}
							else
							{
								fwrite(pix, sizeof(Real), batch, outfile);
							
								fclose(outfile);
							}
						}
						else
#ifdef SINGLE
							MPI_Reduce(pix, NULL, batch, MPI_FLOAT, MPI_SUM, parallel.root(), parallel.lat_world_comm());
#else
							MPI_Reduce(pix, NULL, batch, MPI_DOUBLE, MPI_SUM, parallel.root(), parallel.lat_world_comm());
#endif
					}
					
					pix2vec_ring64(Nside, q, w);
					
					pos[0] = dist * (R[0][0] * w[0] + R[0][1] * w[1] + R[0][2] * w[2]) + vertex[0] * sim.numpts;
					pos[1] = dist * (R[1][0] * w[0] + R[1][1] * w[1] + R[1][2] * w[2]) + vertex[1] * sim.numpts;
					pos[2] = dist * (R[2][0] * w[0] + R[2][1] * w[1] + R[2][2] * w[2]) + vertex[2] * sim.numpts;

					wrap[0] = floor(pos[0]/sim.numpts);
					wrap[1] = floor(pos[1]/sim.numpts);
					wrap[2] = floor(pos[2]/sim.numpts);

					w[0] = modf(pos[0] - wrap[0] * sim.numpts, &temp);
					base_pos[0] = (int) temp;
					w[1] = modf(pos[1] - wrap[1] * sim.numpts, &temp);
					base_pos[1] = (int) temp;
					w[2] = modf(pos[2] - wrap[2] * sim.numpts, &temp);
					base_pos[2] = (int) temp;

					if (x.setCoord(base_pos))
					{
						pos[0] = ((double) base_pos[0]) / ((double) sim.numpts);
						pos[1] = ((double) base_pos[1]) / ((double) sim.numpts);
						pos[2] = ((double) base_pos[2]) / ((double) sim.numpts);

						if (fieldselector == MASK_B)
						{
							pos[(int) (fieldparam[1]-'1')] = (0.5 + (double) base_pos[(int) (fieldparam[1]-'1')]) / ((double) sim.numpts);
							if (w[(int) (fieldparam[1]-'1')] >= 0.5)
							{
								w[(int) (fieldparam[1]-'1')] -= 0.5;
								temp = 1;
							}
							else
							{
								w[(int) (fieldparam[1]-'1')] = 0.5 - w[(int) (fieldparam[1]-'1')];
								temp = 0;
							}
						}

						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
#ifdef TSP_INTERPOLATION
							pix[p] = (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[0])(x);
#else
							pix[p] = (1.-w[0]) * (1.-w[1]) * (1.-w[2]) * (*field[0])(x);
#endif
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
#ifdef TSP_INTERPOLATION
							pix[p] = (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[1])(x);
#else
							pix[p] = (1.-w[0]) * (1.-w[1]) * (1.-w[2]) * (*field[1])(x);
#endif
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
#ifdef TSP_INTERPOLATION
							pix[p] = (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[2])(x);
#else
							pix[p] = (1.-w[0]) * (1.-w[1]) * (1.-w[2]) * (*field[2])(x);
#endif
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}

						if (fieldselector == MASK_B && fieldparam[1] == '1')
						{
							if (temp)
							{
								pos[0] = (0.5 + (double) (base_pos[0]+1)) / ((double) sim.numpts);
								x = x+0;
							}
							else
							{
								pos[0] = (0.5 + (double) (base_pos[0]-1)) / ((double) sim.numpts);
								x = x-0;
							}
						}
						else
						{
							pos[0] = ((double) (base_pos[0]+1)) / ((double) sim.numpts);
							x = x+0;
						}

						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
#ifdef TSP_INTERPOLATION
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[0])(x);
#else
							pix[p] += w[0] * (1.-w[1]) * (1.-w[2]) * (*field[0])(x);
#endif
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
#ifdef TSP_INTERPOLATION
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[1])(x);
#else
							pix[p] += w[0] * (1.-w[1]) * (1.-w[2]) * (*field[1])(x);
#endif
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
#ifdef TSP_INTERPOLATION
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[2])(x);
#else
							pix[p] += w[0] * (1.-w[1]) * (1.-w[2]) * (*field[2])(x);
#endif
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}

						if (fieldselector == MASK_B && fieldparam[1] == '2')
						{
							if (temp)
							{
								pos[1] = (0.5 + (double) (base_pos[1]+1)) / ((double) sim.numpts);
								x = x+1;
							}
							else
							{
								pos[1] = (0.5 + (double) (base_pos[1]-1)) / ((double) sim.numpts);
								x = x-1;
							}
						}
						else
						{
							pos[1] = ((double) (base_pos[1]+1)) / ((double) sim.numpts);
							x = x+1;
						}

						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
#ifdef TSP_INTERPOLATION
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[0])(x);
#else
							pix[p] += w[0] * w[1] * (1.-w[2]) * (*field[0])(x);
#endif
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
#ifdef TSP_INTERPOLATION
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[1])(x);
#else
							pix[p] += w[0] * w[1] * (1.-w[2]) * (*field[1])(x);
#endif
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
#ifdef TSP_INTERPOLATION
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[2])(x);
#else
							pix[p] += w[0] * w[1] * (1.-w[2]) * (*field[2])(x);
#endif
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}

						if (fieldselector == MASK_B && fieldparam[1] == '1')
						{
							pos[0] = (0.5 + (double) (base_pos[0])) / ((double) sim.numpts);
							if (temp)
								x = x-0;
							else
								x = x+0;
						}
						else
						{
							pos[0] = ((double) base_pos[0]) / ((double) sim.numpts);
							x = x-0;
						}

						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
#ifdef TSP_INTERPOLATION
							pix[p] += (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[0])(x);
#else
							pix[p] += (1.-w[0]) * w[1] * (1.-w[2]) * (*field[0])(x);
#endif
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
#ifdef TSP_INTERPOLATION
							pix[p] += (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[1])(x);
#else
							pix[p] += (1.-w[0]) * w[1] * (1.-w[2]) * (*field[1])(x);
#endif
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
#ifdef TSP_INTERPOLATION
							pix[p] += (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[2])(x);
#else
							pix[p] += (1.-w[0]) * w[1] * (1.-w[2]) * (*field[2])(x);
#endif
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}

						if (fieldselector == MASK_B && fieldparam[1] == '3')
						{
							if (temp)
							{
								pos[2] = (0.5 + (double) (base_pos[2]+1)) / ((double) sim.numpts);
								x = x+2;
							}
							else
							{
								pos[2] = (0.5 + (double) (base_pos[2]-1)) / ((double) sim.numpts);
								x = x-2;
							}
						}
						else
						{
							pos[2] = ((double) (base_pos[2]+1)) / ((double) sim.numpts);
							x = x+2;
						}

						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
#ifdef TSP_INTERPOLATION
							pix[p] += (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[0])(x);
#else
							pix[p] += (1.-w[0]) * w[1] * w[2] * (*field[0])(x);
#endif
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
#ifdef TSP_INTERPOLATION
							pix[p] += (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[1])(x);
#else
							pix[p] += (1.-w[0]) * w[1] * w[2] * (*field[1])(x);
#endif
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
#ifdef TSP_INTERPOLATION
							pix[p] += (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[2])(x);
#else
							pix[p] += (1.-w[0]) * w[1] * w[2] * (*field[2])(x);
#endif
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}

						if (fieldselector == MASK_B && fieldparam[1] == '1')
						{
							if (temp)
							{
								pos[0] = (0.5 + (double) (base_pos[0]+1)) / ((double) sim.numpts);
								x = x+0;
							}
							else
							{
								pos[0] = (0.5 + (double) (base_pos[0]-1)) / ((double) sim.numpts);
								x = x-0;
							}
						}
						else
						{
							pos[0] = ((double) (base_pos[0]+1)) / ((double) sim.numpts);
							x = x+0;
						}

						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
#ifdef TSP_INTERPOLATION
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[0])(x);
#else
							pix[p] += w[0] * w[1] * w[2] * (*field[0])(x);
#endif
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
#ifdef TSP_INTERPOLATION
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[1])(x);
#else
							pix[p] += w[0] * w[1] * w[2] * (*field[1])(x);
#endif
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
#ifdef TSP_INTERPOLATION
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[2])(x);
#else
							pix[p] += w[0] * w[1] * w[2] * (*field[2])(x);
#endif
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}

						if (fieldselector == MASK_B && fieldparam[1] == '2')
						{
							pos[1] = (0.5 + (double) (base_pos[1])) / ((double) sim.numpts);
							if (temp)
								x = x-1;
							else
								x = x+1;
						}
						else
						{
							pos[1] = ((double) base_pos[1]) / ((double) sim.numpts);
							x = x-1;
						}

						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
#ifdef TSP_INTERPOLATION
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[0])(x);
#else
							pix[p] += w[0] * (1.-w[1]) * w[2] * (*field[0])(x);
#endif
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
#ifdef TSP_INTERPOLATION
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[1])(x);
#else
							pix[p] += w[0] * (1.-w[1]) * w[2] * (*field[1])(x);
#endif
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
#ifdef TSP_INTERPOLATION
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[2])(x);
#else
							pix[p] += w[0] * (1.-w[1]) * w[2] * (*field[2])(x);
#endif
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}

						if (fieldselector == MASK_B && fieldparam[1] == '1')
						{
							pos[0] = (0.5 + (double) (base_pos[0])) / ((double) sim.numpts);
							if (temp)
								x = x-0;
							else
								x = x+0;
						}
						else
						{
							pos[0] = ((double) base_pos[0]) / ((double) sim.numpts);
							x = x-0;
						}

						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
#ifdef TSP_INTERPOLATION
							pix[p] += (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[0])(x);
#else
							pix[p] += (1.-w[0]) * (1.-w[1]) * w[2] * (*field[0])(x);
#endif
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
#ifdef TSP_INTERPOLATION
							pix[p] += (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[1])(x);
#else
							pix[p] += (1.-w[0]) * (1.-w[1]) * w[2] * (*field[1])(x);
#endif
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
#ifdef TSP_INTERPOLATION
							pix[p] += (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[2])(x);
#else
							pix[p] += (1.-w[0]) * (1.-w[1]) * w[2] * (*field[2])(x);
#endif
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}
						
#ifdef TSP_INTERPOLATION
						if (w[2] > 0.5)
						{
							pos[2] = ((double) (base_pos[2]+2)) / ((double) sim.numpts);
							x = x+2;
						}
						else
						{
							pos[2] = ((double) (base_pos[2]-1)) / ((double) sim.numpts);
							x = x-2;
							x = x-2;
						}
							
						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[0])(x);
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[1])(x);
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[2])(x);
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}
							
						pos[0] = ((double) (base_pos[0]+1)) / ((double) sim.numpts);
						x = x+0;
							
						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[0])(x);
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[1])(x);
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[2])(x);
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}
							
						pos[1] = ((double) (base_pos[1]+1)) / ((double) sim.numpts);
						x = x+1;
							
						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[0])(x);
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[1])(x);
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[2])(x);
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}
							
						pos[0] = ((double) (base_pos[0])) / ((double) sim.numpts);
						x = x-0;
							
						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[0])(x);
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[1])(x);
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[2])(x);
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}
							
						if (w[1] > 0.5)
						{
							pos[1] = ((double) (base_pos[1]+2)) / ((double) sim.numpts);
							x = x+1;
						}
						else
						{
							pos[1] = ((double) (base_pos[1]-1)) / ((double) sim.numpts);
							x = x-1;
							x = x-1;
						}
						
							
						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * ((w[1]*(w[1]-1.)+0.25)/2.) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[0])(x);
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * ((w[1]*(w[1]-1.)+0.25)/2.) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[1])(x);
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * ((w[1]*(w[1]-1.)+0.25)/2.) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[2])(x);
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}
								
						pos[0] = ((double) (base_pos[0]+1)) / ((double) sim.numpts);
						x = x+0;
							
						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * ((w[1]*(w[1]-1.)+0.25)/2.) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[0])(x);
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * ((w[1]*(w[1]-1.)+0.25)/2.) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[1])(x);
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * ((w[1]*(w[1]-1.)+0.25)/2.) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[2])(x);
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}
							
						if (w[0] > 0.5)
						{
							pos[0] = ((double) (base_pos[0]+2)) / ((double) sim.numpts);
							x = x+0;
						}
						else
						{
							pos[0] = ((double) (base_pos[0]-1)) / ((double) sim.numpts);
							x = x-0;
							x = x-0;
						}
								
						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * ((w[1]*(w[1]-1.)+0.25)/2.) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[0])(x);
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * ((w[1]*(w[1]-1.)+0.25)/2.) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[1])(x);
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * ((w[1]*(w[1]-1.)+0.25)/2.) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[2])(x);
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}
						
						pos[1] = ((double) (base_pos[1])) / ((double) sim.numpts);
						if (w[1] > 0.5)
						{
							x = x-1;
							x = x-1;
						}
						else
							x = x+1;
							
						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[0])(x);
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[1])(x);
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[2])(x);
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}
						
						pos[1] = ((double) (base_pos[1]+1)) / ((double) sim.numpts);
						x = x+1;
						
						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[0])(x);
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[1])(x);
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * ((w[2]*(w[2]-1.)+0.25)/2.) * (*field[2])(x);
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}
						
						pos[2] = ((double) (base_pos[2])) / ((double) sim.numpts);
						if (w[2] > 0.5)
						{
							x = x-2;
							x = x-2;
						}
						else
							x = x+2;
							
						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[0])(x);
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[1])(x);
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[2])(x);
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}
						
						pos[1] = ((double) (base_pos[1])) / ((double) sim.numpts);
						x = x-1;
						
						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[0])(x);
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[1])(x);
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[2])(x);
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}
						
						pos[2] = ((double) (base_pos[2]+1)) / ((double) sim.numpts);
						x = x+2;
						
						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[0])(x);
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[1])(x);
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * (w[1]>0.5 ? (2.25+w[1]*(w[1]-3.))/2. : 0.75-w[1]*w[1]) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[2])(x);
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}
						
						pos[1] = ((double) (base_pos[1]+1)) / ((double) sim.numpts);
						x = x+1;
						
						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[0])(x);
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[1])(x);
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * (w[1]>0.5 ? w[1]*(2.-w[1])-0.25 : (w[1]*(w[1]+1.)+0.25)/2.) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[2])(x);
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}
						
						if (w[1] > 0.5)
						{
							pos[1] = ((double) (base_pos[1]+2)) / ((double) sim.numpts);
							x = x+1;
						}
						else
						{
							pos[1] = ((double) (base_pos[1]-1)) / ((double) sim.numpts);
							x = x-1;
							x = x-1;
						}
						
						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * ((w[1]*(w[1]-1.)+0.25)/2.) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[0])(x);
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * ((w[1]*(w[1]-1.)+0.25)/2.) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[1])(x);
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * ((w[1]*(w[1]-1.)+0.25)/2.) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[2])(x);
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}
						
						pos[2] = ((double) (base_pos[2])) / ((double) sim.numpts);
						x = x-2;
						
						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * ((w[1]*(w[1]-1.)+0.25)/2.) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[0])(x);
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * ((w[1]*(w[1]-1.)+0.25)/2.) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[1])(x);
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += ((w[0]*(w[0]-1.)+0.25)/2.) * ((w[1]*(w[1]-1.)+0.25)/2.) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[2])(x);
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}
						
						pos[0] = ((double) (base_pos[0])) / ((double) sim.numpts);
						if (w[0] > 0.5)
						{
							x = x-0;
							x = x-0;
						}
						else
							x = x+0;
							
						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * ((w[1]*(w[1]-1.)+0.25)/2.) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[0])(x);
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * ((w[1]*(w[1]-1.)+0.25)/2.) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[1])(x);
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * ((w[1]*(w[1]-1.)+0.25)/2.) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[2])(x);
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}
						
						pos[0] = ((double) (base_pos[0]+1)) / ((double) sim.numpts);
						x = x+0;
						
						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * ((w[1]*(w[1]-1.)+0.25)/2.) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[0])(x);
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * ((w[1]*(w[1]-1.)+0.25)/2.) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[1])(x);
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * ((w[1]*(w[1]-1.)+0.25)/2.) * (w[2]>0.5 ? (2.25+w[2]*(w[2]-3.))/2. : 0.75-w[2]*w[2]) * (*field[2])(x);
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}
						
						pos[2] = ((double) (base_pos[2]+1)) / ((double) sim.numpts);
						x = x+2;
						
						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * ((w[1]*(w[1]-1.)+0.25)/2.) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[0])(x);
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * ((w[1]*(w[1]-1.)+0.25)/2.) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[1])(x);
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? w[0]*(2.-w[0])-0.25 : (w[0]*(w[0]+1.)+0.25)/2.) * ((w[1]*(w[1]-1.)+0.25)/2.) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[2])(x);
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}
						
						pos[0] = ((double) (base_pos[0])) / ((double) sim.numpts);
						x = x-0;
						
						if (pointInShell(pos, sim.lightcone[i], d01_info[i][c-start_cycle[i]], di_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * ((w[1]*(w[1]-1.)+0.25)/2.) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[0])(x);
						else if (pointInShell(pos, sim.lightcone[i], d12_info[i][c-start_cycle[i]], d01_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * ((w[1]*(w[1]-1.)+0.25)/2.) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[1])(x);
						else if (pointInShell(pos, sim.lightcone[i], do_info[i][c-start_cycle[i]], d12_info[i][c-start_cycle[i]], vertex, wrap))
							pix[p] += (w[0]>0.5 ? (2.25+w[0]*(w[0]-3.))/2. : 0.75-w[0]*w[0]) * ((w[1]*(w[1]-1.)+0.25)/2.) * (w[2]>0.5 ? w[2]*(2.-w[2])-0.25 : (w[2]*(w[2]+1.)+0.25)/2.) * (*field[2])(x);
						else
						{
							pix[p] = Healpix_undef;
							continue;
						}
#endif						
					}
					else pix[p] = 0;
				}
				
				if (Npix % batch != 0) batch = Npix % batch;
				
				//parallel.sum(pix, batch);
					
				if (parallel.isRoot())
				{
#ifdef SINGLE
					MPI_Reduce(MPI_IN_PLACE, pix, batch, MPI_FLOAT, MPI_SUM, parallel.root(), parallel.lat_world_comm());
#else
					MPI_Reduce(MPI_IN_PLACE, pix, batch, MPI_DOUBLE, MPI_SUM, parallel.root(), parallel.lat_world_comm());
#endif
					sprintf(filename, "%s%s_%04d_%s.map", sim.output_path, sim.basename_lightcone, c, fieldparam);

					outfile = fopen(filename, "ab");

					if (outfile == NULL)
					{
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to open file " << filename << " for output!" << endl;
					}
					else
					{
						fwrite(pix, sizeof(Real), batch, outfile);

						blocksize = Npix * sizeof(Real);
						fwrite(&blocksize, sizeof(uint32_t), 1, outfile);
							
						fclose(outfile);
					}
				}
				else
#ifdef SINGLE
					MPI_Reduce(pix, NULL, batch, MPI_FLOAT, MPI_SUM, parallel.root(), parallel.lat_world_comm());
#else
					MPI_Reduce(pix, NULL, batch, MPI_DOUBLE, MPI_SUM, parallel.root(), parallel.lat_world_comm());
#endif

				//free(pix);
			}
		}
	}

	free(pix);

	delete field[0];
	delete field[1];
	delete field[2];

	for (int i = 0; i < sim.num_lightcone; i++)
	{
		if (use_lightcone[i])
		{
			free(di_info[i]);
			free(d01_info[i]);
			free(d12_info[i]);
			free(do_info[i]);
		}
	}

	COUT << endl << " normal completion." << endl;

	return 0;
}
