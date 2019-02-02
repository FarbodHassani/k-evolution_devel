#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include "metadata.hpp"
#include "background.hpp"

// g++ sample_lc.cpp -o sample_lc -std=c++11 -O3 -lgsl -lgslcblas

using namespace std;

#define BATCHSIZE 4194304

int main(int argc, char **argv)
{
	char * filebase = NULL;
	char filename[1024];
	double frac = 0.;
	float offset = 0.;

	gadget2_header hdr;
	float * posbatch = NULL;
	float * velbatch = NULL;
	uint32_t batch;
	uint32_t blocksize;
	long backtrack;
	long blockoffset;
	FILE * infile = NULL;
	int nfile = 0;
	//cosmology cosmo;
	//double fourpiG;
	//double a, chi;

	gsl_rng * rng = gsl_rng_alloc(gsl_rng_ranlux);
	gsl_rng_set(rng, 0);

	for (int i = 1 ; i < argc ; i++ ){
		if ( argv[i][0] != '-' )
			continue;
		switch(argv[i][1]) {
			case 'f':
				filebase = argv[++i]; // input file name
				break;
			case 'o':
				offset = atof(argv[++i]) / GADGET_LENGTH_CONVERSION; // light cone vertex offset in Mpc/h
				break;
			case 'e':
				frac = atof(argv[++i]); // fraction of particles extracted
		}
	}

	if (filename == NULL)
	{
		cerr << " error! no input file base specified." << endl;
		return -1;
	}

	posbatch = (float *) malloc(3 * BATCHSIZE * sizeof(float));
	velbatch = (float *) malloc(3 * BATCHSIZE * sizeof(float));

	cout << "  X  Y  Z  VX  VY  VZ" << endl;

	do
	{
		sprintf(filename, "%s.%d", filebase, nfile);
		infile = fopen(filename, "rb");

		cerr << " reading " << filename << endl;

		if (infile == NULL)
		{
			cerr << " error! unable to open file " << filename << endl;
			return -1;
		}

		fread(&blocksize, sizeof(uint32_t), 1, infile);

		if (blocksize != sizeof(hdr))
		{
			cerr << " error! unknown file format " << filename << endl;
			fclose(infile);
			return -1;
		}

		if(fread(&hdr, sizeof(hdr), 1, infile) != 1)
		{
			cerr << " error! unable to read header block from " << filename << "!" << endl;
			fclose(infile);
			return -1;
		}

		/*cosmo.Omega_cdm = hdr.Omega0;
		cosmo.Omega_b = 0.;
		cosmo.Omega_m = cosmo.Omega_cdm;
		cosmo.Omega_Lambda = hdr.OmegaLambda;
		cosmo.Omega_g = 1. - (hdr.Omega0 + hdr.OmegaLambda);
		cosmo.Omega_ur = 0.;
		cosmo.Omega_rad = cosmo.Omega_g;
		cosmo.num_ncdm = 0;
		cosmo.h = hdr.HubbleParam;
		fourpiG = 1.5 * hdr.BoxSize * hdr.BoxSize * GADGET_LENGTH_CONVERSION * GADGET_LENGTH_CONVERSION / C_SPEED_OF_LIGHT / C_SPEED_OF_LIGHT;*/

		fread(&blocksize, sizeof(uint32_t), 1, infile);
		fread(&blocksize, sizeof(uint32_t), 1, infile);

		blockoffset = 3l * sizeof(float) * (long) hdr.npart[1] + 2l * sizeof(uint32_t);
		
		for (uint32_t n = 0; n < hdr.npart[1]; n += batch)
		{
			batch = (hdr.npart[1] - n >= BATCHSIZE) ? BATCHSIZE : (hdr.npart[1] - n);

			if (fread(posbatch, sizeof(float), 3l*batch, infile) != 3l*batch)
			{
				cerr << " error! unable to read position data from " << filename << endl;
				fclose(infile);
				return -1;
			}

			backtrack = ftell(infile);

			if (fseek(infile, blockoffset - 3l * batch * sizeof(float), SEEK_CUR))
			{
				cerr << " error! unable to fast forward to velocities block in " << filename << endl;
				fclose(infile);
				return -1;
			}

			if (fread(velbatch, sizeof(float), 3l*batch, infile) != 3l*batch)
			{
				cerr << " error! unable to read velocity data from " << filename << endl;
				fclose(infile);
				return -1;
			}

			if (fseek(infile, backtrack, SEEK_SET))
			{
				cerr << " error! unable to backtrack to positions block in " << filename << endl;
				fclose(infile);
				return -1;
			}

			for (long i = 0; i < batch; i++)
			{
				if (gsl_rng_uniform(rng) < frac)
				{
					//chi = sqrt((posbatch[3*i]-offset)*(posbatch[3*i]-offset) + (posbatch[3*i+1]-offset)*(posbatch[3*i+1]-offset) + (posbatch[3*i+2]-offset)*(posbatch[3*i+2]-offset)) / hdr.BoxSize;
					//a = 1.;//hdr.time;
					//for (int t = 0; t < 20; t++) rungekutta4bg(a, fourpiG, cosmo, -chi/20.);
					cout << " " << (posbatch[3*i]-offset)*GADGET_LENGTH_CONVERSION << "  " << (posbatch[3*i+1]-offset)*GADGET_LENGTH_CONVERSION << "  " << (posbatch[3*i+2]-offset)*GADGET_LENGTH_CONVERSION << "  " << /* GADGET_VELOCITY_CONVERSION * sqrt(a) */ velbatch[3*i] << "  " << /* GADGET_VELOCITY_CONVERSION * sqrt(a) */ velbatch[3*i+1] << "  " << /* GADGET_VELOCITY_CONVERSION * sqrt(a) */ velbatch[3*i+2] << endl;
				}
			}
		}

		fclose(infile);

		nfile++;
	}
	while (nfile < hdr.num_files);

	free(posbatch);
	free(velbatch);

	gsl_rng_free(rng);

	return 0;
}

