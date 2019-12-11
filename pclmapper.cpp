#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include <chealpix.h>
#include <cstdint>

using namespace std;

struct io_header_1
{
  uint32_t npart[6];
  double mass[6];
  double time;
  double redshift;
  int32_t flag_sfr;
  int32_t flag_feedback;
  uint32_t npartTotal[6];
  int32_t flag_cooling;
  uint32_t num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  int32_t flag_age;
  int32_t flag_metals;
  uint32_t npartTotalHW[6];
  char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8 - 2 * 4 - 6 * 4];       /* fills to 256 Bytes */
} header1;

#define BATCHSIZE 134217728

int main (int argc, char **argv)
{
	long nside = 0;
	long npix = 0;
	long ipix;
	double vec[3];
	int  nmaps = 1;
	int  nfile = 1;
	char * infilebase = NULL;
	char * outfilebase = NULL;
	char filename[1024];
	uint32_t ** data;
	float * pcldata;
	float r;
	long npart = 0;
	long batch;
	FILE * pclfile;
	FILE * mapfile;
	uint32_t blocksize1, blocksize2;

	if (argc < 2)
	{
		cout << "use option --help for help." << endl;
		cout << "don't know what to do, exiting..." << endl;
		return 0;
	}

	for (int i = 1; i < argc; i++)
	{
		if (!strcmp(argv[i], "--help"))
		{
			cout << "usage: parameters are passed as --<option>=<setting>" << endl << endl;
			cout << "options:" << endl << endl;
			cout << " --nside          Nside of HEALPix maps" << endl;
			cout << " --input          input file basename (Gadget2 files)" << endl;
			cout << " --output         output file basename (binary files containing HEALPix map)" << endl;
			cout << " --nmaps          number of HEALPix maps" << endl;
			cout << " --nfile          number of Gadget2 files" << endl << endl;
			return 0;
		}
		else if (!strncmp(argv[i], "--nside=", 8))
		{
			nside = atol(argv[i]+8);
			cout << "Nside = " << nside << endl;
			if (nside < 2)
			{
				cout << "Nside must be larger than 1!" << endl;
				return 0;
			}
		}
		else if (!strncmp(argv[i], "--output=", 9))
		{
			if (argv[i][9] == '\0')
			{
				cout << "output file basename empty!" << endl;
				return 0;
			}
			else
			{
				outfilebase = argv[i]+9;
				cout << "output file basename = " << outfilebase << endl;
			}
		}
		else if (!strncmp(argv[i], "--input=", 8))
		{
			if (argv[i][8] == '\0')
			{
				cout << "input file basename empty!" << endl;
				return 0;
			}
			else
			{
				infilebase = argv[i]+8;
				cout << "input file basename = " << infilebase << endl;
			}
		}
		else if (!strncmp(argv[i], "--nmaps=", 8))
		{
			nmaps = atol(argv[i]+8);
			cout << "number of maps = " << nmaps << endl;
			if (nmaps < 1)
			{
				cout << "number of maps must be larger than 0!" << endl;
				return 0;
			}
		}
		else if (!strncmp(argv[i], "--nfile=", 8))
		{
			nfile = atol(argv[i]+8);
			cout << "number of files = " << nfile << endl;
			if (nfile < 1)
			{
				cout << "number of files must be larger than 0!" << endl;
				return 0;
			}
		}
		else
		{
			cout << "unrecognized parameter: " << argv[i] << endl << "exiting..." << endl;
			return 0;
		}
	}

	if (nside < 2)
	{
		cout << "Nside not specified!" << endl;
		return 0;
	}
	else if (infilebase == NULL)
	{
		cout << "input file name not specified!" << endl;
		return 0;
	}
	else if (outfilebase == NULL)
	{
		cout << "output file name not specified!" << endl;
		return 0;
	}

	npix = nside2npix(nside);

	data = (uint32_t **) malloc(nmaps * sizeof(uint32_t*));

	for (int i = 0; i < nmaps; i++)
	{
		data[i] = (uint32_t *) malloc(npix * sizeof(uint32_t));
		if (data[i] == NULL)
		{
			cout << " memory error! could not allocate memory for map #" << i << endl;
			return 0;
		}
		for (int j = 0; j < npix; j++) data[i][j] = 0;
	}
	
	pcldata = (float *) malloc(BATCHSIZE * 3l * sizeof(float));
	
	if (pcldata == NULL)
	{
		cout << " memory error! could not allocate memory for reading particles" << endl;
		return 0;
	}

	for (int c = 0; c < nfile; c++)
	{
		sprintf(filename, "%s.%d", infilebase, c);
		cout << " reading " << filename << endl;
		pclfile = fopen(filename, "rb");
		
		if (pclfile == NULL)
  		{
    			cout << "file not found!" << endl;
    			return 0;
  		}
  		
  		fread(&blocksize1, sizeof(uint32_t), 1, pclfile);
    		fread(&header1, sizeof(header1), 1, pclfile);
    		fread(&blocksize2, sizeof(uint32_t), 1, pclfile);

		cout << " file contains " << header1.npart[1] << " particles" << endl;
		
		fread(&blocksize1, sizeof(uint32_t), 1, pclfile);
		
		while (header1.npart[1] > 0)
		{
			batch = (header1.npart[1] > BATCHSIZE) ? BATCHSIZE : header1.npart[1];
			
			fread(pcldata, sizeof(float), 3l * batch, pclfile);
			
#pragma omp parallel for private(ipix,r,vec)
			for (long p = 0; p < batch; p++)
			{
				r = sqrt(pcldata[3*p]*pcldata[3*p]+pcldata[3*p+1]*pcldata[3*p+1]+pcldata[3*p+2]*pcldata[3*p+2]);
				if ((r < 256920.91 && c < 8) || (r >= 256920.91 && c >= 8 && r < 4000000.0))
				{
					vec[0] = pcldata[3*p];
					vec[1] = pcldata[3*p+1];
					vec[2] = pcldata[3*p+2];
					
					vec2pix_ring(nside, vec, &ipix);
#pragma omp critical
{
					data[(int) floor((double) nmaps * r / 4000000.0)][ipix]++;
}
				}
			}
			
			header1.npart[1] -= batch;
			npart += batch;
		}
		
		fclose(pclfile);
	}
	
	free(pcldata);

	cout << " " << npart << " particles read. writing maps..." << endl;

	for (int i = 0; i < nmaps; i++)
	{
		sprintf(filename, "%s%03d.map", outfilebase, i);
		
		mapfile = fopen(filename, "wb");
		
		if (mapfile == NULL)
  		{
    			cout << "unable to open file for output!" << endl;
    			return 0;
  		}
		
		fwrite(data[i], sizeof(uint32_t), npix, mapfile);

		fclose(mapfile);

		free(data[i]);
	}

	cout << " " << nmaps << " maps written. normal completion." << endl;

	free(data);
}

