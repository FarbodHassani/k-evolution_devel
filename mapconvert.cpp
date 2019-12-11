#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "chealpix.h"

// g++ mapconvert.cpp -o mapconvert -std=c++11 -O3 -I/astro/adamek/local/include -L/astro/adamek/local/lib -lchealpix -lcfitsio

using namespace std;

int main(int argc, char **argv)
{
	FILE * infile;
	float * pix;
	char coordsys = 'G';
	uint32_t * data;
	int c = 0;
	float renorm;

	if (argc < 3) return 0;
	
	if (argc > 3) c = atoi(argv[3]);
	
	renorm = 3. * 0.3125 * 0.3125 * 0.3125 * 12 * 1024 * 1024 / (4000. * M_PI * ((c+1)*(c+1)*(c+1) - c*c*c));
	
	infile = fopen(argv[1], "rb");
	
	if (infile == NULL) return 0;
	
	data = (uint32_t *) malloc(12l*1024*1024*sizeof(uint32_t));
	pix = (float *) malloc(12l*1024*1024*sizeof(float));
	
	fread(data, sizeof(uint32_t), 12l*1024*1024, infile);
	
	fclose(infile);
	
	for (int i = 0; i < 12*1024*1024; i++)
		pix[i] = ((float) data[i]) * renorm;
		
	free(data);
	
	write_healpix_map(pix, 1024l, argv[2], 0, &coordsys);
	
	free(pix);
	
	return 0;
}

