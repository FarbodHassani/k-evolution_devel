#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "chealpix.h"

using namespace std;

int main(int argc, char **argv)
{
	FILE * infile = NULL;
	float * map = NULL;
	int Nside = 0;
	char coordsys = 'G';
	
	if (argc < 4) return -1;
	
	Nside = atoi(argv[3]);
	
	cout << " Nside = " << Nside << endl;
	
	map = (float *) malloc(sizeof(float) * 12 * Nside * Nside);
	
	infile = fopen(argv[1], "rb");
	
	if (infile == NULL)
	{
		cout << " error opening file " << argv[1] << " for reading!" << endl;
		return -1;
	}
	
	if (fread((void *) map, sizeof(float), 12 * Nside * Nside, infile) != 12 * Nside * Nside)
	{
		cout << " error reading data from " << argv[1] << "!" << endl;
		return -1;
	}
	
	fclose(infile);
	
	cout << " writing data to " << argv[2] << endl;
	
	write_healpix_map(map, Nside, argv[2], 0, &coordsys);
	
	free(map);
	
	return 0;
}

