#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "LATfield2.hpp"

using namespace std;

using namespace LATfield2;

int main(int argc, char **argv)
{
	H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
	
    FILE * halofile;
    char * filename;
    double pos[3];
    double gradxi[3];
    double w[3];
    double tmp;
    int base[3];
    int count = 0;
    double disp = 0.;
    
    int n = 0, m = 0;

	for (int i = 1; i < argc; i++ ){
		if ( argv[i][0] != '-' )
			continue;
		switch(argv[i][1]) {
			case 'f':
				filename = argv[++i]; // halofilename
				break;
			case 'n':
				n = atoi(argv[++i]); //size of the dim 1 of the processor grid
				break;
			case 'm':
				m =  atoi(argv[++i]); //size of the dim 2 of the processor grid
		}
	}
	
	parallel.initialize(n,m);
	
	int dim = 3;
	int latSize[3] = {960,960,960};
	int halo = 1;
	Lattice lat(dim,latSize,halo);
	Field<Real> xi(lat,1);
	Site x(lat);
	
	xi.loadHDF5("output960/lc_9600box_displacement.h5");
	
	xi.updateHalo();
	
	x.first();
	
	cerr << " task " << parallel.rank() << " read complete, first value = " << xi(x) << endl;
	
	halofile = fopen(filename, "r");
	
	while (true)
	{
		if (fscanf(halofile, "%lf, %lf, %lf", pos, pos+1, pos+2) != 3) break;
		
		for (int i = 0; i < 3; i++)
		{
			w[i] = modf(480.+pos[i]/10., &tmp);
			base[i] = (int) tmp;
		}
		
		if (x.setCoord(base))
		{
			gradxi[0] = (1.-w[1]) * (1.-w[2]) * (xi(x+0) - xi(x));
			gradxi[1] = (1.-w[0]) * (1.-w[2]) * (xi(x+1) - xi(x));
			gradxi[2] = (1.-w[0]) * (1.-w[1]) * (xi(x+2) - xi(x));
			gradxi[0] += w[1] * (1.-w[2]) * (xi(x+1+0) - xi(x+1));
			gradxi[1] += w[0] * (1.-w[2]) * (xi(x+1+0) - xi(x+0));
			gradxi[2] += w[0] * (1.-w[1]) * (xi(x+2+0) - xi(x+0));
			gradxi[0] += (1.-w[1]) * w[2] * (xi(x+2+0) - xi(x+2));
			gradxi[1] += (1.-w[0]) * w[2] * (xi(x+2+1) - xi(x+2));
			gradxi[2] += (1.-w[0]) * w[1] * (xi(x+2+1) - xi(x+1));
			gradxi[0] += w[1] * w[2] * (xi(x+2+1+0) - xi(x+2+1));
			gradxi[1] += w[0] * w[2] * (xi(x+2+1+0) - xi(x+2+0));
			gradxi[2] += w[0] * w[1] * (xi(x+2+1+0) - xi(x+1+0));
			
			for (i = 0; i < 3; i++) pos[i] += gradxi[i] * 960. * 9600.;
			
			if (disp < gradxi[0]*gradxi[0] + gradxi[1]*gradxi[1] + gradxi[2]*gradxi[2])
				disp = gradxi[0]*gradxi[0] + gradxi[1]*gradxi[1] + gradxi[2]*gradxi[2];
				
			printf(" %d\t%.5lf\t%.5lf\t%.5lf\n", count, pos[0], pos[1], pos[2]);
		}
		
		count++;
	}
	
	fclose(halofile);
	
	cerr << " task " << parallel.rank() << " max displacement = " << sqrt(disp)*960.*960. << " lattice units." << endl;
}
	
