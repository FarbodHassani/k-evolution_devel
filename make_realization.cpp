#include <stdlib.h>
#include <stdio.h>
#include <cstdint>
#include "LATfield2.hpp"
#include "metadata.hpp"
#include "tools.hpp"
#include "background.hpp"
#include "Particles_gevolution.hpp"
#include "gevolution.hpp"
#include "ic_basic.hpp"

using namespace std;


int main(int argc, char **argv)
{
	char * filebase = NULL;
	int n = 0;
	int m = 0;
	int numpts = 0, i;
	icsettings ic;
	int box[3];
	double h = P_HUBBLE;
	double boxsize;
	gsl_spline * dspline = NULL;
	gsl_spline * tspline = NULL;
	char filename[1024];
	char * species;
	double * temp;

#ifndef H5_DEBUG
	H5Eset_auto2 (H5E_DEFAULT, NULL, NULL);
#endif
	
	for (i=1 ; i < argc ; i++ ){
		if ( argv[i][0] != '-' )
			continue;
		switch(argv[i][1]) {
			case 'f':
				filebase = argv[++i]; //file base for halo file
				break;
			case 'N':
				numpts = atoi(argv[++i]); //Ngrid
				break;
			case 'x':
				species = argv[++i]; //species name
				break;
			case 's':
				ic.seed = atoi(argv[++i]); //seed
				break;
			case 'b':
				boxsize = atof(argv[++i]); //boxsize
				break;
			case 'h':
				h = atof(argv[++i]); //h
				break;
			case 'n':
				n = atoi(argv[++i]); //size of the dim 1 of the processor grid
				break;
			case 'm':
				m = atoi(argv[++i]); //size of the dim 2 of the processor grid
		}
	}

	parallel.initialize(n,m);

	if (filebase == NULL)
	{
		COUT << " error: no file base for halo file specified! (option -f)" << endl;
		return -1;
	}

	if (numpts < n)
	{
		COUT << " error: Ngrid not set properly! (option -N)" << endl;
		return -1;
	}

	COUT << " Make realization tool" << endl;
	COUT << " Ngrid = " << numpts << "; species = " << species << "; seed = " << ic.seed << "; boxsize = " << boxsize << " Mpc/h; h = " << h << endl << endl;

	box[0] = numpts;
	box[1] = numpts;
	box[2] = numpts;

	ic.A_s = P_SPECTRAL_AMP;
	ic.n_s = P_SPECTRAL_INDEX;
	ic.k_pivot = P_PIVOT_SCALE;

	Lattice lat(3,box,1);
	Lattice latFT;
	latFT.initializeRealFFT(lat,0);

	Site x(lat);

	Field<Real> delta;
	Field<Cplx> scalarFT;

	delta.initialize(lat,1);
	scalarFT.initialize(latFT,1);
	PlanFFT<Cplx> plan_delta(&delta, &scalarFT);

	sprintf(filename, "%s_tk.dat", filebase);

	loadTransferFunctions(filename, dspline, tspline, species, boxsize, h);

	temp = (double *) malloc(dspline->size * sizeof(double));

	for (i = 0; i < dspline->size; i++)
		temp[i] = -dspline->y[i] * M_PI * sqrt(Pk_primordial(dspline->x[i] * h / boxsize, ic) / dspline->x[i]) / dspline->x[i];
	
	gsl_spline_free(dspline);
	dspline = gsl_spline_alloc(gsl_interp_cspline, tspline->size);
	gsl_spline_init(dspline, tspline->x, temp, tspline->size);
	gsl_spline_free(tspline);

	generateRealization(scalarFT, 0., dspline, (unsigned int) ic.seed, ICFLAG_KSPHERE);
	
	gsl_spline_free(dspline);
	free(temp);

	plan_delta.execute(FFT_BACKWARD);

	string h5filename(filebase);
	delta.saveHDF5(h5filename + "_" + species + ".h5");
}

