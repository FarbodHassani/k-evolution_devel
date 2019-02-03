#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <gsl/gsl_rng.h>
//#include <gsl/gsl_spline.h>

// g++ create_random_catalog.cpp -o create_random_catalog -std=c++11 -O3 -lgsl -lgslcblas

using namespace std;


int main(int argc, char **argv)
{
	FILE * infile = NULL;
	gsl_rng * rng = gsl_rng_alloc(gsl_rng_ranlxs1);
	long count = 0;
	unsigned long s = 0;
//	std::multiset<float> zvalues;
//	gsl_interp_accel * acc;
//	gsl_spline * spline;
//	double * zdata;
//	double * hdata;
	std::vector<float> zvalues;
	float z;
	float mincos = -1.;

	if (argc < 3) return -1;
	
	count = atol(argv[2]);
	
	if (argc > 3) mincos = atof(argv[3]);
	
	if (argc > 4) s = atol(argv[4]);
	
	gsl_rng_set(rng, s);

	infile = fopen(argv[1], "r");

	if (infile == NULL) return -1;
	
	while (true)
	{
		if (fscanf(infile, " %f %*f %*f \n", &z) != 1) break;
		//zvalues.insert(z);
		zvalues.push_back(z);
	}
	
	fclose(infile);

	if (zvalues.size() < 2) return -1;

	qsort((void *) zvalues.data(), zvalues.size(), sizeof(float), [](const void * a, const void * b)->int{ if (*((float *) a) < *((float *) b)) return -1; else if (*((float *) a) > *((float *) b)) return 1; else return 0; });

/*	zdata = (double *) malloc(zvalues.size()*sizeof(double));
	hdata = (double *) malloc(zvalues.size()*sizeof(double));

	s = 0;
	for (std::multiset<float>::iterator it = zvalues.begin(); it != zvalues.end(); it++)
	{
		zdata[s] = (*it);
		hdata[s] = (double) s / (double) (zvalues.size()-1);
		s++;
	}

	zvalues.clear();

	acc = gsl_interp_accel_alloc();
	spline = gsl_spline_alloc(gsl_interp_linear, s);
	gsl_spline_init(spline, hdata, zdata, s);

	free(hdata);
	free(zdata);*/
	
	for (long l = 0; l < count; l++)
	{
		//cout << " " << gsl_spline_eval(spline, gsl_rng_uniform_pos(rng), acc) << "  " << mincos + (1.-mincos) * gsl_rng_uniform_pos(rng) << "  " << 2. * M_PI * gsl_rng_uniform(rng) << endl;
		z = gsl_rng_uniform_pos(rng) * (zvalues.size()-1);
		s = (unsigned long) floor(z);
		z -= s;
		printf(" %f  %.10f  %.8f\n", (1.-z)*zvalues[s]+z*zvalues[s+1], mincos + (1.-mincos) * gsl_rng_uniform_pos(rng), 2. * M_PI * gsl_rng_uniform(rng));
	}

	gsl_rng_free(rng);
//	gsl_spline_free(spline);
//	gsl_interp_accel_free(acc);
	
	return 0;
}

