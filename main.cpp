//////////////////////////
// Copyright (c) 2015-2016 Julian Adamek (Université de Genève)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESSED OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//////////////////////////

//////////////////////////
// main.cpp
//////////////////////////
//
// main control sequence of Geneva N-body code with evolution of metric perturbations (gevolution)
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris)
//
// Last modified: February 2017
//
//////////////////////////
#include <stdlib.h>
#ifdef BACKREACTION_TEST
#include <iomanip>
#include<fstream>
#include<sstream>
#endif
#ifdef HAVE_CLASS
#include "class.h"
#undef MAX			// due to macro collision this has to be done BEFORE including LATfield2 headers!
#undef MIN
#endif
#include "LATfield2.hpp"
#include "metadata.hpp"
#include "class_tools.hpp"
#include "background.hpp"
#include "Particles_gevolution.hpp"
#include "gevolution.hpp"
#include "ic_basic.hpp"
#include "ic_read.hpp"
#ifdef ICGEN_PREVOLUTION
#include "ic_prevolution.hpp"
#endif
#ifdef ICGEN_FALCONIC
#include "fcn/togevolution.hpp"
#endif
#include "radiation.hpp"
#include "parser.hpp"
#include "tools.hpp"
#include "output.hpp"
#include "hibernation.hpp"

using namespace std;

using namespace LATfield2;

int main(int argc, char **argv)
{

#ifdef BENCHMARK
	//benchmarking variables

	double ref_time, ref2_time, cycle_start_time;

	double initialization_time;
	double run_time;
  double kessence_update_time=0; // How much time is put on updating the kessence field
	double cycle_time=0;
	double projection_time = 0;
	double snapshot_output_time = 0;
	double spectra_output_time = 0;
	double gravity_solver_time = 0;
	double fft_time = 0;
	int fft_count = 0;
	double update_q_time = 0;
	int update_q_count = 0;
	double moveParts_time = 0;
	int  moveParts_count =0;
	//kessence
	double a_kess;


#endif  //BENCHMARK

	int n = 0, m = 0;
	int io_size = 0;
	int io_group_size = 0;

	int i, j, cycle = 0, snapcount = 0, pkcount = 0, restartcount = 0, usedparams, numparam = 0, numsteps, numspecies;
	int numsteps_ncdm[MAX_PCL_SPECIES-2];
	long numpts3d;
	int box[3];
	double dtau, dtau_old, dx, tau, a, fourpiG, tau_Lambda, tmp, start_time;
	double maxvel[MAX_PCL_SPECIES];
	FILE * outfile;
	char filename[2*PARAM_MAX_LENGTH+24];
	string h5filename;
	char * settingsfile = NULL;
	parameter * params = NULL;
	metadata sim;
	cosmology cosmo;
	icsettings ic;
	gadget2_header hdr;
	Real T00hom;

#ifndef H5_DEBUG
	H5Eset_auto2 (H5E_DEFAULT, NULL, NULL);
#endif

	for (i=1 ; i < argc ; i++ ){
		if ( argv[i][0] != '-' )
			continue;
		switch(argv[i][1]) {
			case 's':
				settingsfile = argv[++i]; //settings file name
				break;
			case 'n':
				n = atoi(argv[++i]); //size of the dim 1 of the processor grid
				break;
			case 'm':
				m =  atoi(argv[++i]); //size of the dim 2 of the processor grid
				break;
			case 'i':
#ifndef EXTERNAL_IO
				cout << "EXTERNAL_IO needs to be set at compilation to use the I/O server"<<endl;
				exit(-1000);
#endif
				io_size =  atoi(argv[++i]);
				break;
			case 'g':
#ifndef EXTERNAL_IO
				cout << "EXTERNAL_IO needs to be set at compilation to use the I/O server"<<endl;
				exit(-1000);
#endif
				io_group_size = atoi(argv[++i]);
		}
	}

#ifndef EXTERNAL_IO
	parallel.initialize(n,m);
#else
	parallel.initialize(n,m,io_size,io_group_size);
	if(parallel.isIO()) ioserver.start();
	else
	{
#endif

	COUT << COLORTEXT_WHITE << endl;
	COUT << "| /  "<<endl;
	COUT << "|/    _      _   _  _ , _" << endl;
	COUT << "|\\  (-' \\/ (_) (_ (_| (  ( (_) /\\/	version 1.0    running on " << n*m << " cores." << endl;
  COUT << "| \\" << endl << COLORTEXT_RESET << endl;

	if (settingsfile == NULL)
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no settings file specified!" << endl;
		parallel.abortForce();
	}

	COUT << " initializing..." << endl;

	start_time = MPI_Wtime();

	numparam = loadParameterFile(settingsfile, params);

	usedparams = parseMetadata(params, numparam, sim, cosmo, ic);

	COUT << " parsing of settings file completed. " << numparam << " parameters found, " << usedparams << " were used." << endl;

	sprintf(filename, "%s%s_settings_used.ini", sim.output_path, sim.basename_generic);
	saveParameterFile(filename, params, numparam);

	free(params);

#ifdef HAVE_CLASS
	background class_background;
  	perturbs class_perturbs;
  	spectra class_spectra;
#endif

	h5filename.reserve(2*PARAM_MAX_LENGTH);
	h5filename.assign(sim.output_path);
	h5filename += sim.basename_snapshot;

	box[0] = sim.numpts;
	box[1] = sim.numpts;
	box[2] = sim.numpts;

	Lattice lat(3,box,1);
	Lattice latFT;
	latFT.initializeRealFFT(lat,0);

	Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> pcls_cdm;
	Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> pcls_b;
	Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> pcls_ncdm[MAX_PCL_SPECIES-2];
	Field<Real> * update_cdm_fields[3];
	Field<Real> * update_b_fields[3];
	Field<Real> * update_ncdm_fields[3];
	double f_params[5];

	Field<Real> phi;
	Field<Real> source;
	Field<Real> chi;
	Field<Real> chi_old;
	Field<Real> Sij;
	Field<Real> Bi;
	Field<Cplx> scalarFT;
	Field<Cplx> SijFT;
	Field<Cplx> BiFT;
	source.initialize(lat,1);
	phi.initialize(lat,1);

	//kessence
	Field<Real> phi_old;
  //phi at two step before to compute phi'(n+1/2)
	Field<Real> phi_prime;
	Field<Real> pi_k;
	// Field<Real> zeta_integer;
  Field<Real> zeta_half;
	Field<Real> T00_Kess;
	Field<Real> T0i_Kess;
	Field<Real> Tij_Kess;
	Field<Cplx> scalarFT_phi_old;
	Field<Cplx> phi_prime_scalarFT;
	Field<Cplx> scalarFT_chi_old;
	Field<Cplx> scalarFT_pi;
	// Field<Cplx> scalarFT_zeta_integer;
  Field<Cplx> scalarFT_zeta_half;
	Field<Cplx> T00_KessFT;
	Field<Cplx> T0i_KessFT;
	Field<Cplx> Tij_KessFT;
	chi.initialize(lat,1);
	scalarFT.initialize(latFT,1);
	PlanFFT<Cplx> plan_source(&source, &scalarFT);
	PlanFFT<Cplx> plan_phi(&phi, &scalarFT);
	PlanFFT<Cplx> plan_chi(&chi, &scalarFT);
	Sij.initialize(lat,3,3,symmetric);
	SijFT.initialize(latFT,3,3,symmetric);
	PlanFFT<Cplx> plan_Sij(&Sij, &SijFT);
	Bi.initialize(lat,3);
	BiFT.initialize(latFT,3);
	PlanFFT<Cplx> plan_Bi(&Bi, &BiFT);
#ifdef CHECK_B
	Field<Real> Bi_check;
	Field<Cplx> BiFT_check;
	Bi_check.initialize(lat,3);
	BiFT_check.initialize(latFT,3);
	PlanFFT<Cplx> plan_Bi_check(&Bi_check, &BiFT_check);
#endif
	//Kessence end

	//Kessence part initializing
	//Phi_old
	phi_old.initialize(lat,1);
	scalarFT_phi_old.initialize(latFT,1);
	PlanFFT<Cplx> plan_phi_old(&phi_old, &scalarFT_phi_old);
	//Phi'
	phi_prime.initialize(lat,1);
	phi_prime_scalarFT.initialize(latFT,1);
	PlanFFT<Cplx> phi_prime_plan(&phi_prime, &phi_prime_scalarFT);
	//pi_k kessence
	pi_k.initialize(lat,1);
	scalarFT_pi.initialize(latFT,1);
	PlanFFT<Cplx> plan_pi_k(&pi_k, &scalarFT_pi);
	//zeta_integer_k kessence
	// zeta_half.initialize(lat,1);
	// scalarFT_zeta_half.initialize(latFT,1);
	// PlanFFT<Cplx> plan_zeta_half(&zeta_half, &scalarFT_zeta_half);
  //zeta_half_k kessence
  zeta_half.initialize(lat,1);
  scalarFT_zeta_half.initialize(latFT,1);
  PlanFFT<Cplx> plan_zeta_half(&zeta_half, &scalarFT_zeta_half);
	//chi_old initialize
	chi_old.initialize(lat,1);
	scalarFT_chi_old.initialize(latFT,1);
	PlanFFT<Cplx> plan_chi_old(&chi_old, &scalarFT_chi_old);
	//Stress tensor initializing
	T00_Kess.initialize(lat,1);
	T00_KessFT.initialize(latFT,1);
	PlanFFT<Cplx> plan_T00_Kess(&T00_Kess, &T00_KessFT);
	// T00_Kess.alloc();  // It seems we don't need it!
	T0i_Kess.initialize(lat,3);
	T0i_KessFT.initialize(latFT,3);
	PlanFFT<Cplx> plan_T0i_Kess(&T0i_Kess, &T0i_KessFT);
	// T0i_Kess.alloc();
	Tij_Kess.initialize(lat,3,3,symmetric);
	Tij_KessFT.initialize(latFT,3,3,symmetric);
	PlanFFT<Cplx> plan_Tij_Kess(&Tij_Kess, &Tij_KessFT);
	// Tij_Kess.alloc();
	// kessence end


	update_cdm_fields[0] = &phi;
	update_cdm_fields[1] = &chi;
	update_cdm_fields[2] = &Bi;

	update_b_fields[0] = &phi;
	update_b_fields[1] = &chi;
	update_b_fields[2] = &Bi;

	update_ncdm_fields[0] = &phi;
	update_ncdm_fields[1] = &chi;
	update_ncdm_fields[2] = &Bi;

	Site x(lat);
	rKSite kFT(latFT);

	dx = 1.0 / (double) sim.numpts;
	numpts3d = (long) sim.numpts * (long) sim.numpts * (long) sim.numpts;

	for (i = 0; i < 3; i++) // particles may never move farther than to the adjacent domain
	{
		if (lat.sizeLocal(i)-1 < sim.movelimit)
			sim.movelimit = lat.sizeLocal(i)-1;
	}
	parallel.min(sim.movelimit);

	fourpiG = 1.5 * sim.boxsize * sim.boxsize / C_SPEED_OF_LIGHT / C_SPEED_OF_LIGHT; // Just a definition to make Friedmann equation simplofied! and working with normak numbers
	a = 1. / (1. + sim.z_in);
	tau = particleHorizon(a, fourpiG, cosmo);
	tau_Lambda = -1.0;

	if (sim.Cf * dx < sim.steplimit / Hconf(a, fourpiG, cosmo))
		// dtau = sim.Cf * dx / sim.nKe_numsteps;
    dtau = sim.Cf * dx;

	else
		// dtau = sim.steplimit / Hconf(a, fourpiG, cosmo) / sim.nKe_numsteps;
    dtau = sim.steplimit / Hconf(a, fourpiG, cosmo);

	dtau_old = 0.;

	if (ic.generator == ICGEN_BASIC)
		generateIC_basic(sim, ic, cosmo, fourpiG, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &pi_k, &zeta_half, &chi, &Bi, &source, &Sij, &scalarFT, &scalarFT_pi, &scalarFT_zeta_half, &BiFT, &SijFT, &plan_phi, &plan_pi_k, &plan_zeta_half, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
	// generates ICs on the fly
	else if (ic.generator == ICGEN_READ_FROM_DISK)
		readIC(sim, ic, cosmo, fourpiG, a, tau, dtau, dtau_old, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, cycle, snapcount, pkcount, restartcount);
#ifdef ICGEN_PREVOLUTION
	else if (ic.generator == ICGEN_PREVOLUTION)
		generateIC_prevolution(sim, ic, cosmo, fourpiG, a, tau, dtau, dtau_old, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
#endif
#ifdef ICGEN_FALCONIC
	else if (ic.generator == ICGEN_FALCONIC)
		maxvel[0] = generateIC_FalconIC(sim, ic, cosmo, fourpiG, dtau, &pcls_cdm, pcls_ncdm, maxvel+1, &phi, &source, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_source, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
#endif
	else
	{
		COUT << " error: IC generator not implemented!" << endl;
		parallel.abortForce();
	}

	if (sim.baryon_flag > 1)
	{
		COUT << " error: baryon_flag > 1 after IC generation, something went wrong in IC generator!" << endl;
		parallel.abortForce();
	}

	numspecies = 1 + sim.baryon_flag + cosmo.num_ncdm;
	parallel.max<double>(maxvel, numspecies);

	if (sim.gr_flag > 0)
	{
		for (i = 0; i < numspecies; i++)
			maxvel[i] /= sqrt(maxvel[i] * maxvel[i] + 1.0);
	}

#ifdef CHECK_B
	if (sim.vector_flag == VECTOR_ELLIPTIC)
	{
		for (kFT.first(); kFT.test(); kFT.next())
		{
			BiFT_check(kFT, 0) = BiFT(kFT, 0);
			BiFT_check(kFT, 1) = BiFT(kFT, 1);
			BiFT_check(kFT, 2) = BiFT(kFT, 2);
		}
	}
#endif

	for (i = 0; i < 6; i++)
	{
		hdr.npart[i] = 0;
		hdr.npartTotal[i] = 0;
		hdr.mass[i] = 0.;
	}
	hdr.num_files = 1;
	hdr.Omega0 = cosmo.Omega_m;
	hdr.OmegaLambda = cosmo.Omega_Lambda;
	hdr.HubbleParam = cosmo.h;
	hdr.BoxSize = sim.boxsize / GADGET_LENGTH_CONVERSION;
	hdr.flag_sfr = 0;
	hdr.flag_cooling = 0;
	hdr.flag_feedback = 0;
	for (i = 0; i < 256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8; i++)
		hdr.fill[i] = 0;

#ifdef BENCHMARK
	initialization_time = MPI_Wtime() - start_time;
	parallel.sum(initialization_time);
	COUT << COLORTEXT_GREEN << " initialization complete." << COLORTEXT_RESET << " BENCHMARK: " << hourMinSec(initialization_time) << endl << endl;
#else
	COUT << COLORTEXT_GREEN << " initialization complete." << COLORTEXT_RESET << endl << endl;
#endif

#ifdef HAVE_CLASS
	if (sim.radiation_flag > 0)
	{
		initializeCLASSstructures(sim, ic, cosmo, class_background, class_perturbs, class_spectra);
		if (sim.gr_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.) && (ic.generator == ICGEN_BASIC || (ic.generator == ICGEN_READ_FROM_DISK && cycle == 0)))
		{
			prepareFTchiLinear(class_background, class_perturbs, class_spectra, scalarFT, sim, ic, cosmo, fourpiG, a);
			plan_source.execute(FFT_BACKWARD);
			for (x.first(); x.test(); x.next())
				chi(x) += source(x);
			chi.updateHalo();
		}
	}
#endif


//INITIAL CONDITION If U wanna set it yourself
// for (x.first(); x.test(); x.next())
//   {
//     zeta_half(x)=1.e-3;
//     pi_k(x)=1.e-3;
//   }
// for (kFT.first(); kFT.test(); kFT.next())
// {
//   phi(k)=0;
// }

#ifdef BACKREACTION_TEST
//   //****************************
//   //****SAVE DATA To test Backreaction
//   //****************************
  FILE* Result_avg;
  FILE* Result_real;
  FILE* Result_fourier;
  FILE* Result_max;


  char filename_avg[60];
  char filename_real[60];
  char filename_fourier[60];
  char filename_max[60];


  snprintf(filename_avg, sizeof(filename_avg),"./output/Result_avg.txt");
  snprintf(filename_real, sizeof(filename_real),"./output/Result_real.txt");
  snprintf(filename_fourier, sizeof(filename_fourier),"./output/Result_fourier.txt");
  snprintf(filename_max, sizeof(filename_max),"./output/Results_max.txt");

  // ofstream out(filename_avg,ios::out);
  ofstream out_avg(filename_avg,ios::out);
  ofstream out_real(filename_real,ios::out);
  ofstream out_fourier(filename_fourier,ios::out);
  ofstream out_max(filename_max,ios::out);


  Result_avg=fopen(filename_avg,"w");
  Result_real=fopen(filename_real,"w");
  Result_fourier=fopen(filename_fourier,"w");
  Result_max=fopen(filename_max,"w");


  out_avg<<"### The result of the verage over time \n### d tau = "<< dtau<<endl;
  out_avg<<"### number of kessence update = "<<  sim.nKe_numsteps <<endl;
  out_avg<<"### initial time = "<< tau <<endl;
  out_avg<<"### 1- tau\t2- average(H pi_k)\t3- average (zeta)\t 4- average (phi)\t5-z(redshift)   " <<endl;


  out_max<<"### The result of the maximum over time \n### d tau = "<< dtau<<endl;
  out_max<<"### number of kessence update = "<<  sim.nKe_numsteps <<endl;
  out_max<<"### initial time = "<< tau <<endl;
  out_max<<"### 1- tau\t2- max(H pi_k)\t3- max (zeta)\t 4- max (phi)   " <<endl;


  out_real<<"### The result of the verage over time \n### d tau = "<< dtau<<endl;
  out_real<<"### number of kessence update = "<<  sim.nKe_numsteps <<endl;
  out_real<<"### initial time = "<< tau <<endl;
  out_real<<"### 1- tau\t2- pi_k(x)\t3-zeta(x)\t 4-x" <<endl;


  out_fourier<<"### The result of the verage over time \n### d tau = "<< dtau<<endl;
  out_fourier<<"### number of kessence update = "<<  sim.nKe_numsteps <<endl;
  out_fourier<<"### initial time = "<< tau <<endl;
  out_fourier<<"### 1- tau\t 2- pi_k(k)\t\t3-zeta(k)\t\t4-|k|\t\t 5-vec{k} \t 6-|k|^2"<<endl;

//defining the average
double avg_pi = 0.;
double avg_zeta = 0.;
double avg_phi = 0.;

double max_pi = 0.;
double max_zeta = 0.;
double max_phi = 0.;

int norm_kFT_squared = 0.;
#endif

	//******************************************************************
	//Write spectra check!
	// Kessence projection Tmunu Test IC
	//******************************************************************
	//  	if (sim.vector_flag == VECTOR_ELLIPTIC)
	// 		{
	// 			projection_Tmunu_kessence( T00_Kess,T0i_Kess,Tij_Kess, dx, a, phi, phi_old, chi, pi_k, zeta_integer_k, cosmo.Omega_kessence, cosmo.w_kessence, cosmo.cs2_kessence, Hconf(a, fourpiG, cosmo), fourpiG, 1 );
	// 		}
	//  	else
	// 		{
	// 			projection_Tmunu_kessence( T00_Kess,T0i_Kess,Tij_Kess, dx, a, phi, phi_old, chi, pi_k, zeta_integer_k, cosmo.Omega_kessence, cosmo.w_kessence, cosmo.cs2_kessence, Hconf(a, fourpiG, cosmo), fourpiG, 0 );
	// 		}
	//
// writeSpectra(sim, cosmo, fourpiG, a, pkcount, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &pi_k, &zeta_half, &chi, &Bi, &T00_Kess, &T0i_Kess, &Tij_Kess, &source, &Sij, &scalarFT ,&scalarFT_pi, &scalarFT_zeta_half, &BiFT, &T00_KessFT, &T0i_KessFT, &Tij_KessFT, &SijFT, &plan_phi, &plan_pi_k, &plan_zeta_half, &plan_chi, &plan_Bi, &plan_T00_Kess, &plan_T0i_Kess, &plan_Tij_Kess, &plan_source, &plan_Sij);

// writeSpectra_phi_prime(sim, cosmo, fourpiG, a, pkcount, &phi_prime, &phi_prime_scalarFT, &phi_prime_plan);

// writeSpectra(sim, cosmo, fourpiG, a, pkcount, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &pi_k, &zeta_half, &chi, &Bi, &T00_Kess, &T0i_Kess, &Tij_Kess, &source, &Sij, &scalarFT ,&scalarFT_pi, &scalarFT_zeta_half, &BiFT, &T00_KessFT, &T0i_KessFT, &Tij_KessFT, &SijFT, &plan_phi, &plan_pi_k, &plan_zeta_half, &plan_chi, &plan_Bi, &plan_T00_Kess, &plan_T0i_Kess, &plan_Tij_Kess, &plan_source, &plan_Sij);


	while (true)    // main loop
	{
    //Kessence
  	for (x.first(); x.test(); x.next())
  		{
        // cout<<"tau: "<<tau<<" z: "<<1./(a) -1.<<endl;
        //\Phi(n-1) = \Phi(old) and \Phi(n) which will be updated in this loops
        // Just note that in the first 2-3 steps it does not work since we
  			phi_old(x) =phi(x);
  			chi_old(x) =chi(x);
         // if(x.coord(0)==32 && x.coord(1)==12 && x.coord(2)==32) cout<<"zeta_half: "<<zeta_half(x)<<endl;
  		}
#ifdef BACKREACTION_TEST
      //****************************
      //****PRINTING AVERAGE OVER TIME
      //****************************
      // check_field(  zeta_half, 1. , " H pi_k", numpts3d);
      avg_pi =average(  pi_k, Hconf(a, fourpiG, cosmo), numpts3d ) ;
      avg_zeta =average(  zeta_half,1., numpts3d ) ;
      avg_phi =average(  phi , 1., numpts3d ) ;

      max_pi =maximum(  pi_k, Hconf(a, fourpiG, cosmo), numpts3d ) ;
      max_zeta =maximum(  zeta_half,1., numpts3d ) ;
      max_phi =maximum(  phi , 1., numpts3d ) ;

      COUT << scientific << setprecision(8);
      // if(parallel.isRoot())
      // {
        // fprintf(Result_avg,"\n %20.20e %20.20e ", tau, avg ) ;
      out_avg<<setw(9) << tau <<"\t"<< setw(9) << avg_pi<<"\t"<< setw(9) << avg_zeta<<"\t"<< setw(9) << avg_phi<<"\t"<< setw(9) << 1./a -1.<<endl;

        out_max<<setw(9) << tau <<"\t"<< setw(9) << max_pi<<"\t"<< setw(9) << max_zeta<<"\t"<< setw(9) << max_phi<<endl;

      // }
      //****************************
      //****PRINTING REAL SPACE INFO
      //****************************
      for (x.first(); x.test(); x.next())
    	{
          //NL_test, Printing out average
        if(x.coord(0)==32 && x.coord(1)==32 && x.coord(2)==32)
        {
          // if(parallel.isRoot())
          // {
          out_real<<setw(9) << tau <<"\t"<< setw(9) <<pi_k (x)<<"\t"<< setw(9)<<zeta_half (x)<<"\t"<<x<<endl;
          // }
        }
    	}
      //****************************
      //FOURIER PRINTING
      //****************************
      for(kFT.first();kFT.test();kFT.next())
      {
        norm_kFT_squared= kFT.coord(0)*kFT.coord(0) + kFT.coord(1) * kFT.coord(1) + kFT.coord(2) * kFT.coord(2);
        if(norm_kFT_squared == 1)
        {
          out_fourier<<setw(9) << tau <<"\t"<< setw(9) << scalarFT_pi(kFT)<<"\t"<< setw(9)<<scalarFT_zeta_half (kFT)<<"\t"<<kFT<<"\t"<<norm_kFT_squared<<endl;
        }
      }
      //**********************
      //END ADDED************
      //**********************
#endif

#ifdef BENCHMARK
		cycle_start_time = MPI_Wtime();
#endif
		// construct stress-energy tensor
		projection_init(&source);
#ifdef HAVE_CLASS
		if (sim.radiation_flag > 0)
			projection_T00_project(class_background, class_perturbs, class_spectra, source, scalarFT, &plan_source, sim, ic, cosmo, fourpiG, a);
#endif
		if (sim.gr_flag > 0)
		{
			projection_T00_project(&pcls_cdm, &source, a, &phi);
			if (sim.baryon_flag)
				projection_T00_project(&pcls_b, &source, a, &phi);
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (a >= 1. / (sim.z_switch_deltancdm[i] + 1.))
					projection_T00_project(pcls_ncdm+i, &source, a, &phi);
				else if (sim.radiation_flag == 0)
				{
					tmp = bg_ncdm(a, cosmo, i);
					for(x.first(); x.test(); x.next())
						source(x) += tmp;
				}
			}
		}
		else
		{
			scalarProjectionCIC_project(&pcls_cdm, &source);
			if (sim.baryon_flag)
				scalarProjectionCIC_project(&pcls_b, &source);
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (a >= 1. / (sim.z_switch_deltancdm[i] + 1.))
					scalarProjectionCIC_project(pcls_ncdm+i, &source);
			}
		}
		projection_T00_comm(&source);

		if (sim.vector_flag == VECTOR_ELLIPTIC)
		{
			projection_init(&Bi);
			projection_T0i_project(&pcls_cdm, &Bi, &phi);
			if (sim.baryon_flag)
				projection_T0i_project(&pcls_b, &Bi, &phi);
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (a >= 1. / (sim.z_switch_Bncdm[i] + 1.))
					projection_T0i_project(pcls_ncdm+i, &Bi, &phi);
			}
			projection_T0i_comm(&Bi);
		}

		projection_init(&Sij);
		projection_Tij_project(&pcls_cdm, &Sij, a, &phi);
		if (sim.baryon_flag)
			projection_Tij_project(&pcls_b, &Sij, a, &phi);
		if (a >= 1. / (sim.z_switch_linearchi + 1.))
		{
			for (i = 0; i < cosmo.num_ncdm; i++)
				projection_Tij_project(pcls_ncdm+i, &Sij, a, &phi);
		}
		projection_Tij_comm(&Sij);

#ifdef BENCHMARK
		projection_time += MPI_Wtime() - cycle_start_time;
		ref_time = MPI_Wtime();
#endif

if (sim.Kess_source_gravity==1)
{
// Kessence projection Tmunu
// In the projection zeta_integer comes, since synched with particles..
 	if (sim.vector_flag == VECTOR_ELLIPTIC)
		{
			projection_Tmunu_kessence( T00_Kess,T0i_Kess,Tij_Kess, dx, a, phi, phi_old, 	chi, pi_k, zeta_half, cosmo.Omega_kessence, cosmo.w_kessence, cosmo.cs2_kessence, Hconf(a, fourpiG, cosmo), fourpiG, sim.NL_kessence ,1 );
		}
 	else
		{
			projection_Tmunu_kessence( T00_Kess,T0i_Kess,Tij_Kess, dx, a, phi, phi_old, 	chi, pi_k, zeta_half, cosmo.Omega_kessence, cosmo.w_kessence, cosmo.cs2_kessence, Hconf(a, fourpiG, cosmo), fourpiG, sim.NL_kessence, 0 );
		}

		for (x.first(); x.test(); x.next())
		{
			// The coefficient is because it wanted to to be source according to eq C.2 of Gevolution paper
			// Note that it is multiplied to dx^2 and is divived by -a^3 because of definition of T00 which is scaled by a^3
			// We have T00 and Tij according to code's units, but source is important to calculate potentials and moving particles.
			// There is coefficient between Tij and Sij as source.
			source(x) += T00_Kess(x);
			if (sim.vector_flag == VECTOR_ELLIPTIC)for(int 	c=0;c<3;c++)Bi(x,c)+= (2. * fourpiG * dx * dx / a) * T0i_Kess(x,c);
			for(int c=0;c<6;c++)Sij(x,c)+=(2.) * Tij_Kess(x,c);
		}
}
#ifdef BENCHMARK
		kessence_update_time += MPI_Wtime() - ref_time;
		ref_time = MPI_Wtime();
#endif
// Kessence projection Tmunu end

		if (sim.gr_flag > 0)
		{
			T00hom = 0.;
			for (x.first(); x.test(); x.next())
				T00hom += source(x);
			parallel.sum<Real>(T00hom);
			T00hom /= (Real) numpts3d;

			if (cycle % CYCLE_INFO_INTERVAL == 0)
			{
				COUT << " cycle " << cycle << ", background information: z = " << (1./a) - 1. << ", average T00 = " << T00hom << ", background model = " << cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo) << endl;
			}

			if (dtau_old > 0.)
			{
				prepareFTsource<Real>(phi, chi, source, cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo), source, 3. * Hconf(a, fourpiG, cosmo) * dx * dx / dtau_old, fourpiG * dx * dx / a, 3. * Hconf(a, fourpiG, cosmo) * Hconf(a, fourpiG, cosmo) * dx * dx);  // prepare nonlinear source for phi update

#ifdef BENCHMARK
				ref2_time= MPI_Wtime();
#endif
				plan_source.execute(FFT_FORWARD);  // go to k-space
#ifdef BENCHMARK
				fft_time += MPI_Wtime() - ref2_time;
				fft_count++;
#endif

				solveModifiedPoissonFT(scalarFT, scalarFT, 1. / (dx * dx), 3. * Hconf(a, fourpiG, cosmo) / dtau_old);  // phi update (k-space)

#ifdef BENCHMARK
				ref2_time= MPI_Wtime();
#endif
				plan_phi.execute(FFT_BACKWARD);	 // go back to position space
#ifdef BENCHMARK
				fft_time += MPI_Wtime() - ref2_time;
				fft_count++;
#endif
			}
		}
		else
		{
#ifdef BENCHMARK
			ref2_time= MPI_Wtime();
#endif
			plan_source.execute(FFT_FORWARD);  // Newton: directly go to k-space
#ifdef BENCHMARK
			fft_time += MPI_Wtime() - ref2_time;
			fft_count++;
#endif

			solveModifiedPoissonFT(scalarFT, scalarFT, fourpiG / a);  // Newton: phi update (k-space)

#ifdef BENCHMARK
			ref2_time= MPI_Wtime();
#endif
			plan_phi.execute(FFT_BACKWARD);	 // go back to position space
#ifdef BENCHMARK
			fft_time += MPI_Wtime() - ref2_time;
			fft_count++;
#endif
		}

		phi.updateHalo();  // communicate halo values

		// record some background data
		if (kFT.setCoord(0, 0, 0))
		{
			sprintf(filename, "%s%s_background.dat", sim.output_path, sim.basename_generic);
			outfile = fopen(filename, "a");
			if (outfile == NULL)
			{
				cout << " error opening file for background output!" << endl;
			}
			else
			{
				if (cycle == 0)
					fprintf(outfile, "# background statistics\n# cycle   tau/boxsize    a             conformal H/H0         Hconf_prime       phi(k=0)       T00(k=0)\n");
				fprintf(outfile, " %6d   %e   %e   %e   %e   %e   %e\n", cycle, tau, a, Hconf(a, fourpiG, cosmo) / Hconf(1., fourpiG, cosmo),Hconf_prime(a_kess, fourpiG, cosmo), scalarFT(kFT).real(), T00hom);
				fclose(outfile);
			}
		}
		// done recording background data

		prepareFTsource<Real>(phi, Sij, Sij, 2. * fourpiG * dx * dx / a);  // prepare nonlinear source for additional equations

#ifdef BENCHMARK
		ref2_time= MPI_Wtime();
#endif
		plan_Sij.execute(FFT_FORWARD);  // go to k-space
#ifdef BENCHMARK
		fft_time += MPI_Wtime() - ref2_time;
		fft_count += 6;
#endif

#ifdef HAVE_CLASS
		if (sim.radiation_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.))
		{
			prepareFTchiLinear(class_background, class_perturbs, class_spectra, scalarFT, sim, ic, cosmo, fourpiG, a);
			projectFTscalar(SijFT, scalarFT, 1);
		}
		else
#endif
		projectFTscalar(SijFT, scalarFT);  // construct chi by scalar projection (k-space)

#ifdef BENCHMARK
		ref2_time= MPI_Wtime();
#endif
		plan_chi.execute(FFT_BACKWARD);	 // go back to position space
#ifdef BENCHMARK
		fft_time += MPI_Wtime() - ref2_time;
		fft_count++;
#endif
		chi.updateHalo();  // communicate halo values

		if (sim.vector_flag == VECTOR_ELLIPTIC)
		{
#ifdef BENCHMARK
			ref2_time= MPI_Wtime();
#endif
			plan_Bi.execute(FFT_FORWARD);
#ifdef BENCHMARK
			fft_time += MPI_Wtime() - ref2_time;
			fft_count++;
#endif
			projectFTvector(BiFT, BiFT, fourpiG * dx * dx); // solve B using elliptic constraint (k-space)
#ifdef CHECK_B
			evolveFTvector(SijFT, BiFT_check, a * a * dtau_old);
#endif
		}
		else
			evolveFTvector(SijFT, BiFT, a * a * dtau_old);  // evolve B using vector projection (k-space)

		if (sim.gr_flag > 0)
		{
#ifdef BENCHMARK
			ref2_time= MPI_Wtime();
#endif
			plan_Bi.execute(FFT_BACKWARD);  // go back to position space
#ifdef BENCHMARK
			fft_time += MPI_Wtime() - ref2_time;
			fft_count += 3;
#endif
			Bi.updateHalo();  // communicate halo values
		}

#ifdef BENCHMARK
		gravity_solver_time += MPI_Wtime() - ref_time;
		ref_time = MPI_Wtime();
#endif




		// snapshot output
		if (snapcount < sim.num_snapshot && 1. / a < sim.z_snapshot[snapcount] + 1.)
		{
			COUT << COLORTEXT_CYAN << " writing snapshot" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;

#ifdef CHECK_B
			//kessence included
			writeSnapshots(sim, cosmo, fourpiG, hdr, a, snapcount, h5filename, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &pi_k, &zeta_half, &chi, &Bi, &T00_Kess, &T0i_Kess, &Tij_Kess, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, &Bi_check, &BiFT_check, &plan_Bi_check);

#else
			//kessence included
			writeSnapshots(sim, cosmo, fourpiG, hdr, a, snapcount, h5filename, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &pi_k, &zeta_half, &chi, &Bi, &T00_Kess, &T0i_Kess, &Tij_Kess, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
			#endif

			snapcount++;
		}

#ifdef BENCHMARK
		snapshot_output_time += MPI_Wtime() - ref_time;
		ref_time = MPI_Wtime();
#endif

		// power spectra
		if (pkcount < sim.num_pk && 1. / a < sim.z_pk[pkcount] + 1.)
		{
			COUT << COLORTEXT_CYAN << " writing power spectra" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;

#ifdef CHECK_B
			//kessence included
			writeSpectra(sim, cosmo, fourpiG, a, pkcount, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &pi_k, &zeta_half, &chi, &Bi, &T00_Kess, &T0i_Kess, &Tij_Kess, &source, &Sij, &scalarFT ,&scalarFT_pi, &scalarFT_zeta_half, &BiFT, &T00_KessFT, &T0i_KessFT, &Tij_KessFT, &SijFT, &plan_phi, &plan_pi_k, &plan_zeta_half, &plan_chi, &plan_Bi, &plan_T00_Kess, &plan_T0i_Kess, &plan_Tij_Kess, &plan_source, &plan_Sij, &Bi_check, &BiFT_check, &plan_Bi_check);
#else
			//kessence included
			writeSpectra(sim, cosmo, fourpiG, a, pkcount, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &pi_k, &zeta_half, &chi, &Bi, &T00_Kess, &T0i_Kess, &Tij_Kess, &source, &Sij, &scalarFT ,&scalarFT_pi, &scalarFT_zeta_half, &BiFT, &T00_KessFT, &T0i_KessFT, &Tij_KessFT, &SijFT, &plan_phi, &plan_pi_k, &plan_zeta_half, &plan_chi, &plan_Bi, &plan_T00_Kess, &plan_T0i_Kess, &plan_Tij_Kess, &plan_source, &plan_Sij);

      writeSpectra_phi_prime(sim, cosmo, fourpiG, a, pkcount, &phi_prime, &phi_prime_scalarFT, &phi_prime_plan);


#endif

			pkcount++;
		}

    // cout<<"EXACT_OUTPUT_REDSHIFTS: "<<EXACT_OUTPUT_REDSHIFTS<<endl;
    #ifdef EXACT_OUTPUT_REDSHIFTS
    		tmp = a;
    		rungekutta4bg(tmp, fourpiG, cosmo, 0.5 * dtau);
    		rungekutta4bg(tmp, fourpiG, cosmo, 0.5 * dtau);

    		if (pkcount < sim.num_pk && 1. / tmp < sim.z_pk[pkcount] + 1.)
    		{
    			writeSpectra(sim, cosmo, fourpiG, a, pkcount,
    #ifdef HAVE_CLASS
    					class_background, class_perturbs, class_spectra, ic,
    #endif
    					&pcls_cdm, &pcls_b, pcls_ncdm, &phi,&pi_k, &zeta_half, &chi, &Bi,&T00_Kess, &T0i_Kess, &Tij_Kess, &source, &Sij, &scalarFT ,&scalarFT_pi, &scalarFT_zeta_half, &BiFT, &T00_KessFT, &T0i_KessFT, &Tij_KessFT, &SijFT, &plan_phi, &plan_pi_k, &plan_zeta_half, &plan_chi, &plan_Bi, &plan_T00_Kess, &plan_T0i_Kess, &plan_Tij_Kess, &plan_source, &plan_Sij
    #ifdef CHECK_B
    					, &Bi_check, &BiFT_check, &plan_Bi_check
    #endif
		    );
    writeSpectra_phi_prime(sim, cosmo, fourpiG, a, pkcount, &phi_prime, &phi_prime_scalarFT, &phi_prime_plan);


    		}
    #endif // EXACT_OUTPUT_REDSHIFTS


#ifdef BENCHMARK
		spectra_output_time += MPI_Wtime() - ref_time;
#endif

		if (pkcount >= sim.num_pk && snapcount >= sim.num_snapshot) break; // simulation complete

		// compute number of step subdivisions for particle updates
		numsteps = 1;
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			if (dtau * maxvel[i+1+sim.baryon_flag] > dx * sim.movelimit)
				numsteps_ncdm[i] = (int) ceil(dtau * maxvel[i+1+sim.baryon_flag] / dx / sim.movelimit);
			else numsteps_ncdm[i] = 1;

			if (numsteps < numsteps_ncdm[i]) numsteps = numsteps_ncdm[i];
		}
		if (numsteps > 1 && numsteps % 2 > 0) numsteps++;   // if >1, make it an even number

		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			if (numsteps / numsteps_ncdm[i] <= 1) numsteps_ncdm[i] = numsteps;
			else if (numsteps_ncdm[i] > 1) numsteps_ncdm[i] = numsteps / 2;
		}

		if (cycle % CYCLE_INFO_INTERVAL == 0)
		{
			COUT << " cycle " << cycle << ", time integration information: max |v| = " << maxvel[0] << " (cdm Courant factor = " << maxvel[0] * dtau / dx;
			if (sim.baryon_flag)
			{
				COUT << "), baryon max |v| = " << maxvel[1] << " (Courant factor = " << maxvel[1] * dtau / dx;
			}

			COUT << "), time step / Hubble time = " << Hconf(a, fourpiG, cosmo) * dtau;

			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (i == 0)
				{
					COUT << endl << " time step subdivision for ncdm species: ";
				}
				COUT << numsteps_ncdm[i] << " (max |v| = " << maxvel[i+1+sim.baryon_flag] << ")";
				if (i < cosmo.num_ncdm-1)
				{
					COUT << ", ";
				}
			}

			COUT << endl;
		}

		//Kessence
#ifdef BENCHMARK
		ref_time = MPI_Wtime();
#endif
        for (x.first(); x.test(); x.next())
    		{
    			phi_prime(x) =(phi(x)-phi_old(x))/(dtau);
    		}
        // We just need to update halo when we want to calculate spatial derivative or use some neibours at the same time! So here wo do not nee to update halo for phi_prime!
//**********************
//Kessence - LeapFrog:START
//**********************
  double a_kess=a;
  //First we update zeta_half to have it at -1/2 just in the first loop
  if(cycle==0)
  {
    for (i=0;i<sim.nKe_numsteps;i++)
    {
      //computing zeta_half(-1/2) and zeta_int(-1) but we do not work with zeta(-1)
      update_zeta(-dtau/ (2. * sim.nKe_numsteps) , dx, a_kess, phi, phi_old, chi, chi_old, pi_k, zeta_half, cosmo.Omega_kessence, cosmo.w_kessence, cosmo.cs2_kessence, Hconf(a_kess, fourpiG, cosmo), Hconf_prime(a_kess, fourpiG, cosmo), sim.NL_kessence);
      // zeta_integer.updateHalo();
      zeta_half.updateHalo();
    }
  }

 //Then fwe start the main loop zeta is updated to get zeta(n+1/2) from pi(n) and zeta(n-1/2)
	for (i=0;i<sim.nKe_numsteps;i++)
	{

    //********************************************************************************
    //Updating zeta_integer to get zeta_integer(n+1/2) and zeta_integer(n+1), in the first loop is getting zeta_integer(1/2) and zeta_integer(1)
    // In sum: zeta_integer(n+1/2) = zeta_integer(n-1/2)+ zeta_integer'(n)dtau which needs background to be at n with then
    //Note that here for zeta_integer'(n) we need background to be at n and no need to update it.
    //\zeta_integer(n+1/2) = \zeta_integer(n-1/2) + \zeta_integer'(n)  dtau
    //We also update zeta_int from n to n+1
    //********************************************************************************
    update_zeta(dtau/ sim.nKe_numsteps, dx, a_kess, phi, phi_old, chi, chi_old, pi_k, zeta_half, cosmo.Omega_kessence, cosmo.w_kessence, cosmo.cs2_kessence, Hconf(a_kess, fourpiG, cosmo), Hconf_prime(a_kess, fourpiG, cosmo), sim.NL_kessence);
    // zeta_integer.updateHalo();
    zeta_half.updateHalo();
    //********************************************************************************
    //Since we have pi(n+1)=pi(n) + pi'(n+1/2), and in pi'(n+1/2) we have H(n+1/2) we update the background before updating the pi to have H(n+1/2), Moreover zeta(n+1) = zeta(n+1/2) + zeta'(n+1/2), so we put zeta_int updating in the pi updating!
    //********************************************************************************
    rungekutta4bg(a_kess, fourpiG, cosmo,  dtau  / sim.nKe_numsteps / 2.0);
    //********************************************************************************
    //we update pi to have it at n+1 (at first loop from the value at (0) and the value of zeta_integer at 1/2 and H(n+1/2) we update pi at (1))
    //In the pi update we also update zeta_int because we need the values of a_kess and H_kess at step n+1/2
    //By the below update we get pi(n+1) and zeta(n+1)
    //********************************************************************************
    update_pi_k(dtau/ sim.nKe_numsteps, dx, a_kess, phi, phi_old, chi, chi_old, pi_k, zeta_half, cosmo.Omega_kessence, cosmo.w_kessence, cosmo.cs2_kessence, Hconf(a_kess, fourpiG, cosmo), Hconf_prime(a_kess, fourpiG, cosmo), sim.NL_kessence); // H_old is updated here in the function
		pi_k.updateHalo();

    //********************************************************************************
    // Now we have pi(n+1) and a_kess(n+1/2) so we update background by halfstep to have a_kess(n+1)
    //********************************************************************************
    rungekutta4bg(a_kess, fourpiG, cosmo,  dtau  / sim.nKe_numsteps / 2.0 );

	}
#ifdef BENCHMARK
    kessence_update_time += MPI_Wtime() - ref_time;
    ref_time = MPI_Wtime();
#endif
//**********************
//Kessence - LeapFrog: End
//**********************

		for (j = 0; j < numsteps; j++) // particle update
		{
#ifdef BENCHMARK
			ref2_time = MPI_Wtime();
#endif
			f_params[0] = a;
			f_params[1] = a * a * sim.numpts;
			if (j == 0)
			{
				if (sim.gr_flag > 0)
				{
					maxvel[0] = pcls_cdm.updateVel(update_q, (dtau + dtau_old) / 2., update_cdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
					if (sim.baryon_flag)
						maxvel[1] = pcls_b.updateVel(update_q, (dtau + dtau_old) / 2., update_b_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
				}
				else
				{
					maxvel[0] = pcls_cdm.updateVel(update_q_Newton, (dtau + dtau_old) / 2., update_cdm_fields, ((sim.radiation_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);
					if (sim.baryon_flag)
						maxvel[1] = pcls_b.updateVel(update_q_Newton, (dtau + dtau_old) / 2., update_b_fields, ((sim.radiation_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);
				}

#ifdef BENCHMARK
				update_q_count++;
#endif
			}

			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (j % (numsteps / numsteps_ncdm[i]) == 0)
				{
					if (sim.gr_flag > 0)
						maxvel[i+1+sim.baryon_flag] = pcls_ncdm[i].updateVel(update_q, (dtau + dtau_old) / 2. / numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
					else
						maxvel[i+1+sim.baryon_flag] = pcls_ncdm[i].updateVel(update_q_Newton, (dtau + dtau_old) / 2. / numsteps_ncdm[i], update_ncdm_fields, ((sim.radiation_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);

#ifdef BENCHMARK
					update_q_count++;
#endif
				}
			}
#ifdef BENCHMARK
			update_q_time += MPI_Wtime() - ref2_time;
			ref2_time = MPI_Wtime();
#endif

			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (numsteps > 1 && ((numsteps_ncdm[i] == 1 && j == numsteps / 2) || (numsteps_ncdm[i] == numsteps / 2 && j % 2 > 0)))
				{
					if (sim.gr_flag > 0)
						pcls_ncdm[i].moveParticles(update_pos, dtau / numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
					else
						pcls_ncdm[i].moveParticles(update_pos_Newton, dtau / numsteps_ncdm[i], NULL, 0, f_params);
#ifdef BENCHMARK
						moveParts_count++;
						moveParts_time += MPI_Wtime() - ref2_time;
						ref2_time = MPI_Wtime();
#endif
				}
			}

			if (numsteps == 1)
				rungekutta4bg(a, fourpiG, cosmo, 0.5 * dtau / numsteps);  // evolve background by half a time step

			f_params[0] = a;
			f_params[1] = a * a * sim.numpts;
			if (numsteps == 1 || j == numsteps / 2)
			{
				if (sim.gr_flag > 0)
				{
					pcls_cdm.moveParticles(update_pos, dtau, update_cdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 0), f_params);
					if (sim.baryon_flag)
						pcls_b.moveParticles(update_pos, dtau, update_b_fields, (1. / a < ic.z_relax + 1. ? 3 : 0), f_params);
				}
				else
				{
					pcls_cdm.moveParticles(update_pos_Newton, dtau, NULL, 0, f_params);
					if (sim.baryon_flag)
						pcls_b.moveParticles(update_pos_Newton, dtau, NULL, 0, f_params);
				}

#ifdef BENCHMARK
				moveParts_count++;
				moveParts_time += MPI_Wtime() - ref2_time;
				ref2_time = MPI_Wtime();
#endif
			}

			if (numsteps != 1)
				rungekutta4bg(a, fourpiG, cosmo, 0.5 * dtau / numsteps);  // evolve background by half a time step

			f_params[0] = a;
			f_params[1] = a * a * sim.numpts;
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (numsteps_ncdm[i] == numsteps)
				{
					if (sim.gr_flag > 0)
						pcls_ncdm[i].moveParticles(update_pos, dtau / numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
					else
						pcls_ncdm[i].moveParticles(update_pos_Newton, dtau / numsteps_ncdm[i], NULL, 0, f_params);
#ifdef BENCHMARK
						moveParts_count++;
						moveParts_time += MPI_Wtime() - ref2_time;
						ref2_time = MPI_Wtime();
#endif
				}
			}

			rungekutta4bg(a, fourpiG, cosmo, 0.5 * dtau / numsteps);  // evolve background by half a time step
		}   // particle update done

		parallel.max<double>(maxvel, numspecies);

		if (sim.gr_flag > 0)
		{
			for (i = 0; i < numspecies; i++)
				maxvel[i] /= sqrt(maxvel[i] * maxvel[i] + 1.0);
		}

		tau += dtau;

		if (tau_Lambda < 0. && (cosmo.Omega_m / a / a / a) < cosmo.Omega_Lambda)
		{
			tau_Lambda = tau;
			COUT << "matter-dark energy equality at z=" << ((1./a) - 1.) << endl;
		}

		if (sim.wallclocklimit > 0.)   // check for wallclock time limit
		{
			tmp = MPI_Wtime() - start_time;
			parallel.max(tmp);
			if (tmp > sim.wallclocklimit)   // hibernate
			{
				COUT << COLORTEXT_YELLOW << " reaching hibernation wallclock limit, hibernating..." << COLORTEXT_RESET << endl;
				COUT << COLORTEXT_CYAN << " writing hibernation point" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;
				if (sim.vector_flag == VECTOR_PARABOLIC && sim.gr_flag == 0)
					plan_Bi.execute(FFT_BACKWARD);
#ifdef CHECK_B
				if (sim.vector_flag == VECTOR_ELLIPTIC)
				{
					plan_Bi_check.execute(FFT_BACKWARD);
					hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, pi_k, zeta_half, chi, Bi_check, a, tau, dtau, cycle);
				}
				else
#endif
				hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, pi_k, zeta_half, chi, Bi, a, tau, dtau, cycle);
				break;
			}
		}

		if (restartcount < sim.num_restart && 1. / a < sim.z_restart[restartcount] + 1.)
		{
			COUT << COLORTEXT_CYAN << " writing hibernation point" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;
			if (sim.vector_flag == VECTOR_PARABOLIC && sim.gr_flag == 0)
				plan_Bi.execute(FFT_BACKWARD);
#ifdef CHECK_B
			if (sim.vector_flag == VECTOR_ELLIPTIC)
			{
				plan_Bi_check.execute(FFT_BACKWARD);
				hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, pi_k, zeta_half, chi, Bi, a, tau, dtau, cycle, restartcount);
			}
			else
#endif
			hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, pi_k, zeta_half, chi, Bi, a, tau, dtau, cycle, restartcount);
			restartcount++;
		}

		dtau_old = dtau;

		if (sim.Cf * dx < sim.steplimit / Hconf(a, fourpiG, cosmo))
			dtau = sim.Cf * dx;
		else
			dtau = sim.steplimit / Hconf(a, fourpiG, cosmo);

		cycle++;

#ifdef BENCHMARK
		cycle_time += MPI_Wtime()-cycle_start_time;
#endif
	}

	COUT << COLORTEXT_GREEN << " simulation complete." << COLORTEXT_RESET << endl;

#ifdef HAVE_CLASS
	if (sim.radiation_flag > 0)
		freeCLASSstructures(class_background, class_perturbs, class_spectra);
#endif

#ifdef BENCHMARK
	run_time = MPI_Wtime() - start_time;

	parallel.sum(run_time);
	parallel.sum(cycle_time);
	parallel.sum(projection_time);
	parallel.sum(snapshot_output_time);
	parallel.sum(spectra_output_time);
	parallel.sum(gravity_solver_time);
  parallel.sum(kessence_update_time);
	parallel.sum(fft_time);
	parallel.sum(update_q_time);
	parallel.sum(moveParts_time);

	COUT << endl << "BENCHMARK" << endl;
	COUT << "total execution time  : "<<hourMinSec(run_time) << endl;
	COUT << "total number of cycles: "<< cycle << endl;
	COUT << "time consumption breakdown:" << endl;
	COUT << "initialization   : "  << hourMinSec(initialization_time) << " ; " << 100. * initialization_time/run_time <<"%."<<endl;
	COUT << "main loop        : "  << hourMinSec(cycle_time) << " ; " << 100. * cycle_time/run_time <<"%."<<endl;

	COUT << "----------- main loop: components -----------"<<endl;

	COUT << "projections                : "<< hourMinSec(projection_time) << " ; " << 100. * projection_time/cycle_time <<"%."<<endl;
  //Kessence update
  COUT << "Kessence_update                : "<< hourMinSec(kessence_update_time) << " ; " << 100. * kessence_update_time/cycle_time <<"%."<<endl;
	COUT << "snapshot outputs           : "<< hourMinSec(snapshot_output_time) << " ; " << 100. * snapshot_output_time/cycle_time <<"%."<<endl;
	COUT << "power spectra outputs      : "<< hourMinSec(spectra_output_time) << " ; " << 100. * spectra_output_time/cycle_time <<"%."<<endl;
	COUT << "update momenta (count: "<<update_q_count <<"): "<< hourMinSec(update_q_time) << " ; " << 100. * update_q_time/cycle_time <<"%."<<endl;
	COUT << "move particles (count: "<< moveParts_count <<"): "<< hourMinSec(moveParts_time) << " ; " << 100. * moveParts_time/cycle_time <<"%."<<endl;
	COUT << "gravity solver             : "<< hourMinSec(gravity_solver_time) << " ; " << 100. * gravity_solver_time/cycle_time <<"%."<<endl;
	COUT << "-- thereof Fast Fourier Transforms (count: " << fft_count <<"): "<< hourMinSec(fft_time) << " ; " << 100. * fft_time/gravity_solver_time <<"%."<<endl;
#endif

#ifdef EXTERNAL_IO
		ioserver.stop();
	}
#endif

	return 0;
}
