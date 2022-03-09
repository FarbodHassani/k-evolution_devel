//////////////////////////
// Copyright (c) 2015-2019 Julian Adamek
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
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: May 2019

//
//////////////////////////
#include <stdlib.h>
#include <set>
#include <vector>
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
#include "tools.hpp"
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
#include "output.hpp"
#include "hibernation.hpp"
#ifdef VELOCITY
#include "velocity.hpp"
#endif

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
	double lightcone_output_time = 0;
	double gravity_solver_time = 0;
	double fft_time = 0;
	int fft_count = 0;
	double update_q_time = 0;
	int update_q_count = 0;
	double moveParts_time = 0;
	int  moveParts_count =0;
	//kessence
	double a_kess;
  double Hc;


#endif  //BENCHMARK

	int n = 0, m = 0;
	int io_size = 0;
	int io_group_size = 0;

	int i, j, cycle = 0, snapcount = 0, snapcount_b=1, pkcount = 0, restartcount = 0, usedparams, numparam = 0, numsteps, numspecies, done_hij;
	int numsteps_ncdm[MAX_PCL_SPECIES-2];
	long numpts3d;
	int box[3];
	double dtau, dtau_old, dx, tau, a, fourpiG, tmp, start_time;
	double maxvel[MAX_PCL_SPECIES];
	FILE * outfile;
	char filename[2*PARAM_MAX_LENGTH+24];
	string h5filename;
	char * settingsfile = NULL;
	char * precisionfile = NULL;
	parameter * params = NULL;
	metadata sim;
	cosmology cosmo;
	icsettings ic;
	double T00hom;
  #ifdef HAVE_CLASS_BG
	gsl_interp_accel * acc = gsl_interp_accel_alloc();
	//Background variables EFTevolution //TODO_EB: add as many as necessary
	gsl_spline * H_spline = NULL;
  gsl_spline * phi_smg = NULL;
  gsl_spline * phi_smg_prime = NULL;

  // gsl_spline * phi_smg_prime_prime = NULL;
	gsl_spline * cs2_spline = NULL;
	gsl_spline * rho_smg_spline = NULL;
	gsl_spline * p_smg_spline = NULL;
	gsl_spline * rho_crit_spline = NULL;
	#endif

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
			case 'p':
#ifndef HAVE_CLASS
				cout << "HAVE_CLASS needs to be set at compilation to use CLASS precision files" << endl;
				exit(-100);
#endif
				precisionfile = argv[++i];
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
	if (!io_size || !io_group_size)
	{
		cout << "invalid number of I/O tasks and group sizes for I/O server (-DEXTERNAL_IO)" << endl;
		exit(-1000);
	}
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
  thermo class_thermo;
  perturbs class_perturbs;

  	if (precisionfile != NULL)
	  	numparam = loadParameterFile(precisionfile, params);
	else
#endif
		numparam = 0;

#ifdef HAVE_CLASS_BG
  //TODO_EB:add BG functions here
  initializeCLASSstructures(sim, ic, cosmo, class_background, class_thermo, class_perturbs, params, numparam);
  loadBGFunctions(class_background, H_spline, "H [1/Mpc]", sim.z_in);
  loadBGFunctions(class_background, phi_smg, "phi_smg", sim.z_in);
  loadBGFunctions(class_background, phi_smg_prime, "phi_prime_smg", sim.z_in);
  // loadBGFunctions(class_background, phi_smg_prime_prime, "phi_smg_prime_prime", sim.z_in);
  loadBGFunctions(class_background, cs2_spline, "c_s^2", sim.z_in);
  loadBGFunctions(class_background, rho_smg_spline, "(.)rho_smg", sim.z_in);
  loadBGFunctions(class_background, p_smg_spline, "(.)p_smg", sim.z_in);
  loadBGFunctions(class_background, rho_crit_spline, "(.)rho_crit", sim.z_in);

#endif

	h5filename.reserve(2*PARAM_MAX_LENGTH);
	h5filename.assign(sim.output_path);

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
	set<long> IDbacklog[MAX_PCL_SPECIES];

	Field<Real> phi;
	Field<Real> source;
	Field<Real> chi;
	Field<Real> chi_old;
  Field<Real> phi_old;
	Field<Real> Sij;
	Field<Real> Bi;
	Field<Cplx> scalarFT;
	Field<Cplx> SijFT;
	Field<Cplx> BiFT;
	source.initialize(lat,1);
	phi.initialize(lat,1);

	//kessence
  //phi at two step before to compute phi'(n+1/2)
  #ifdef BACKREACTION_TEST
  Field<Real> short_wave;
  Field<Real> relativistic_term;
  Field<Real> stress_tensor;
  Field<Real> pi_old;
  Field<Real> det_gamma;
  Field<Real> Gradpi_Gradpi;
  Field<Real> det_gamma_old;
  Field<Real> pi_prime_old;
  Field<Real> pi_prime_dot;
  Field<Real> cs2_full;

  Field<Cplx> scalarFT_pi_old;
  Field<Cplx> scalarFT_pi_prime_old;
  Field<Cplx> scalarFT_pi_prime_dot;
  Field<Cplx> scalarFT_det_gamma;
  Field<Cplx> scalarFT_Gradpi_Gradpi;
  Field<Cplx> scalarFT_det_gamma_old;
  Field<Cplx> scalarFT_cs2_full;
  #endif

  Field<Real> pi;
  Field<Real> pi_prime;

  // Runge Kutta 4th order coefficients
  Field<Real> l0; // Update  pi_v
  Field<Real> k0; // Update  pi
  Field<Real> l1;
  Field<Real> k1;
  Field<Real> l2;
  Field<Real> k2;
  Field<Real> l3;
  Field<Real> k3;
  //
	Field<Real> T00_Kess;
	Field<Real> T0i_Kess;
	Field<Real> Tij_Kess;
	// Field<Cplx> scalarFT_phi_old;
  #ifdef BACKREACTION_TEST
  Field<Cplx> short_wave_scalarFT;
  Field<Cplx> relativistic_term_scalarFT;
  Field<Cplx> stress_tensor_scalarFT;
  #endif
	// Field<Cplx> scalarFT_chi_old;
  Field<Cplx> scalarFT_pi;
  // Field<Cplx> scalarFT_zeta;
  Field<Cplx> scalarFT_pi_prime;
  Field<Cplx> scalarFT_phi_old;
  Field<Cplx> scalarFT_chi_old;


  Field<Cplx> scalarFT_l0;
  Field<Cplx> scalarFT_l1;
  Field<Cplx> scalarFT_l2;
  Field<Cplx> scalarFT_l3;
  Field<Cplx> scalarFT_k0;
  Field<Cplx> scalarFT_k1;
  Field<Cplx> scalarFT_k2;
  Field<Cplx> scalarFT_k3;

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
  #ifdef VELOCITY
  	Field<Real> vi;
  	Field<Cplx> viFT;
  	vi.initialize(lat,3);
  	viFT.initialize(latFT,3);
  	PlanFFT<Cplx> plan_vi(&vi, &viFT);
  	double a_old;
  #endif
  #ifdef BACKREACTION_TEST
  if(parallel.isRoot()) cout << "\033[1;31m The blowup tests are requested\033[0m\n";
    pi_old.initialize(lat,1);
  	scalarFT_pi_old.initialize(latFT,1);
  	PlanFFT<Cplx> plan_pi_old(&pi_old, &scalarFT_pi_old);

    pi_prime_old.initialize(lat,1);
    scalarFT_pi_prime_old.initialize(latFT,1);
    PlanFFT<Cplx> plan_pi_prime_old(&pi_prime_old, &scalarFT_pi_prime_old);

    // zeta.initialize(lat,1);
    // scalarFT_zeta.initialize(latFT,1);
    // PlanFFT<Cplx> plan_zeta(&zeta, &scalarFT_zeta);

    pi_prime_dot.initialize(lat,1);
    scalarFT_pi_prime_dot.initialize(latFT,1);
    PlanFFT<Cplx> plan_pi_prime_dot(&pi_prime_dot, &scalarFT_pi_prime_dot);

    det_gamma.initialize(lat,1);
    scalarFT_det_gamma.initialize(latFT,1);
    PlanFFT<Cplx> plan_det_gamma(&det_gamma, &scalarFT_det_gamma);


    Gradpi_Gradpi.initialize(lat,1);
    scalarFT_Gradpi_Gradpi.initialize(latFT,1);
    PlanFFT<Cplx> plan_Gradpi_Gradpi(&Gradpi_Gradpi, &scalarFT_Gradpi_Gradpi);

    det_gamma_old.initialize(lat,1);
    scalarFT_det_gamma_old.initialize(latFT,1);
    PlanFFT<Cplx> plan_det_gamma_old(&det_gamma_old, &scalarFT_det_gamma_old);

    cs2_full.initialize(lat,1);
    scalarFT_cs2_full.initialize(latFT,1);
    PlanFFT<Cplx> plan_cs2_full(&cs2_full, &scalarFT_cs2_full);

  #endif

  phi_old.initialize(lat,1);
	scalarFT_phi_old.initialize(latFT,1);
	PlanFFT<Cplx> plan_phi_old(&phi_old, &scalarFT_phi_old);

  chi_old.initialize(lat,1);
  scalarFT_chi_old.initialize(latFT,1);
  PlanFFT<Cplx> plan_chi_old(&chi_old, &scalarFT_chi_old);

	//Kessence part initializing
	// phi_old.initialize(lat,1);
	// scalarFT_phi_old.initialize(latFT,1);
	// PlanFFT<Cplx> plan_phi_old(&phi_old, &scalarFT_phi_old);
  //Relativistic corrections
  // #ifdef BACKREACTION_TEST
  // short_wave.initialize(lat,1);
  // short_wave_scalarFT.initialize(latFT,1);
  // PlanFFT<Cplx> short_wave_plan(&short_wave, &short_wave_scalarFT);
  // relativistic_term.initialize(lat,1);
  // relativistic_term_scalarFT.initialize(latFT,1);
  // PlanFFT<Cplx> relativistic_term_plan(&relativistic_term, &relativistic_term_scalarFT);
  // stress_tensor.initialize(lat,1);
  // stress_tensor_scalarFT.initialize(latFT,1);
  // PlanFFT<Cplx> stress_tensor_plan(&stress_tensor, &stress_tensor_scalarFT);
  // #endif
	//pi kessence
	pi.initialize(lat,1);
	scalarFT_pi.initialize(latFT,1);
	PlanFFT<Cplx> plan_pi(&pi, &scalarFT_pi);
  pi_prime.initialize(lat,1);
  scalarFT_pi_prime.initialize(latFT,1);
  PlanFFT<Cplx> plan_pi_prime(&pi_prime, &scalarFT_pi_prime);

  Field<Real> deltaX;
  Field<Cplx> scalarFT_deltaX;
  deltaX.initialize(lat,1); // deltaX = X - Xhat
  scalarFT_deltaX.initialize(latFT,1);
  PlanFFT<Cplx> plan_deltaX(&deltaX, &scalarFT_deltaX);

  l0.initialize(lat,1);
  scalarFT_l0.initialize(latFT,1);
  PlanFFT<Cplx> plan_l0(&l0, &scalarFT_l0);

  k0.initialize(lat,1);
  scalarFT_k0.initialize(latFT,1);
  PlanFFT<Cplx> plan_k0(&k0, &scalarFT_k0);

  l1.initialize(lat,1);
  scalarFT_l1.initialize(latFT,1);
  PlanFFT<Cplx> plan_l1(&l1, &scalarFT_l1);

  k1.initialize(lat,1);
  scalarFT_k1.initialize(latFT,1);
  PlanFFT<Cplx> plan_k1(&k1, &scalarFT_k1);

  l2.initialize(lat,1);
  scalarFT_l2.initialize(latFT,1);
  PlanFFT<Cplx> plan_l2(&l2, &scalarFT_l2);

  k2.initialize(lat,1);
  scalarFT_k2.initialize(latFT,1);
  PlanFFT<Cplx> plan_k2(&k2, &scalarFT_k2);

  l3.initialize(lat,1);
  scalarFT_l3.initialize(latFT,1);
  PlanFFT<Cplx> plan_l3(&l3, &scalarFT_l3);

  k3.initialize(lat,1);
  scalarFT_k3.initialize(latFT,1);
  PlanFFT<Cplx> plan_k3(&k3, &scalarFT_k3);
	//chi_old initialize
	// chi_old.initialize(lat,1);
	// scalarFT_chi_old.initialize(latFT,1);
	// PlanFFT<Cplx> plan_chi_old(&chi_old, &scalarFT_chi_old);
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
	fourpiG = 1.5 * sim.boxsize * sim.boxsize / C_SPEED_OF_LIGHT / C_SPEED_OF_LIGHT;
	a = 1. / (1. + sim.z_in);
  tau = particleHorizon(a, fourpiG,
    #ifdef HAVE_CLASS_BG
    gsl_spline_eval(H_spline, 1., acc), class_background
    #else
    cosmo
    #endif
  );

  if (sim.Cf * dx < sim.steplimit / Hconf(a, fourpiG,//TODO_EB
		#ifdef HAVE_CLASS_BG
			H_spline, acc
		#else
			cosmo
		#endif
	) )
    dtau = sim.Cf * dx;

	else
	    dtau = sim.steplimit / 	Hconf(a, fourpiG,//TODO_EB
			#ifdef HAVE_CLASS_BG
				H_spline, acc
			#else
				cosmo
			#endif
			);
	dtau_old = 0.;

	if (ic.generator == ICGEN_BASIC)
		generateIC_basic(sim, ic, cosmo, fourpiG, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &pi, &pi_prime, &chi, &Bi, &source, &Sij, &scalarFT, &scalarFT_pi, &scalarFT_pi_prime, &BiFT, &SijFT, &plan_phi, &plan_pi, &plan_pi_prime, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, params, numparam);
	// generates ICs on the fly
	else if (ic.generator == ICGEN_READ_FROM_DISK)
  {
    // readIC(sim, ic, cosmo, fourpiG, a, tau, dtau, dtau_old, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &chi, &Bi, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, cycle, snapcount, pkcount, restartcount, IDbacklog, params, numparam);
    COUT << " error: IC generator is wrongly chosen!- For this code you need to use the basic initial conditions (IC_basic) " << endl;
	  parallel.abortForce();
  }
#ifdef ICGEN_PREVOLUTION
COUT << " error: IC generator is wrongly chosen!- For this code you need to use the basic initial conditions (IC_basic) " << endl;
parallel.abortForce();
#endif
#ifdef ICGEN_FALCONIC
COUT << " error: IC generator is wrongly chosen!- For this code you need to use the basic initial conditions (IC_basic) " << endl;
	parallel.abortForce();
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
#ifdef VELOCITY
	a_old = a;
	projection_init(&vi);
#endif

#ifdef BENCHMARK
	initialization_time = MPI_Wtime() - start_time;
	parallel.sum(initialization_time);
	COUT << COLORTEXT_GREEN << " initialization complete." << COLORTEXT_RESET << " BENCHMARK: " << hourMinSec(initialization_time) << endl << endl;
#else
	COUT << COLORTEXT_GREEN << " initialization complete." << COLORTEXT_RESET << endl << endl;
#endif

// #ifdef HAVE_CLASS
// 	if (sim.radiation_flag > 0 || sim.fluid_flag > 0)
// 	{
//     initializeCLASSstructures(sim, ic, cosmo, class_background, class_thermo, class_perturbs, params, numparam);
// 		if (sim.gr_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.) && (ic.generator == ICGEN_BASIC || (ic.generator == ICGEN_READ_FROM_DISK && cycle == 0)))
// 		{
// 			prepareFTchiLinear(class_background, class_perturbs, scalarFT, sim, ic, cosmo, fourpiG, a);
// 			plan_source.execute(FFT_BACKWARD);
// 			for (x.first(); x.test(); x.next())
// 				chi(x) += source(x);
// 			chi.updateHalo();
// 		}
// 	}
// 	if (numparam > 0) free(params);
// #endif


#ifdef BACKREACTION_TEST

// In case we want to initialize the /IC ourselves
//FH
//Initial condition
if (cosmo.MGtheory == 1) // Only if we use the full theory we add the BG
{
  if (cosmo.IC_scf == 1) // hiclass initial codnitions for the scalar field
  {
    if(cycle==0)
    {
    if(parallel.isRoot())  cout << "\033[1;34mThe initial condition for the scalar field is made using hiclass: \033[0m"<<1<<endl;
    if(parallel.isRoot() & (cosmo.IC_scf==0) )  cout << "\033[1;34mThe amplitude of the initial condition is set as field(x) = phi_bg(hiclass) + Amplitude * field(x) where the Amplitude is \033[0m"<<cosmo.IC_amplitude<<endl;
    }
    for (x.first(); x.test(); x.next())
      {
        Gradpi_Gradpi(x)= 0.25 * (pi(x + 0)  - pi(x - 0)) * (pi(x + 0) - pi(x - 0)) / (dx * dx); // Gradpi_Gradpi
        Gradpi_Gradpi(x)+=0.25 * (pi(x + 1)  - pi(x - 1)) * (pi(x + 1) - pi(x - 1)) / (dx * dx); // Gradpi_Gradpi
        Gradpi_Gradpi(x)+=0.25 * (pi(x + 2)  - pi(x - 2)) * (pi(x + 2) - pi(x - 2)) / (dx * dx); // Gradpi_Gradpi
      }
    for (x.first(); x.test(); x.next())
      {
        pi_prime(x) += gsl_spline_eval(phi_smg_prime, 1. / (1. + sim.z_in), acc);
        deltaX(x) =(pi_prime(x) * pi_prime(x) - Gradpi_Gradpi(x))/2./a/a  - cosmo.X_hat; // definition of deltaX_ini from phi_ini = phi_bar_ini + delta phi_ini and d^ipi d_i pi
        pi(x) += gsl_spline_eval(phi_smg, 1. / (1. + sim.z_in), acc) * gsl_spline_eval(H_spline, 1.0, acc)/
        #ifdef HAVE_CLASS_BG
        Hconf(1.0, fourpiG, H_spline, acc)
        #else
        Hconf(1.0, fourpiG, cosmo)
        #endif
        ;// phi has dimension of time so we multiply by H0_class/H_0 gevolution
        // pi_prime(x) +=  gsl_spline_eval(phi_smg_prime, 1. / (1. + sim.z_in), acc);
        det_gamma(x) = 0.;
      }
    pi_prime.updateHalo();  // communicate halo values
    deltaX.updateHalo();  // communicate halo values
    pi.updateHalo();  // communicate halo values
  }

  else if (cosmo.IC_scf == 0) // If the initial conditions are assumed to be small initially
  {
    if(cycle==0)
    {
    if(parallel.isRoot())  cout << "\033[1;34mThe initial condition for the scalar field is made using hiclass: \033[0m"<<1<<endl;
    if(parallel.isRoot() & (cosmo.IC_scf==0) )  cout << "\033[1;34mThe amplitude of the initial condition is set as p_prime(x) = 0 and field(x) = phi_bg(hiclass) + Amplitude * field(x) where the Amplitude is \033[0m"<<cosmo.IC_amplitude<<endl;
    }

    for (x.first(); x.test(); x.next())
      {
        Gradpi_Gradpi(x)= 0.25 * (pi(x + 0)  - pi(x - 0)) * (pi(x + 0) - pi(x - 0)) / (dx * dx); // Gradpi_Gradpi
        Gradpi_Gradpi(x)+=0.25 * (pi(x + 1)  - pi(x - 1)) * (pi(x + 1) - pi(x - 1)) / (dx * dx); // Gradpi_Gradpi
        Gradpi_Gradpi(x)+=0.25 * (pi(x + 2)  - pi(x - 2)) * (pi(x + 2) - pi(x - 2)) / (dx * dx); // Gradpi_Gradpi
      }
    for (x.first(); x.test(); x.next())
      {
        // We multiplu the initial power by an amplitude
        pi_prime(x) = pi_prime(x) * cosmo.IC_amplitude + gsl_spline_eval(phi_smg_prime, 1. / (1. + sim.z_in), acc);
        deltaX(x) =(pi_prime(x) * pi_prime(x) - Gradpi_Gradpi(x))/2./a/a  - cosmo.X_hat; // definition of deltaX_ini from phi_ini = phi_bar_ini + delta phi_ini and d^ipi d_i pi
        // We multiplu the initial power by an amplitude
        pi(x) = pi(x) * cosmo.IC_amplitude + gsl_spline_eval(phi_smg, 1. / (1. + sim.z_in), acc) * gsl_spline_eval(H_spline, 1.0, acc)/
        #ifdef HAVE_CLASS_BG
        Hconf(1.0, fourpiG, H_spline, acc)
        #else
        Hconf(1.0, fourpiG, cosmo)
        #endif
        ;
        det_gamma(x) = 0.;
      }
    pi_prime.updateHalo();  // communicate halo values
    deltaX.updateHalo();  // communicate halo values
    pi.updateHalo();  // communicate halo values
  }
}
//   //****************************
//   //****SAVE DATA To test Backreaction
//   //****************************
  FILE* Result_avg;
  FILE* Result_real;
  FILE* Result_fourier;
  FILE* Result_max;
  FILE* Result_min;
  FILE* Redshifts;
  FILE* snapshots_file;

  char filename_avg[60];
  char filename_real[60];
  char filename_fourier[60];
  char filename_max[60];
  char filename_min[60];
  char filename_redshift[60];
  char filename_snapshot[60];

  snprintf(filename_avg, sizeof(filename_avg),"./output/Result_avg.txt");
  snprintf(filename_real, sizeof(filename_real),"./output/Result_real.txt");
  snprintf(filename_fourier, sizeof(filename_fourier),"./output/Result_fourier.txt");
  snprintf(filename_max, sizeof(filename_max),"./output/Results_max.txt");
  snprintf(filename_min, sizeof(filename_min),"./output/Results_min.txt");
  snprintf(filename_snapshot, sizeof(filename_snapshot),"./output/snapshots.txt");

  // ofstream out(filename_avg,ios::out);
  ofstream out_avg(filename_avg,ios::out);
  ofstream out_real(filename_real,ios::out);
  ofstream out_fourier(filename_fourier,ios::out);
  ofstream out_max(filename_max,ios::out);
  ofstream out_min(filename_min,ios::out);
  ofstream out_snapshots(filename_snapshot,ios::out);
  snapshots_file=fopen(filename_snapshot,"w");


  Result_avg=fopen(filename_avg,"w");
  Result_real=fopen(filename_real,"w");
  Result_fourier=fopen(filename_fourier,"w");
  Result_max=fopen(filename_max,"w");
  Result_min=fopen(filename_min,"w");
  snapshots_file=fopen(filename_snapshot,"w");

  out_avg<<"### The result of the average over time \n### d tau = "<< dtau<<endl;
  out_avg<<"### number of kessence update = "<<  sim.nKe_numsteps <<endl;
  out_avg<<"### initial time = "<< tau <<endl;
  if (cosmo.MGtheory == 0)
{
  out_avg<<"### 1- z(redshift)\t2- average(H0 pi)\t3- average (zeta)\t 4- w \t 5- average (phi)\t6-c_s^2  \t7-tau " <<endl;
}
else
{
  out_avg<<"### 1- z(redshift)\t2- average(H0 pi)\t3- average (zeta)\t 4- average (deta_gamma)\t 5- average (phi)\t6-avg(c_s2)  \t7-tau " <<endl;
}


  out_max<<"### The result of the maximum over time \n### d tau = "<< dtau<<endl;
  out_max<<"### number of kessence update = "<<  sim.nKe_numsteps <<endl;
  out_max<<"### initial time = "<< tau <<endl;
  out_max<<"### 1- z\t2- max(H pi)\t3- max (zeta)\t 4- max (phi)\t 5- max(c_s^2)   " <<endl;

  out_min<<"### The result of the maximum over time \n### d tau = "<< dtau<<endl;
  out_min<<"### number of kessence update = "<<  sim.nKe_numsteps <<endl;
  out_min<<"### initial time = "<< tau <<endl;
  out_min<<"### 1- z\t2- min(H pi)\t3- min(zeta)\t 4- min(phi)\t 5- min(c_s^2)   " <<endl;

  out_real<<"### The result of the average over time \n### d tau = "<< dtau<<endl;
  out_real<<"### number of kessence update = "<<  sim.nKe_numsteps <<endl;
  out_real<<"### initial time = "<< tau <<endl;
  out_real<<"### 1- tau\t2- pi(x)\t3-zeta(x)\t 4-x" <<endl;


  out_fourier<<"### The result of the average over time \n### d tau = "<< dtau<<endl;
  out_fourier<<"### number of kessence update = "<<  sim.nKe_numsteps <<endl;
  out_fourier<<"### initial time = "<< tau <<endl;
  out_fourier<<"### 1- tau\t 2- pi(k)\t\t3-zeta(k)\t\t4-|k|\t\t 5-vec{k} \t 6-|k|^2"<<endl;


  out_snapshots<<"### The result of the snapshots produced over time for blow-up \n### d tau = "<< dtau<<endl;
  out_snapshots<<"### number of kessence update = "<<  sim.nKe_numsteps <<endl;
  out_snapshots<<"### initial time = "<< tau <<endl;
  out_snapshots<<"### 1- tau\t2- z \t3- a\t 4- zeta_avg\t 5- avg_pi\t 6- avg_phi\t 7- avg_det_gamma\t 8- H_conf/H0 \t 9- snap_count \t 10-cs2_full"<<endl;


//defining the average
double avg_pi = 0.;
double avg_zeta = 0.;
double avg_phi = 0.;
double avg_det_gamma = 0.;
double avg_cs2_full = 1.;
double avg_pi_old = 0.;

double max_pi = 0.;
double max_pi_old = 0.;
double max_zeta = 0.;
double max_cs2 = 0.;
double max_phi = 0.;

double min_pi = 0.;
double min_pi_old = 0.;
double min_zeta = 0.;
double min_cs2 = 0.;
double min_phi = 0.;

int norm_kFT_squared = 0.;

// HDF5 outputs!
string str_filename ;
string str_filename2 ;
string str_filename3 ;
string str_filename4 ;
string str_filename5 ;
#endif


  outfile = fopen(filename, "w");

	while (true)    // main loop
	{
    //Kessence IC

// #ifdef BACKREACTION_TEST
//       //****************************
//       //****PRINTING AVERAGE OVER TIME
//       //****************************
//       // check_field(  pi_prime, 1. , " H pi", numpts3d);
//       avg_pi =average(  pi, Hconf(1.0, fourpiG,//TODO_EB
// 			#ifdef HAVE_CLASS_BG
// 				H_spline, acc
// 			#else
// 				cosmo
// 			#endif
// 				), numpts3d ) ;
//       avg_zeta =average( pi_prime,1., numpts3d ) ;
//       avg_phi =average(  phi , 1., numpts3d ) ;
//       avg_det_gamma =average(  det_gamma , 1., numpts3d ) ;
//       avg_cs2_full =average(  cs2_full , 1., numpts3d ) ;
//
//       max_pi =maximum(  pi, 1., numpts3d ) ;
//       max_pi_old =maximum(  pi_old, 1., numpts3d ) ;
//       max_zeta =maximum(  pi_prime,1., numpts3d ) ;
//       max_phi =maximum(  phi , 1., numpts3d ) ;
//
//       COUT << scientific << setprecision(8);
//       if(isnan(avg_pi))
//       {
//         parallel.abortForce();
//       }
//         if (cosmo.MGtheory == 0)
//       {
//         cout<<"z = "<<1./a - 1<<" p_smg_spline:"<<gsl_spline_eval(p_smg_spline, a, acc)<<" rho_smg_spline:"<<gsl_spline_eval(rho_smg_spline, a, acc)<<endl;
//         out_avg<<setw(9) << 1./a -1. <<"\t"<< setw(15) <<setprecision(15) << avg_pi<<"\t"<< setw(15)<<setprecision(15) << avg_zeta<<"\t"<< setw(9) << gsl_spline_eval(p_smg_spline, a, acc)/gsl_spline_eval(rho_smg_spline, a, acc) <<"\t"<< setw(9) << avg_phi<<"\t"<< setw(15) <<setprecision(15) <<gsl_spline_eval(cs2_spline, a, acc) <<"\t"<< setw(9) << tau<<endl;
//       }
//       else
//       {
//         out_avg<<setw(9) << 1./a -1. <<"\t"<< setw(15) <<setprecision(15) << avg_pi<<"\t"<< setw(15)<<setprecision(15) << avg_zeta<<"\t"<< setw(9) << avg_det_gamma<<"\t"<< setw(9) << avg_phi<<"\t"<< setw(15) <<setprecision(15) <<avg_cs2_full <<"\t"<< setw(9) << tau<<endl;
//       }
//         out_max<<setw(9) << tau <<"\t"<< setw(9) << max_pi<<"\t"<< setw(9) << max_zeta<<"\t"<< setw(9) << max_phi<<endl;
//
//         for (x.first(); x.test(); x.next())
//     	{
//           //NL_test, Printing out average
//         if(x.coord(0)==32 && x.coord(1)==20 && x.coord(2)==10)
//         {
//           // if(parallel.isRoot())
//           // {
//           out_real<<setw(9) << tau <<"\t"<< setw(9) <<pi (x)<<"\t"<< setw(9)<<pi_prime (x)<<"\t"<<x<<endl;
//           // }
//         }
//     	}
//       //****************************
//       //FOURIER PRINTING
//       //****************************
//       for(kFT.first();kFT.test();kFT.next())
//       {
//         norm_kFT_squared= kFT.coord(0)*kFT.coord(0) + kFT.coord(1) * kFT.coord(1) + kFT.coord(2) * kFT.coord(2);
//         if(norm_kFT_squared == 1)
//         {
//           out_fourier<<setw(9) << tau <<"\t"<< setw(9) << scalarFT_pi(kFT)<<"\t"<< setw(9)<<scalarFT_pi_prime (kFT)<<"\t"<<kFT<<"\t"<<norm_kFT_squared<<endl;
//         }
//       }
//       //**********************
//       //END ADDED************
//       //**********************
// #endif

#ifdef BENCHMARK
		cycle_start_time = MPI_Wtime();
#endif
		// construct stress-energy tensor
		projection_init(&source);
// #ifdef HAVE_CLASS
// 		if (sim.radiation_flag > 0 || sim.fluid_flag > 0)
// 			projection_T00_project(class_background, class_perturbs, source, scalarFT, &plan_source, sim, ic, cosmo, fourpiG, a);
// #endif
		if (sim.gr_flag > 0)
		{
			projection_T00_project(&pcls_cdm, &source, a, &phi);

			if (sim.baryon_flag)
				projection_T00_project(&pcls_b, &source, a, &phi);
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (a >= 1. / (sim.z_switch_deltancdm[i] + 1.) && sim.numpcl[1+sim.baryon_flag+i] > 0)
					projection_T00_project(pcls_ncdm+i, &source, a, &phi);
				else if (sim.radiation_flag == 0 || (a >= 1. / (sim.z_switch_deltancdm[i] + 1.) && sim.numpcl[1+sim.baryon_flag+i] == 0))
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
				if (a >= 1. / (sim.z_switch_deltancdm[i] + 1.) && sim.numpcl[1+sim.baryon_flag+i] > 0)
					scalarProjectionCIC_project(pcls_ncdm+i, &source);
			}
		}
		projection_T00_comm(&source);

    #ifdef VELOCITY
    		if ((sim.out_pk & MASK_VEL) || (sim.out_snapshot & MASK_VEL))
    		{
    			projection_init(&Bi);
                projection_Ti0_project(&pcls_cdm, &Bi, &phi, &chi);
                vertexProjectionCIC_comm(&Bi);
    						compute_vi_rescaled(cosmo, &vi, &source, &Bi, a, a_old
    							#ifdef HAVE_CLASS_BG
    							, H_spline, acc
    							#endif
    						);
                a_old = a;
    		}
    #endif

		if (sim.vector_flag == VECTOR_ELLIPTIC)
		{
			projection_init(&Bi);
			projection_T0i_project(&pcls_cdm, &Bi, &phi);
			if (sim.baryon_flag)
				projection_T0i_project(&pcls_b, &Bi, &phi);
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (a >= 1. / (sim.z_switch_Bncdm[i] + 1.) && sim.numpcl[1+sim.baryon_flag+i] > 0)
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
			{
				if (sim.numpcl[1+sim.baryon_flag+i] > 0)
					projection_Tij_project(pcls_ncdm+i, &Sij, a, &phi);
			}
		}
		projection_Tij_comm(&Sij);

#ifdef BENCHMARK
		projection_time += MPI_Wtime() - cycle_start_time;
		ref_time = MPI_Wtime();
#endif

// NO k-essence projection
// if (sim.Kess_source_gravity==1)
// {
// // Kessence projection Tmunu
// // In the projection zeta_integer comes, since synched with particles..
//  	if (sim.vector_flag == VECTOR_ELLIPTIC)
// 		{
// 			projection_Tmunu_kessence( T00_Kess,T0i_Kess,Tij_Kess, dx, a, phi, phi_old, 	chi, pi, pi_prime,
// 				#ifdef HAVE_CLASS_BG
// 				gsl_spline_eval(rho_smg_spline, a, acc)/gsl_spline_eval(rho_crit_spline, a, acc),
// 				gsl_spline_eval(p_smg_spline, a, acc)/gsl_spline_eval(rho_smg_spline, a, acc),
// 				gsl_spline_eval(cs2_spline, a, acc),
// 				Hconf(a, fourpiG, H_spline, acc)
// 				#else
// 				cosmo.Omega_kessence,
// 				cosmo.w_kessence,
// 				cosmo.cs2_kessence,
// 				Hconf(a, fourpiG, cosmo)
// 				#endif
// 				, fourpiG, sim.NL_kessence ,1 );
// 		}
//  	else
// 		{
// 			projection_Tmunu_kessence( T00_Kess,T0i_Kess,Tij_Kess, dx, a, phi, phi_old, 	chi, pi, pi_prime,
// 				#ifdef HAVE_CLASS_BG
// 				gsl_spline_eval(rho_smg_spline, a, acc)/gsl_spline_eval(rho_crit_spline, a, acc),
// 				gsl_spline_eval(p_smg_spline, a, acc)/gsl_spline_eval(rho_smg_spline, a, acc),
// 				gsl_spline_eval(cs2_spline, a, acc),
// 				Hconf(a, fourpiG, H_spline, acc)
// 				#else
// 				cosmo.Omega_kessence,
// 				cosmo.w_kessence,
// 				cosmo.cs2_kessence,
// 				Hconf(a, fourpiG, cosmo)
// 				#endif
// 				, fourpiG, sim.NL_kessence, 0 );
// 		}
//
// 		for (x.first(); x.test(); x.next())
// 		{
//        // if(x.coord(0)==32 && x.coord(1)==12 && x.coord(2)==32)cout<<"T00: "<< T00_Kess(x)<<" Phi:"<<phi(x)<<endl;
// 			// The coefficient is because it wanted to to be source according to eq C.2 of Gevolution paper
// 			// Note that it is multiplied to dx^2 and is divived by -a^3 because of definition of T00 which is scaled by a^3
// 			// We have T00 and Tij according to code's units, but source is important to calculate potentials and moving particles.
// 			// There is coefficient between Tij and Sij as source.
// 			source(x) += T00_Kess(x);
// 			if (sim.vector_flag == VECTOR_ELLIPTIC)for(int 	c=0;c<3;c++)Bi(x,c)+=  T0i_Kess(x,c);
// 			for(int c=0;c<6;c++)Sij(x,c)+=(2.) * Tij_Kess(x,c);
//       // if(x.coord(0)==32 && x.coord(1)==20 && x.coord(2)==10)
//       // {
//       // cout<<"x"<<x<<"T00_Kess(x): "<<T00_Kess(x)<<endl;
//       // }
// 		}
// }
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
			parallel.sum<double>(T00hom);
			T00hom /= (double) numpts3d;

			if (cycle % CYCLE_INFO_INTERVAL == 0)
			{
				COUT << " cycle " << cycle << ", background information: z = " << (1./a) - 1. << ", average T00 = " << T00hom << ", background model = " << cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo) << endl;
			}

			if (dtau_old > 0.)
			{
        Hc = Hconf(a, fourpiG,//TODO_EB
        #ifdef HAVE_CLASS_BG
          H_spline, acc
        #else
          cosmo
        #endif
          );

        // #ifdef BACKREACTION_TEST
        //
        // prepareFTsource_BackReactionTest<Real>(short_wave, dx, phi, chi, source, cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo), source, 3. * Hc * dx * dx / dtau_old, fourpiG * dx * dx / a, 3. * Hc * Hc * dx * dx, sim.boxsize);  // prepare nonlinear source for phi update
        // #else
        prepareFTsource<Real>(phi, chi, source, cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo), source, 3. * Hc * dx * dx / dtau_old, fourpiG * dx * dx / a, 3. * Hc * Hc * dx * dx);  // prepare nonlinear source for phi update
        // #endif

#ifdef BENCHMARK
				ref2_time= MPI_Wtime();
#endif
				plan_source.execute(FFT_FORWARD);  // go to k-space
#ifdef BENCHMARK
				fft_time += MPI_Wtime() - ref2_time;
				fft_count++;
#endif

    solveModifiedPoissonFT(scalarFT, scalarFT, 1. / (dx * dx), 3. * Hc / dtau_old);  // phi update (k-space)


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
      if (cycle == 0)
        outfile = fopen(filename, "w");
			outfile = fopen(filename, "a");
			if (outfile == NULL)
			{
				cout << " error opening file for background output!" << endl;
			}
			else
			{
        if (cycle == 0)
          fprintf(outfile, "# background statistics\n# cycle   tau/boxsize    a             conformal H/H0  phi(k=0)       T00(k=0)\n");
        fprintf(outfile, " %6d   %e   %e   %e   %e   %e\n", cycle, tau, a, Hconf(a, fourpiG,//TODO_EB
        #ifdef HAVE_CLASS_BG
          H_spline, acc
        #else
          cosmo
        #endif
        ) /
        Hconf(1., fourpiG,//TODO_EB
        #ifdef HAVE_CLASS_BG
          H_spline, acc
        #else
          cosmo
        #endif
        ), scalarFT(kFT).real(), T00hom);
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

// #ifdef HAVE_CLASS
// 		if (sim.radiation_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.))
// 		{
// 			prepareFTchiLinear(class_background, class_perturbs, scalarFT, sim, ic, cosmo, fourpiG, a);
// 			projectFTscalar(SijFT, scalarFT, 1);
// 		}
// 		else
// #endif
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


// lightcone output
if (sim.num_lightcone > 0)
  writeLightcones(sim, cosmo, fourpiG, a, tau, dtau, dtau_old, maxvel[0], cycle, h5filename + sim.basename_lightcone,
		#ifdef HAVE_CLASS_BG
		class_background, H_spline, acc,
		#endif
		&pcls_cdm, &pcls_b, pcls_ncdm, &phi, &chi, &Bi, &Sij, &BiFT, &SijFT, &plan_Bi, &plan_Sij, done_hij, IDbacklog);
else done_hij = 0;

#ifdef BENCHMARK
lightcone_output_time += MPI_Wtime() - ref_time;
ref_time = MPI_Wtime();
#endif

		// snapshot output
		if (snapcount < sim.num_snapshot && 1. / a < sim.z_snapshot[snapcount] + 1.)
		{
			COUT << COLORTEXT_CYAN << " writing snapshot" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;

      writeSnapshots(sim, cosmo, fourpiG, a, dtau_old, done_hij, snapcount, h5filename + sim.basename_snapshot,
        #ifdef HAVE_CLASS_BG
        H_spline, acc,
        #endif
        &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &pi,&pi_prime, &chi, &Bi, &source, &Bi, &Sij, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij
#ifdef CHECK_B
				, &Bi_check, &BiFT_check, &plan_Bi_check
#endif
#ifdef VELOCITY
				, &vi
#endif
			);

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

			writeSpectra(sim, cosmo, fourpiG, a, pkcount,
#ifdef HAVE_CLASS
				class_background, class_perturbs, ic,
#endif
				&pcls_cdm, &pcls_b, pcls_ncdm, &phi, &pi,&pi_prime, &chi, &Bi, &source, &Sij, &Sij ,&source, &Sij, &scalarFT, &scalarFT_pi, &scalarFT_pi_prime, &BiFT, &scalarFT, &BiFT, &SijFT, &SijFT, &plan_phi, &plan_pi , &plan_pi_prime, &plan_chi, &plan_Bi, &plan_source, &plan_Bi, &plan_Sij, &plan_source, &plan_Sij
#ifdef CHECK_B
				, &Bi_check, &BiFT_check, &plan_Bi_check
#endif
#ifdef VELOCITY
				, &vi, &viFT, &plan_vi
#endif
			);

			pkcount++;
		}

    // cout<<"EXACT_OUTPUT_REDSHIFTS: "<<EXACT_OUTPUT_REDSHIFTS<<endl;
    #ifdef EXACT_OUTPUT_REDSHIFTS
    		tmp = a;
				rungekutta4bg(tmp, fourpiG,
					#ifdef HAVE_CLASS_BG
						H_spline, acc,
					#else
						cosmo,
					#endif
					0.5 * dtau);
				rungekutta4bg(tmp, fourpiG,
					#ifdef HAVE_CLASS_BG
						H_spline, acc,
					#else
						cosmo,
					#endif
					0.5 * dtau);

    		if (pkcount < sim.num_pk && 1. / tmp < sim.z_pk[pkcount] + 1.)
    		{
    			writeSpectra(sim, cosmo, fourpiG, a, pkcount,
    #ifdef HAVE_CLASS
    					class_background, class_perturbs, ic,
    #endif
    					&pcls_cdm, &pcls_b, pcls_ncdm, &phi,&pi, &pi_prime, &chi, &Bi,&source, &Bi, &Sij, &source, &Sij, &scalarFT ,&scalarFT_pi, &scalarFT_pi_prime, &BiFT, &scalarFT, &BiFT, &SijFT, &SijFT, &plan_phi, &plan_pi, &plan_pi_prime, &plan_chi, &plan_Bi, &plan_source, &plan_Bi, &plan_Sij, &plan_source, &plan_Sij
    #ifdef CHECK_B
    					, &Bi_check, &BiFT_check, &plan_Bi_check
    #endif
    #ifdef VELOCITY
    				, &vi, &viFT, &plan_vi
    #endif
		    );
    		}
    #endif // EXACT_OUTPUT_REDSHIFTS


#ifdef BENCHMARK
		spectra_output_time += MPI_Wtime() - ref_time;
#endif

		if (pkcount >= sim.num_pk && snapcount >= sim.num_snapshot)
		{
			for (i = 0; i < sim.num_lightcone; i++)
			{
				if (sim.lightcone[i].z + 1. < 1. / a)
					i = sim.num_lightcone + 1;
			}
			if (i == sim.num_lightcone) break; // simulation complete
		}

		// compute number of step subdivisions for ncdm particle updates
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			if (dtau * maxvel[i+1+sim.baryon_flag] > dx * sim.movelimit)
				numsteps_ncdm[i] = (int) ceil(dtau * maxvel[i+1+sim.baryon_flag] / dx / sim.movelimit);
			else numsteps_ncdm[i] = 1;
		}

		if (cycle % CYCLE_INFO_INTERVAL == 0)
		{
			COUT << " cycle " << cycle << ", time integration information: max |v| = " << maxvel[0] << " (cdm Courant factor = " << maxvel[0] * dtau / dx;
			if (sim.baryon_flag)
			{
				COUT << "), baryon max |v| = " << maxvel[1] << " (Courant factor = " << maxvel[1] * dtau / dx;
			}

      COUT << "), time step / Hubble time = " << Hconf(a, fourpiG,//TODO_EB
			#ifdef HAVE_CLASS_BG
				H_spline, acc
			#else
				cosmo
			#endif
			) * dtau;


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

        // We just need to update halo when we want to calculate spatial derivative or use some neibours at the same time! So here wo do not nee to update halo for phi_prime!

 //Then fwe start the main loop zeta is updated to get zeta(n+1/2) from pi(n) and zeta(n-1/2)
 // EFT of k-essence theory
 if (cosmo.MGtheory == 0)
 {

   if(cycle==0)
   {
   if(parallel.isRoot())  cout << "\033[1;34mThe EFT approach equations are being solved!\033[0m\n";
   }
   //**********************
   //Kessence - LeapFrog:START
   //**********************
     double a_kess=a;
     //First we update pi_prime to have it at -1/2 just in the first loop
     if(cycle==0)
     {
       for (i=0;i<sim.nKe_numsteps;i++)
       {
         //computing pi_prime(-1/2) and zeta_int(-1) but we do not work with zeta(-1)
         update_zeta_eft(-dtau/ (2. * sim.nKe_numsteps) , dx, a_kess, phi, phi_old, chi, chi_old, pi, pi_prime,
   				#ifdef HAVE_CLASS_BG
   				gsl_spline_eval(rho_smg_spline, a_kess, acc)/gsl_spline_eval(rho_crit_spline, a_kess, acc),
   				gsl_spline_eval(p_smg_spline, a_kess, acc)/gsl_spline_eval(rho_smg_spline, a_kess, acc),
   				gsl_spline_eval(cs2_spline, a_kess, acc),
          0., // TODO: FH  "s" should be provided from hiclass columns
   				Hconf(a_kess, fourpiG, H_spline, acc),
   				Hconf_prime(a_kess, fourpiG, H_spline, acc)
   				#else
   				cosmo.Omega_kessence,
   				cosmo.w_kessence,
   				cosmo.cs2_kessence,
          ,0., // s parameter (FH:todo) should be read from hiclass as well as other params
   				Hconf(a_kess, fourpiG, cosmo),
   				Hconf_prime(a_kess, fourpiG, cosmo)
   				#endif
   				, sim.NL_kessence);
         // zeta_integer.updateHalo();
         pi_prime.updateHalo();
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
       update_zeta_eft(dtau/ sim.nKe_numsteps, dx, a_kess, phi, phi_old, chi, chi_old, pi, pi_prime,
   			#ifdef HAVE_CLASS_BG
   			gsl_spline_eval(rho_smg_spline, a_kess, acc)/gsl_spline_eval(rho_crit_spline, a_kess, acc),
   			gsl_spline_eval(p_smg_spline, a_kess, acc)/gsl_spline_eval(rho_smg_spline, a_kess, acc),
   			gsl_spline_eval(cs2_spline, a_kess, acc),
          0., // TODO: FH  "s" should be provided from hiclass columns
   			Hconf(a_kess, fourpiG, H_spline, acc),
   			Hconf_prime(a_kess, fourpiG, H_spline, acc)
   			#else
   			cosmo.Omega_kessence,
   			cosmo.w_kessence,
   			cosmo.cs2_kessence,
        ,0.0, // s parameter (FH:todo) should be read from hiclass as well as other params
   			Hconf(a_kess, fourpiG, cosmo),
   			Hconf_prime(a_kess, fourpiG, cosmo)
   			#endif
   			, sim.NL_kessence);
       // zeta_integer.updateHalo();
       pi_prime.updateHalo();
       //********************************************************************************
       //Since we have pi(n+1)=pi(n) + pi'(n+1/2), and in pi'(n+1/2) we have H(n+1/2) we update the background before updating the pi to have H(n+1/2), Moreover zeta(n+1) = zeta(n+1/2) + zeta'(n+1/2), so we put zeta_int updating in the pi updating!
       //********************************************************************************
   		rungekutta4bg(a_kess, fourpiG,
   			#ifdef HAVE_CLASS_BG
   				H_spline, acc,
   			#else
   				cosmo,
   			#endif
   			dtau  / sim.nKe_numsteps / 2.0);
       //********************************************************************************
       //we update pi to have it at n+1 (at first loop from the value at (0) and the value of zeta_integer at 1/2 and H(n+1/2) we update pi at (1))
       //In the pi update we also update zeta_int because we need the values of a_kess and H_kess at step n+1/2
       //By the below update we get pi(n+1) and zeta(n+1)
       //********************************************************************************
       update_pi_eft(dtau/ sim.nKe_numsteps, dx, a_kess, phi, phi_old, chi, chi_old, pi, pi_prime,
   			#ifdef HAVE_CLASS_BG
   			gsl_spline_eval(rho_smg_spline, a_kess, acc)/gsl_spline_eval(rho_crit_spline, a_kess, acc),
   			gsl_spline_eval(p_smg_spline, a_kess, acc)/gsl_spline_eval(rho_smg_spline, a_kess, acc),
   			gsl_spline_eval(cs2_spline, a_kess, acc),
   			Hconf(a_kess, fourpiG, H_spline, acc),
   			Hconf_prime(a_kess, fourpiG, H_spline, acc)
   			#else
   			cosmo.Omega_kessence,
   			cosmo.w_kessence,
   			cosmo.cs2_kessence,
   			Hconf(a_kess, fourpiG, cosmo),
   			Hconf_prime(a_kess, fourpiG, cosmo)
   			#endif
   			, sim.NL_kessence); // H_old is updated here in the function
   		pi.updateHalo();

       //********************************************************************************
       // Now we have pi(n+1) and a_kess(n+1/2) so we update background by halfstep to have a_kess(n+1)
       //********************************************************************************
   		rungekutta4bg(a_kess, fourpiG,
   			#ifdef HAVE_CLASS_BG
   				H_spline, acc,
   			#else
   				cosmo,
   			#endif
   			dtau  / sim.nKe_numsteps / 2.0);

        //*************************
        #ifdef BACKREACTION_TEST
          avg_zeta =average(  pi_prime,1., numpts3d ) ;
          // avg_zeta_old =average(  pi_prime_old,1., numpts3d ) ;
          avg_pi =average(  pi,1., numpts3d ) ;
          avg_phi =average(  phi , 1., numpts3d ) ;
          avg_pi_old =average(  pi_old, 1., numpts3d ) ;
         if (avg_zeta > 4.e-7 && abs(avg_pi/avg_pi_old)>1.015 && snapcount_b< sim.num_snapshot_kess )
         {
         if(parallel.isRoot())  cout << "\033[1;32mThe blowup criteria for EFT equations are met, the requested snapshots being produced\033[0m\n";
           writeSpectra(sim, cosmo, fourpiG, a, snapcount_b,
             #ifdef HAVE_CLASS
             				class_background, class_perturbs, ic,
             #endif
                     &pcls_cdm, &pcls_b, pcls_ncdm, &phi,&pi, &pi_prime, &chi, &Bi,&T00_Kess, &T0i_Kess, &Tij_Kess, &source, &Sij, &scalarFT ,&scalarFT_pi, &scalarFT_pi_prime, &BiFT, &T00_KessFT, &T0i_KessFT, &Tij_KessFT, &SijFT, &plan_phi, &plan_pi, &plan_pi_prime, &plan_chi, &plan_Bi, &plan_T00_Kess, &plan_T0i_Kess, &plan_Tij_Kess, &plan_source, &plan_Sij);
             str_filename =  "./output/pi_" + to_string(snapcount_b) + ".h5";
             str_filename2 = "./output/zeta_" + to_string(snapcount_b) + ".h5";
             str_filename3 = "./output/phi_" + to_string(snapcount_b) + ".h5";
             pi.saveHDF5(str_filename);
             pi_prime.saveHDF5(str_filename2);
             phi.saveHDF5(str_filename3);
             // str_filename =  "./output/pi_" + to_string(snapcount_b-1) + ".h5";
             // str_filename2 = "./output/zeta_" + to_string(snapcount_b-1) + ".h5";
             // pi_old.saveHDF5(str_filename);
             // pi_prime_old.saveHDF5(str_filename2);
             snapcount_b++;

           //****************************
           //****PRINTING snapshots info
           //****************************
             // COUT << scientific << setprecision(8);
             // if(parallel.isRoot())
             // {
             // out_snapshots<<"### 1- tau\t2- z \t3- a\t 4- zeta_avg\t 5- avg_pi\t 6- avg_phi\t 7- tau/boxsize\t 8- H_conf/H0 \t 9- snap_count"<<endl;


             out_snapshots<<setw(9) << tau + dtau/sim.nKe_numsteps <<"\t"<< setw(9) << 1./(a_kess) -1.0 <<"\t"<< setw(9) << a_kess <<"\t"<< setw(9) << avg_zeta <<"\t"<< setw(9) << avg_pi <<"\t"<< setw(9) << avg_phi <<"\t"<< setw(9) <<tau <<"\t"<< setw(9) << Hconf(a_kess, fourpiG,//TODO_EB
   					#ifdef HAVE_CLASS_BG
   						H_spline, acc
   					#else
   						cosmo
   					#endif
   				) / Hconf(1., fourpiG,//TODO_EB
   					#ifdef HAVE_CLASS_BG
   						H_spline, acc
   					#else
   						cosmo
   					#endif
   					) <<"\t"<< setw(9) <<snapcount_b  <<endl;
           }
       #endif  //*****Backreaction test**********
       }

   #ifdef BENCHMARK
       kessence_update_time += MPI_Wtime() - ref_time;
       ref_time = MPI_Wtime();
   #endif
   //**********************
   //Kessence - LeapFrog: End
   //**********************
 }


    // Fundamental k-essence theory  cosmo.MGtheory == 1
  if (cosmo.MGtheory == 1)
  {

      if (cosmo.solver_kessence==0) // Euler method is being used!
        {
          if(cycle==0)
          {
          if(parallel.isRoot())  cout << "\033[1;34mThe fundamental theory is being solved using the Euler method!\033[0m\n";
          if(parallel.isRoot())  cout << "\033[1;34mThe Non-linear terms included: \033[0m"<<cosmo.NL<<endl;
          // if(parallel.isRoot())  cout << "\033[1;34mThe initial condition for the scalar field is made using hiclass: \033[0m"<<cosmo.IC_scf<<endl;
          // if(parallel.isRoot() & (cosmo.IC_scf==0) )  cout << "\033[1;34mThe amplitude of the initial condition is set as p_prime(x) = 0 and pi(x) = phi_bg(hiclass) + Amplitude * Phi(x) where the Amplitude is \033[0m"<<cosmo.IC_amplitude<<endl;

          }
        double a_kess=a;

      //****************************************************
      // Euler algorithm
      //****************************************************
      for (i=0;i<sim.nKe_numsteps;i++)
      {

        #ifdef BACKREACTION_TEST
              //****************************
              //****PRINTING AVERAGE OVER TIME
              //****************************
              avg_pi =average(  pi, Hconf(1.0, fourpiG,//TODO_EB
        			#ifdef HAVE_CLASS_BG
        				H_spline, acc
        			#else
        				cosmo
        			#endif
        				), numpts3d ) ;
              avg_zeta =average( pi_prime,1., numpts3d ) ;
              avg_phi =average(  phi , 1., numpts3d ) ;
              avg_det_gamma =average(  det_gamma , 1., numpts3d ) ;
              avg_cs2_full =average(  cs2_full , 1., numpts3d ) ;

              max_pi =maximum(  pi, Hconf(1.0, fourpiG,//TODO_EB
        			#ifdef HAVE_CLASS_BG
        				H_spline, acc
        			#else
        				cosmo
        			#endif
        				), numpts3d ) ;
              // max_pi_old =maximum(  pi_old, 1., numpts3d ) ;
              max_zeta =maximum(  pi_prime,1., numpts3d ) ;
              max_phi =maximum(  phi , 1., numpts3d ) ;
              max_cs2 = maximum(  cs2_full , 1., numpts3d ) ;
              //
              min_pi =minimum(  pi, Hconf(1.0, fourpiG,//TODO_EB
              #ifdef HAVE_CLASS_BG
                H_spline, acc
              #else
                cosmo
              #endif
                ), numpts3d ) ;
              min_zeta =minimum(  pi_prime, 1., numpts3d ) ;
              min_phi =minimum(  phi , 1., numpts3d ) ;
              min_cs2 = minimum(  cs2_full , 1., numpts3d ) ;

              COUT << scientific << setprecision(8);
              if(isnan(avg_zeta) & isnan(avg_pi))
              {
                if(parallel.isRoot()) cout << "\033[1;34mThe PDE blowus at z=: \033[0m"<<1./a_kess - 1<<endl;
                parallel.abortForce();
              }
                if (cosmo.MGtheory == 0)
              {
                cout<<"z = "<<1./a_kess - 1<<" p_smg_spline:"<<gsl_spline_eval(p_smg_spline, a_kess, acc)<<" rho_smg_spline:"<<gsl_spline_eval(rho_smg_spline, a_kess, acc)<<endl;
                out_avg<<setw(9) << 1./a_kess -1. <<"\t"<< setw(15) <<setprecision(15) << avg_pi<<"\t"<< setw(15)<<setprecision(15) << avg_zeta<<"\t"<< setw(9) << gsl_spline_eval(p_smg_spline, a_kess, acc)/gsl_spline_eval(rho_smg_spline, a_kess, acc) <<"\t"<< setw(9) << avg_phi<<"\t"<< setw(15) <<setprecision(15) <<gsl_spline_eval(cs2_spline, a_kess, acc) <<"\t"<< setw(9) << tau<<endl;
              }
              else
              {
                out_avg<<setw(9) << 1./a_kess -1. <<"\t"<< setw(9) << avg_pi<<"\t"<< setw(9) << avg_zeta<<"\t"<< setw(9) << avg_det_gamma<<"\t"<< setw(9) << avg_phi<<"\t"<< setw(9) <<avg_cs2_full <<"\t"<< setw(9) << tau<<endl;
              }
                out_max<<setw(9) << 1./a_kess - 1 <<"\t"<< setw(9) << max_pi<<"\t"<< setw(9) << max_zeta<<"\t"<< setw(9) << max_phi<<"\t"<< setw(9) << max_cs2 <<endl;

                out_min<< setw(9) << 1./a_kess - 1 <<"\t"<< setw(9) << min_pi<<"\t"<< setw(9) << min_zeta<<"\t"<< setw(9) << min_phi<<"\t"<< setw(9) << min_cs2 <<endl;
        #endif
      update_pi_prime_euler(dtau/ sim.nKe_numsteps, dx, a_kess, pi, deltaX, pi_prime, det_gamma, cs2_full, cosmo.X_hat, cosmo.g0, cosmo.g2, cosmo.g4,
       #ifdef HAVE_CLASS_BG
       Hconf(a_kess, fourpiG, H_spline, acc)
       #else
       Hconf(a_kess, fourpiG, cosmo)
       #endif
        , cosmo.NL);
      deltaX.updateHalo();
      pi_prime.updateHalo();
      pi.updateHalo();
      // Although it's an Euler algorithm we can update the BG part twice in each step to increase the precision!
      rungekutta4bg(a_kess, fourpiG,
  			#ifdef HAVE_CLASS_BG
  				H_spline, acc,
  			#else
  				cosmo,
  			#endif
  			dtau  / sim.nKe_numsteps );// We don't use the BG in the update pi, so we need to update it once!
    }
  }
    //Euler method end!


    else if (cosmo.solver_kessence==1) // Leap-frog method is being used!
      {
        double a_kess=a;

        if(cycle==0)
        {
          if(parallel.isRoot())  cout << "\033[1;34mThe fundamental theory is being solved using leap-frog method!\033[0m\n";

          // for (i=0;i<sim.nKe_numsteps;i++)
          // {
            //computing pi_prime(-1/2) and zeta_int(-1) but we do not work with zeta(-1)
            update_pi_prime_leap_frog(-dtau/ (2.0 * sim.nKe_numsteps), dx, a_kess, pi, pi_prime, pi_prime_dot, det_gamma, cs2_full, cosmo.X_hat, cosmo.g0, cosmo.g2, cosmo.g4,
             #ifdef HAVE_CLASS_BG
             Hconf(a_kess, fourpiG, H_spline, acc)
             #else
             Hconf(a_kess, fourpiG, cosmo)
             #endif
              );            // zeta_integer.updateHalo();
            pi_prime.updateHalo(); // phi' (n-1/2)
        }

      //****************************************************
      // Leap-frod algorithm
      //****************************************************
      for (i=0;i<sim.nKe_numsteps;i++)
      {

        update_pi_prime_leap_frog(dtau/ sim.nKe_numsteps, dx, a_kess,pi, pi_prime, pi_prime_dot, det_gamma, cs2_full, cosmo.X_hat, cosmo.g0, cosmo.g2, cosmo.g4,
         #ifdef HAVE_CLASS_BG
         Hconf(a_kess, fourpiG, H_spline, acc)
         #else
         Hconf(a_kess, fourpiG, cosmo)
         #endif
       );  // phi' (n-1/2) [phi(n), phi(n)] -->  phi' (n+1/2)
        pi_prime.updateHalo();
        rungekutta4bg(a_kess, fourpiG,
          #ifdef HAVE_CLASS_BG
            H_spline, acc,
          #else
            cosmo,
          #endif
          dtau  / sim.nKe_numsteps / 2.0);

        update_pi_full(dtau/ sim.nKe_numsteps, pi, pi_prime); // H_old is updated here in the function
        pi.updateHalo();
        // Although it's an Euler algorithm we can update the BG part twice in each step to increase the precision!
        rungekutta4bg(a_kess, fourpiG,
          #ifdef HAVE_CLASS_BG
            H_spline, acc,
          #else
            cosmo,
          #endif
          dtau  / sim.nKe_numsteps / 2.0);
      }
    // End of Leap-Frog algorithm
  }


  else if (cosmo.solver_kessence==2) // RK2 method is being used!
    {
      double a_kess=a;

      if(cycle==0)
      {
        if(parallel.isRoot()) cout << "\033[1;34mThe k-essence theory is being solved using RK2 method!\033[0m\n";
      }
    //****************************************************
    // RK2 algorithm - We don't update BG in the middle steps! Todo!
    //****************************************************
    for (i=0;i<sim.nKe_numsteps;i++)
    {
      update_RK2(dtau/ sim.nKe_numsteps, dx, a_kess,pi, pi_prime, k1, k2, l1, l2, cs2_full, det_gamma, cosmo.X_hat, cosmo.g0, cosmo.g2, cosmo.g4,
       #ifdef HAVE_CLASS_BG
       Hconf(a_kess, fourpiG, H_spline, acc)
       #else
       Hconf(a_kess, fourpiG, cosmo)
       #endif
     , cosmo.NL);
      pi.updateHalo();
      pi_prime.updateHalo();

      rungekutta4bg(a_kess, fourpiG,
        #ifdef HAVE_CLASS_BG
          H_spline, acc,
        #else
          cosmo,
        #endif
        dtau  / sim.nKe_numsteps / 1.0);
    }
  // End of RK2 algorithm
}


else if (cosmo.solver_kessence==3) // RK4 method is being used!
  {
    double a_kess=a;

    if(cycle==0)
    {
      if(parallel.isRoot())  cout << "\033[1;34mThe fundamental theory is being solved using RK4 method!\033[0m\n";

    }
}

  //****************************************************
  // RK4 algorithm
  //****************************************************
//   for (i=0;i<sim.nKe_numsteps;i++)
//   {
//   //   update_pi_RK2(dtau/ sim.nKe_numsteps, pi, pi_prime, pi_prime_dot);  // phi' (n-1/2) [phi(n), phi(n)] -->  phi' (n+1/2)
//   //     pi.updateHalo();
//   //
//   //  rungekutta4bg(a_kess, fourpiG,
//   //    #ifdef HAVE_CLASS_BG
//   //      H_spline, acc,
//   //    #else
//   //      cosmo,
//   //    #endif
//   //    dtau  / sim.nKe_numsteps / 2.0);
//   //
//   //   update_pi_prime_RK2(dtau/ sim.nKe_numsteps, dx, a_kess,pi, pi_prime, pi_prime_dot, det_gamma, cs2_full, cosmo.X_hat, cosmo.g0, cosmo.g2, cosmo.g4,
//   //    #ifdef HAVE_CLASS_BG
//   //    Hconf(a_kess, fourpiG, H_spline, acc)
//   //    #else
//   //    Hconf(a_kess, fourpiG, cosmo)
//   //    #endif
//   //  );  // phi' (n-1/2) [phi(n), phi(n)] -->  phi' (n+1/2)
//   //   pi_prime.updateHalo();
//   //   rungekutta4bg(a_kess, fourpiG,
//   //     #ifdef HAVE_CLASS_BG
//   //       H_spline, acc,
//   //     #else
//   //       cosmo,
//   //     #endif
//   //     dtau  / sim.nKe_numsteps / 2.0);
//   // }
// // End of RK4
// }
   // Fundamental k-essence theory  END
}

#ifdef BENCHMARK
    kessence_update_time += MPI_Wtime() - ref_time;
    ref_time = MPI_Wtime();
#endif


    //
		// for (j = 0; j < numsteps; j++) // particle update
		// {
#ifdef BENCHMARK
		ref2_time = MPI_Wtime();
#endif
		for (i = 0; i < cosmo.num_ncdm; i++) // non-cold DM particle update
		{
			if (sim.numpcl[1+sim.baryon_flag+i] == 0) continue;

			tmp = a;

			for (j = 0; j < numsteps_ncdm[i]; j++)
			{
				f_params[0] = tmp;
				f_params[1] = tmp * tmp * sim.numpts;
				if (sim.gr_flag > 0)
					maxvel[i+1+sim.baryon_flag] = pcls_ncdm[i].updateVel(update_q, (dtau + dtau_old) / 2. / numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
				else
					maxvel[i+1+sim.baryon_flag] = pcls_ncdm[i].updateVel(update_q_Newton, (dtau + dtau_old) / 2. / numsteps_ncdm[i], update_ncdm_fields, ((sim.radiation_flag + sim.fluid_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);

#ifdef BENCHMARK
				update_q_count++;
				update_q_time += MPI_Wtime() - ref2_time;
				ref2_time = MPI_Wtime();
#endif

        rungekutta4bg(tmp, fourpiG,
          #ifdef HAVE_CLASS_BG
            H_spline, acc,
          #else
            cosmo,
          #endif
          0.5 * dtau / numsteps_ncdm[i]);
        f_params[0] = tmp;
        f_params[1] = tmp * tmp * sim.numpts;

				if (sim.gr_flag > 0)
					pcls_ncdm[i].moveParticles(update_pos, dtau / numsteps_ncdm[i], update_ncdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
				else
					pcls_ncdm[i].moveParticles(update_pos_Newton, dtau / numsteps_ncdm[i], NULL, 0, f_params);
#ifdef BENCHMARK
				moveParts_count++;
				moveParts_time += MPI_Wtime() - ref2_time;
				ref2_time = MPI_Wtime();
#endif
      rungekutta4bg(tmp, fourpiG,
        #ifdef HAVE_CLASS_BG
          H_spline, acc,
        #else
          cosmo,
        #endif
        0.5 * dtau / numsteps_ncdm[i]);
			}
		}

		// cdm and baryon particle update
		f_params[0] = a;
		f_params[1] = a * a * sim.numpts;
		if (sim.gr_flag > 0)
		{
			maxvel[0] = pcls_cdm.updateVel(update_q, (dtau + dtau_old) / 2., update_cdm_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
			if (sim.baryon_flag)
				maxvel[1] = pcls_b.updateVel(update_q, (dtau + dtau_old) / 2., update_b_fields, (1. / a < ic.z_relax + 1. ? 3 : 2), f_params);
		}
		else
		{
			maxvel[0] = pcls_cdm.updateVel(update_q_Newton, (dtau + dtau_old) / 2., update_cdm_fields, ((sim.radiation_flag + sim.fluid_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);
			if (sim.baryon_flag)
				maxvel[1] = pcls_b.updateVel(update_q_Newton, (dtau + dtau_old) / 2., update_b_fields, ((sim.radiation_flag + sim.fluid_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.)) ? 2 : 1), f_params);
		}

#ifdef BENCHMARK
		update_q_count++;
		update_q_time += MPI_Wtime() - ref2_time;
		ref2_time = MPI_Wtime();
#endif

        rungekutta4bg(a, fourpiG,
        #ifdef HAVE_CLASS_BG
          H_spline, acc,
        #else
          cosmo,
        #endif
        0.5 * dtau);  // evolve background by half a time step


		f_params[0] = a;
		f_params[1] = a * a * sim.numpts;
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
#endif

  rungekutta4bg(a, fourpiG,
    #ifdef HAVE_CLASS_BG
      H_spline, acc,
    #else
      cosmo,
    #endif
    0.5 * dtau);  // evolve background by half a time step

		parallel.max<double>(maxvel, numspecies);

		if (sim.gr_flag > 0)
		{
			for (i = 0; i < numspecies; i++)
				maxvel[i] /= sqrt(maxvel[i] * maxvel[i] + 1.0);
		}
		// done particle update

		tau += dtau;

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
					hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, pi, pi_prime, chi, Bi_check, a, tau, dtau, cycle);
				}
				else
#endif
				hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, pi, pi_prime, chi, Bi, a, tau, dtau, cycle);
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
				hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, pi, pi_prime, chi, Bi, a, tau, dtau, cycle, restartcount);
			}
			else
#endif
			hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, pi, pi_prime, chi, Bi, a, tau, dtau, cycle, restartcount);
			restartcount++;
		}

    dtau_old = dtau;

		if (sim.Cf * dx < sim.steplimit / Hconf(a, fourpiG,//TODO_EB
		#ifdef HAVE_CLASS_BG
			H_spline, acc
		#else
			cosmo
		#endif
		))
			dtau = sim.Cf * dx;
		else
			dtau = sim.steplimit / Hconf(a, fourpiG,//TODO_EB
			#ifdef HAVE_CLASS_BG
				H_spline, acc
			#else
				cosmo
			#endif
			);

		cycle++;

#ifdef BENCHMARK
		cycle_time += MPI_Wtime()-cycle_start_time;
#endif
	}

	COUT << COLORTEXT_GREEN << " simulation complete." << COLORTEXT_RESET << endl;

#ifdef BENCHMARK
		ref_time = MPI_Wtime();
#endif

#ifdef HAVE_CLASS
	if (sim.radiation_flag > 0 || sim.fluid_flag > 0)
    freeCLASSstructures(class_background, class_thermo, class_perturbs);
#endif

#ifdef HAVE_CLASS_BG
// Free interpolation structures
gsl_interp_accel_free(acc);
gsl_spline_free(H_spline);
gsl_spline_free(p_smg_spline);
gsl_spline_free(cs2_spline);
gsl_spline_free(rho_smg_spline);
gsl_spline_free(rho_crit_spline);
COUT << endl << " hiclass interpolation structures correctly deallocated." << endl << endl;
#endif

#ifdef BENCHMARK
	lightcone_output_time += MPI_Wtime() - ref_time;
	run_time = MPI_Wtime() - start_time;

	parallel.sum(run_time);
	parallel.sum(cycle_time);
	parallel.sum(projection_time);
	parallel.sum(snapshot_output_time);
	parallel.sum(spectra_output_time);
	parallel.sum(lightcone_output_time);
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
	COUT << "lightcone outputs          : "<< hourMinSec(lightcone_output_time) << " ; " << 100. * lightcone_output_time/cycle_time <<"%."<<endl;
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
