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
// Last modified: April 2017
//
//////////////////////////

#include <stdlib.h>
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

using namespace std;

using namespace LATfield2;

int main(int argc, char **argv)
{

#ifdef BENCHMARK
	//benchmarking variables

	double ref_time, ref2_time, cycle_start_time;

	double initialization_time;
	double run_time;
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

#endif  //BENCHMARK

	int n = 0, m = 0;
	int io_size = 0;
	int io_group_size = 0;

	int i, j, cycle = 0, snapcount = 0, pkcount = 0, restartcount = 0, usedparams, numparam = 0, numsteps, numspecies, done_hij;
	int numsteps_ncdm[MAX_PCL_SPECIES-2];
	long numpts3d;
	int box[3];
	double dtau, dtau_old, dtau_older, dx, tau, a, fourpiG, tau_Lambda, tmp, start_time;
	double maxvel[MAX_PCL_SPECIES];
	FILE * outfile;
	char filename[2*PARAM_MAX_LENGTH+24];
	string h5filename;
	char * settingsfile = NULL;
	parameter * params = NULL;
	metadata sim;
	cosmology cosmo;
	icsettings ic;
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
	COUT << "  _   _      _         __ ,  _" << endl;
	COUT << " (_| (-' \\/ (_) (_ (_| (  ( (_) /\\/	version 1.2 beta    running on " << n*m << " cores." << endl;
	COUT << "  -'" << endl << COLORTEXT_RESET << endl;

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
	Field<Real> phi_old;
	//Kessence
	Field<Real> pi_k;
	Field<Real> pi_v_k;
	Field<Real> T00_Kess;
	Field<Real> T0i_Kess;
	Field<Real> Tij_Kess;
	Field<Real> source;
	Field<Real> chi;
	Field<Real> chi_old;
	Field<Real> Sij;
	Field<Real> Bi;
	Field<Cplx> scalarFT;
	Field<Cplx> scalarFT_phi_old;
	Field<Cplx> scalarFT_chi_old;
	Field<Cplx> scalarFT_pi;
	Field<Cplx> scalarFT_pi_v;
	Field<Cplx> SijFT;
	Field<Cplx> BiFT;
	source.initialize(lat,1);
	phi.initialize(lat,1);
	phi_old.initialize(lat,1);
	//Kessence part
	pi_k.initialize(lat,1);
	pi_v_k.initialize(lat,1);
	T00_Kess.initialize(lat,1);
	T00_Kess.alloc();
	T0i_Kess.initialize(lat,3);
	T0i_Kess.alloc();
	Tij_Kess.initialize(lat,3,3,symmetric);
	Tij_Kess.alloc();
	//
	chi.initialize(lat,1);
	chi_old.initialize(lat,1);
	scalarFT.initialize(latFT,1);
	scalarFT_phi_old.initialize(latFT,1);
	scalarFT_chi_old.initialize(latFT,1);
	//Kessence part
	scalarFT_pi.initialize(latFT,1);
	scalarFT_pi_v.initialize(latFT,1);
	//
	PlanFFT<Cplx> plan_source(&source, &scalarFT);
	PlanFFT<Cplx> plan_phi(&phi, &scalarFT);
	PlanFFT<Cplx> plan_chi(&chi, &scalarFT);
	PlanFFT<Cplx> plan_phi_old(&phi_old, &scalarFT_phi_old);
	PlanFFT<Cplx> plan_chi_old(&chi_old, &scalarFT_chi_old);
	//Kessence part
	PlanFFT<Cplx> plan_pi_k(&pi_k, &scalarFT_pi);
	PlanFFT<Cplx> plan_pi_v_k(&pi_v_k, &scalarFT_pi_v);
	//
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

	Field<Real> * lcbuffer[LIGHTCONE_MAX_FIELDS*LIGHTCONE_THICKNESS];
	Lattice * lclat[LIGHTCONE_MAX_FIELDS];
#ifdef LIGHTCONE_DOWNGRADE
	int dgbox[3];

	for (i = 0; i < 3; i++) dgbox[i] = box[i] / LIGHTCONE_DOWNGRADE;
#endif

	if (sim.num_lightcone > 0)
	{
		for (i = 0, j = 0; i < sim.num_lightcone; i++)
		{
			j |= sim.out_lightcone[i];
		}

		if (j & MASK_PHI)
		{
			lclat[LIGHTCONE_PHI_OFFSET] = new Lattice(3,box,0);
			for (i = 0; i < LIGHTCONE_THICKNESS; i++)
				lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_PHI_OFFSET+i] = new Field<Real>((*lclat[LIGHTCONE_PHI_OFFSET]), 1);
		}

		if (j & MASK_CHI)
		{
			lclat[LIGHTCONE_CHI_OFFSET] = new Lattice(3,box,0);
			for (i = 0; i < LIGHTCONE_THICKNESS; i++)
				lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_CHI_OFFSET+i] = new Field<Real>((*lclat[LIGHTCONE_CHI_OFFSET]), 1);
		}

		if (j & MASK_B)
		{
			lclat[LIGHTCONE_B_OFFSET] = new Lattice(3,box,0);
			for (i = 0; i < 3*LIGHTCONE_THICKNESS; i++)
				lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_B_OFFSET+i] = new Field<Real>((*lclat[LIGHTCONE_B_OFFSET]), 1);
		}

		if (j & MASK_HIJ)
		{
#ifdef LIGHTCONE_DOWNGRADE
			lclat[LIGHTCONE_HIJ_OFFSET] = new Lattice(3,dgbox,0);
#else
			lclat[LIGHTCONE_HIJ_OFFSET] = new Lattice(3,box,0);
#endif
			for (i = 0; i < 5*LIGHTCONE_THICKNESS; i++)
				lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+i] = new Field<Real>((*lclat[LIGHTCONE_HIJ_OFFSET]), 1);
		}
	}

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
	tau = particleHorizon(a, fourpiG, cosmo);
	tau_Lambda = -1.0;
	if (sim.Cf * dx < sim.steplimit / Hconf(a, fourpiG, cosmo))
		dtau = sim.Cf * dx / sim.nKe_numsteps;
	else
		dtau = sim.steplimit / Hconf(a, fourpiG, cosmo) / sim.nKe_numsteps;
	dtau_old = 0.;
	dtau_older = 0.;

	if (ic.generator == ICGEN_BASIC)
		generateIC_basic(sim, ic, cosmo, fourpiG, &pcls_cdm, &pcls_b, pcls_ncdm, maxvel, &phi, &pi_k, &pi_v_k, &chi, &Bi, &source, &Sij, &scalarFT, &scalarFT_pi,&scalarFT_pi_v, &BiFT, &SijFT, &plan_phi, &plan_pi_k,&plan_pi_v_k, &plan_chi, &plan_Bi, &plan_source, &plan_Sij); // generates ICs on the fly

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

#ifdef BENCHMARK
	initialization_time = MPI_Wtime() - start_time;
	parallel.sum(initialization_time);
	COUT << COLORTEXT_GREEN << " initialization complete." << COLORTEXT_RESET << " BENCHMARK: " << hourMinSec(initialization_time) << endl << endl;
#else
	COUT << COLORTEXT_GREEN << " initialization complete." << COLORTEXT_RESET << endl << endl;
#endif

#ifdef HAVE_CLASS
#ifdef COSIRA_HACK
	initializeCLASSstructures(sim, ic, cosmo, class_background, class_perturbs, class_spectra);
#endif // COSIRA_HACK
	if (sim.radiation_flag > 0)
	{
#ifndef COSIRA_HACK
		initializeCLASSstructures(sim, ic, cosmo, class_background, class_perturbs, class_spectra);
#endif // COSIRA_HACK
#ifndef MULTISTEP_PROJECTION
		if (sim.gr_flag > 0 && a < 1. / (sim.z_switch_linearchi + 1.) && (ic.generator == ICGEN_BASIC || (ic.generator == ICGEN_READ_FROM_DISK && cycle == 0)))
		{
			prepareFTchiLinear(class_background, class_perturbs, class_spectra, scalarFT, sim, ic, cosmo, fourpiG, a);
			plan_source.execute(FFT_BACKWARD);
			for (x.first(); x.test(); x.next())
				chi(x) += source(x);
			chi.updateHalo();
		}
#endif
	}
#endif
 // double distance;
 // for (x.first(); x.test(); x.next())
 // {
 //  distance = (x.coord(0)-lat.size(0)/2.)*(x.coord(0)-lat.size(0)/2.);
 //  distance +=(x.coord(1)-lat.size(1)/2.)*(x.coord(1)-lat.size(1)/2.);
 //  distance += (x.coord(2)-lat.size(2)/2.)*(x.coord(2)-lat.size(2)/2.);
 //  pi_k(x)=3*exp(-distance/(4.));
 //  pi_v_k(x)=-distance;
 // }
	while (true)    // main loop
	{
	for (x.first(); x.test(); x.next())
		{
			phi_old(x) =phi(x);
			chi_old(x) =chi(x);

		}
	phi_old.updateHalo();
	chi_old.updateHalo();

#ifdef BENCHMARK
		cycle_start_time = MPI_Wtime();
#endif
		// construct stress-energy tensor
#ifndef MULTISTEP_PROJECTION	// multistep projection: projection_init is called at the beginning of the particle update
		projection_init(&source);
#ifdef HAVE_CLASS
		if (sim.radiation_flag > 0)
			projection_T00_project(class_background, class_perturbs, class_spectra, source, scalarFT, &plan_source, sim, ic, cosmo, fourpiG, a);
#endif
#endif
		if (sim.gr_flag > 0)
		{
			projection_T00_project(&pcls_cdm, &source, a, &phi);
			if (sim.baryon_flag)
				projection_T00_project(&pcls_b, &source, a, &phi);
#ifndef MULTISTEP_PROJECTION	// multistep projection: non-cold species are projected at each of their update steps
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
#endif
		}
		else
		{
			scalarProjectionCIC_project(&pcls_cdm, &source);
			if (sim.baryon_flag)
				scalarProjectionCIC_project(&pcls_b, &source);
#ifndef MULTISTEP_PROJECTION	// multistep projection: non-cold species are projected at each of their update steps
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (a >= 1. / (sim.z_switch_deltancdm[i] + 1.))
					scalarProjectionCIC_project(pcls_ncdm+i, &source);
			}
#endif
		}
		projection_T00_comm(&source);

#ifndef VECTOREXTRA
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
#endif

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

 		// Kessence projection Tmunu
		 if (sim.vector_flag == VECTOR_ELLIPTIC)
		 	{
		 		projection_Tmunu_kessence( T00_Kess,T0i_Kess,Tij_Kess, dx, a, phi, phi_old, chi, pi_k, pi_v_k, cosmo.Omega_kessence, cosmo.w_kessence, cosmo.cs2_kessence, Hconf(a, fourpiG, cosmo), fourpiG, 1 );
		 	}
		 else
		 	{
		 		projection_Tmunu_kessence( T00_Kess,T0i_Kess,Tij_Kess, dx, a, phi, phi_old, chi, pi_k, pi_v_k, cosmo.Omega_kessence, cosmo.w_kessence, cosmo.cs2_kessence, Hconf(a, fourpiG, cosmo), fourpiG, 0 );
		 	}

			for (x.first(); x.test(); x.next())
			{
				source(x) += T00_Kess(x);
				if (sim.vector_flag == VECTOR_ELLIPTIC)for(int c=0;c<3;c++)Bi(x,c)+=T0i_Kess(x,c);
				for(int c=0;c<6;c++)Sij(x,c)+=Tij_Kess(x,c);
			}

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
					fprintf(outfile, "# background statistics\n# cycle   tau/boxsize    a             conformal H/H0  phi(k=0)       T00(k=0)\n");
				fprintf(outfile, " %6d   %e   %e   %e   %e   %e\n", cycle, tau, a, Hconf(a, fourpiG, cosmo) / Hconf(1., fourpiG, cosmo), scalarFT(kFT).real(), T00hom);
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
#ifdef VECTOREXTRA
			prepareFTvector(phi, Bi, Sij, 0.5 / (fourpiG * dx * dx));
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
			for (x.first(); x.test(); x.next())
			{
				Bi(x, 0) += Sij(x, 0);
				Bi(x, 1) += Sij(x, 1);
				Bi(x, 2) += Sij(x, 2);
			}
#endif
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
			writeLightcones(sim, cosmo, fourpiG, a, tau, dtau, dtau_old, dtau_older, maxvel[0], cycle, h5filename + sim.basename_lightcone, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &chi, &Bi, &Sij, &BiFT, &SijFT, &plan_Bi, &plan_Sij, lcbuffer, lclat, done_hij);
		else done_hij = 0;

#ifdef BENCHMARK
		lightcone_output_time += MPI_Wtime() - ref_time;
		ref_time = MPI_Wtime();
#endif

		// snapshot output
		if (snapcount < sim.num_snapshot && 1. / a < sim.z_snapshot[snapcount] + 1.)
		{
			COUT << COLORTEXT_CYAN << " writing snapshot" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;



#ifdef CHECK_B
			writeSnapshots(sim, cosmo, fourpiG, a, done_hij, snapcount, h5filename + sim.basename_snapshot, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &pi_k, &pi_v_k, &chi, &Bi, &T00_Kess, &T0i_Kess, &Tij_Kess, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, &Bi_check, &BiFT_check, &plan_Bi_check);
#else
			writeSnapshots(sim, cosmo, fourpiG, a, done_hij, snapcount, h5filename + sim.basename_snapshot, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &pi_k, &pi_v_k, &chi, &Bi, &T00_Kess, &T0i_Kess, &Tij_Kess, &source, &Sij, &scalarFT, &BiFT, &SijFT, &plan_phi, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
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
#ifdef COSIRA_HACK
			if (sim.gr_flag > 0)
			{
				COUT << COLORTEXT_CYAN << " writing power spectra" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;

#ifdef CHECK_B
				writeSpectra(sim, cosmo, fourpiG, a, pkcount, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &pi_k, &pi_v_k, &chi, &Bi, &source, &Sij, &scalarFT, &scalarFT_pi, &scalarFT_pi_v, &BiFT, &SijFT, &plan_phi, &plan_pi_k, &plan_pi_v_k, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, &Bi_check, &BiFT_check, &plan_Bi_check);
#else
				writeSpectra(sim, cosmo, fourpiG, a, pkcount, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &pi_k, &pi_v_k, &chi, &Bi, &source, &Sij, &scalarFT, &scalarFT_pi, &scalarFT_pi_v, &BiFT, &SijFT, &plan_phi, &plan_pi_k, &plan_pi_v_k, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, &Bi_check, &BiFT_check, &plan_Bi_check);
#endif

				pkcount++;

				COUT << " COSIRA hack: computing gauge transformation..." << endl;

				long tmp2;
				float * pcldata = NULL;
				gsl_spline * tk_d1 = NULL;
				gsl_spline * tk_d2 = NULL;
				gsl_spline * tk_t1 = NULL;
				gsl_spline * tk_t2 = NULL;
				double * temp1 = NULL;
				Field<Real> * psource = &source;

				loadTransferFunctions(class_background, class_perturbs, class_spectra, tk_d1, tk_t1, "tot", sim.boxsize, (1. / a) - 1., cosmo.h);
				loadTransferFunctions(class_background, class_perturbs, class_spectra, tk_d2, tk_t2, NULL, sim.boxsize, (1. / a) - 1., cosmo.h);

				temp1 = (double *) malloc(tk_d1->size * sizeof(double));

				for (i = 0; i < tk_d2->size; i++)
					temp1[i] = -tk_d2->y[i] * M_PI * sqrt(Pk_primordial(tk_d2->x[i] * cosmo.h / sim.boxsize, ic) / tk_d2->x[i]) / tk_d2->x[i];
				gsl_spline_free(tk_t2);
				tk_t2 = gsl_spline_alloc(gsl_interp_cspline, tk_d2->size);
				gsl_spline_init(tk_t2, tk_d2->x, temp1, tk_d2->size);
				gsl_spline_free(tk_d2);

				tmp = 3. * Hconf(a, fourpiG, cosmo)  * M_PI;

				for (i = 0; i < tk_d1->size; i++)
					temp1[i] = 3. * tk_t2->y[i] - tmp * tk_t1->y[i] * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i] / tk_d1->x[i] / tk_d1->x[i];

				gsl_spline_free(tk_d1);
				gsl_spline_free(tk_t1);
				tk_d1 = gsl_spline_alloc(gsl_interp_cspline, tk_t2->size);
				gsl_spline_init(tk_d1, tk_t2->x, temp1, tk_t2->size);
				gsl_spline_free(tk_t2);
				free(temp1);

				if (ic.flags & ICFLAG_CORRECT_DISPLACEMENT)
				{
					loadHomogeneousTemplate(ic.pclfile[0], tmp2, pcldata);
					generateCICKernel(source, tmp2, pcldata, ic.numtile[0]);
					free(pcldata);
				}
				else
					generateCICKernel(source);
				plan_source.execute(FFT_FORWARD);
				generateDisplacementField(scalarFT, 0., tk_d1, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE);
				gsl_spline_free(tk_d1);
				plan_source.execute(FFT_BACKWARD);
				source.updateHalo();

				i = MAX;
				pcls_cdm.moveParticles(displace_pcls_ic_basic, 1., &psource, 1, NULL, &tmp, &i, 1);
				COUT << " Poisson gauge -> N-body gauge, cdm particles maximum displacement = " << tmp * sim.numpts << " lattice units." << endl;
				sim.gr_flag = 0;
			}
#endif // COSIRA_HACK
#ifdef NMGAUGE_HACK
			if (sim.gr_flag == 0)
			{
				COUT << COLORTEXT_CYAN << " writing power spectra" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;

#ifdef CHECK_B
				writeSpectra(sim, cosmo, fourpiG, a, pkcount, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &pi_k, &pi_v_k, &chi, &Bi, &source, &Sij, &scalarFT, &scalarFT_pi, &scalarFT_pi_v, &BiFT, &SijFT, &plan_phi, &plan_pi_k, &plan_pi_v_k, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, &Bi_check, &BiFT_check, &plan_Bi_check);
#else
				writeSpectra(sim, cosmo, fourpiG, a, pkcount, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &pi_k, &pi_v_k, &chi, &Bi, &source, &Sij, &scalarFT, &scalarFT_pi, &scalarFT_pi_v, &BiFT, &SijFT, &plan_phi, &plan_pi_k, &plan_pi_v_k, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, &Bi_check, &BiFT_check, &plan_Bi_check);
#endif

				pkcount++;

				COUT << " Nm-gauge hack: computing gauge transformation..." << endl;

				long tmp2;
				float * pcldata = NULL;
				gsl_spline * tk_d1 = NULL;
				gsl_spline * tk_t1 = NULL;
				double * temp1 = NULL;
				Field<Real> * psource = &source;

				loadTransferFunctions("gaugetrafo.dat", tk_d1, tk_t1, "L", sim.boxsize, cosmo.h);

				temp1 = (double *) malloc(tk_d1->size * sizeof(double));

				for (i = 0; i < tk_d1->size; i++)
					temp1[i] = tk_d1->y[i] * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) * cosmo.h / sim.boxsize;
				gsl_spline_free(tk_t1);
				tk_t1 = gsl_spline_alloc(gsl_interp_cspline, tk_d1->size);
				gsl_spline_init(tk_t1, tk_d1->x, temp1, tk_d1->size);
				gsl_spline_free(tk_d1);
				free(temp1);

				if (ic.flags & ICFLAG_CORRECT_DISPLACEMENT)
				{
					loadHomogeneousTemplate(ic.pclfile[0], tmp2, pcldata);
					generateCICKernel(source, tmp2, pcldata, ic.numtile[0]);
					free(pcldata);
				}
				else
					generateCICKernel(source);
				plan_source.execute(FFT_FORWARD);
				generateDisplacementField(scalarFT, 0., tk_t1, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE);
				gsl_spline_free(tk_t1);
				plan_source.execute(FFT_BACKWARD);
				source.updateHalo();

				i = MAX;
				pcls_cdm.moveParticles(displace_pcls_ic_basic, 1., &psource, 1, NULL, &tmp, &i, 1);
				COUT << " Newtonian motion gauge -> N-body gauge, cdm particles maximum displacement = " << tmp * sim.numpts << " lattice units." << endl;
			}
#endif // NMGAUGE_HACK
			COUT << COLORTEXT_CYAN << " writing power spectra" << COLORTEXT_RESET << " at z = " << ((1./a) - 1.) <<  " (cycle " << cycle << "), tau/boxsize = " << tau << endl;

#ifdef CHECK_B
			writeSpectra(sim, cosmo, fourpiG, a, pkcount, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &pi_k, &pi_v_k, &chi, &Bi, &source, &Sij, &scalarFT, &scalarFT_pi, &scalarFT_pi_v, &BiFT, &SijFT, &plan_phi, &plan_pi_k, &plan_pi_v_k, &plan_chi, &plan_Bi, &plan_source, &plan_Sij, &Bi_check, &BiFT_check, &plan_Bi_check);
#else
			writeSpectra(sim, cosmo, fourpiG, a, pkcount, &pcls_cdm, &pcls_b, pcls_ncdm, &phi, &pi_k, &pi_v_k, &chi, &Bi, &source, &Sij, &scalarFT, &scalarFT_pi, &scalarFT_pi_v, &BiFT, &SijFT, &plan_phi, &plan_pi_k, &plan_pi_v_k, &plan_chi, &plan_Bi, &plan_source, &plan_Sij);
#endif

			pkcount++;
		}

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

#ifdef MULTISTEP_PROJECTION		// multistep projection: non-cold species are projected at each of their update steps
		projection_init(&source);
#ifdef HAVE_CLASS
		if (sim.radiation_flag > 0)
		{
			tmp = a;
			for (j = 0; j < 2 * numsteps; j++)
				rungekutta4bg(tmp, fourpiG, cosmo, 0.5 * dtau / numsteps);
			projection_T00_project(class_background, class_perturbs, class_spectra, source, scalarFT, &plan_source, sim, ic, cosmo, fourpiG, tmp);
		}
#endif
#endif


double a_kess=a;
double Hconf_old= Hconf(a_kess, fourpiG, cosmo);
// if(cycle==0)
// {
// 	for (i=0;i<sim.nKe_numsteps;i++)
// 	{
// 		dtau_old=0.1;
// 		update_pi_k_v( 0.5 * dtau_old/ sim.nKe_numsteps, dx, a_kess, phi, phi_old, chi, chi_old, pi_k, pi_v_k,cosmo.Omega_kessence, cosmo.w_kessence, cosmo.cs2_kessence, Hconf(a_kess, fourpiG, cosmo),Hconf_old);
// 		pi_v_k.updateHalo();
// 	}
// }
if(dtau_old>0)
{
	for (i=0;i<sim.nKe_numsteps;i++)
	{
		update_pi_k( dtau_old  / sim.nKe_numsteps,phi, pi_k, pi_v_k);
		pi_k.updateHalo();
		rungekutta4bg(a_kess, fourpiG, cosmo,  dtau  / sim.nKe_numsteps / 2.0);
		update_pi_k_v(dtau_old/ sim.nKe_numsteps ,
									dx,a_kess,phi,phi_old,chi,chi_old,pi_k, pi_v_k,cosmo.Omega_kessence,cosmo.w_kessence,cosmo.cs2_kessence,
									Hconf(a_kess, fourpiG, cosmo),Hconf_old);
		pi_v_k.updateHalo();
		rungekutta4bg(a_kess, fourpiG, cosmo,  dtau  / sim.nKe_numsteps / 2.0 );
		
	}
}
//end of Kessence + background:x
//:::::::::::::::::

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
#ifdef MULTISTEP_PROJECTION		// multistep projection: non-cold species are projected at each of their update steps
					if (sim.radiation_flag == 0 || a >= 1. / (sim.z_switch_deltancdm[i] + 1.))
					{
						if (sim.gr_flag > 0)
						{
							tmp = a;
							rungekutta4bg(tmp, fourpiG, cosmo, 0.5 * dtau / numsteps);
							projection_T00_project(pcls_ncdm+i, &source, tmp, &phi, 1. / (double) numsteps_ncdm[i]);
						}
						else
						{
							tmp = pcls_ncdm[i].parts_info()->mass;
							pcls_ncdm[i].parts_info()->mass /= (double) numsteps_ncdm[i];
							scalarProjectionCIC_project(pcls_ncdm+i, &source);
							pcls_ncdm[i].parts_info()->mass = tmp;
						}
					}
#ifdef BENCHMARK
						projection_time += MPI_Wtime() - ref2_time;
						ref2_time = MPI_Wtime();
#endif
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
#ifdef MULTISTEP_PROJECTION		// multistep projection: non-cold species are projected at each of their update steps
					if (sim.radiation_flag == 0 || a >= 1. / (sim.z_switch_deltancdm[i] + 1.))
					{
						if (sim.gr_flag > 0)
						{
							tmp = a;
							rungekutta4bg(tmp, fourpiG, cosmo, 0.5 * dtau / numsteps);
							projection_T00_project(pcls_ncdm+i, &source, tmp, &phi, 1. / (double) numsteps_ncdm[i]);
						}
						else
						{
							tmp = pcls_ncdm[i].parts_info()->mass;
							pcls_ncdm[i].parts_info()->mass /= (double) numsteps_ncdm[i];
							scalarProjectionCIC_project(pcls_ncdm+i, &source);
							pcls_ncdm[i].parts_info()->mass = tmp;
						}
					}
#ifdef BENCHMARK
						projection_time += MPI_Wtime() - ref2_time;
						ref2_time = MPI_Wtime();
#endif
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
					hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, pi_k, pi_v_k, chi, Bi_check, a, tau, dtau, cycle);
				}
				else
#endif
				hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, pi_k, pi_v_k, chi, Bi, a, tau, dtau, cycle);
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
				hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi,pi_k, pi_v_k, chi, Bi_check, a, tau, dtau, cycle, restartcount);
			}
			else
#endif
			hibernate(sim, ic, cosmo, &pcls_cdm, &pcls_b, pcls_ncdm, phi, pi_k, pi_v_k, chi, Bi, a, tau, dtau, cycle, restartcount);
			restartcount++;
		}

		dtau_older = dtau_old;
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

#ifdef BENCHMARK
		ref_time = MPI_Wtime();
#endif

	if (sim.num_lightcone > 0) // metric lightcone output
	{
		for (i = 0, j = 0; i < sim.num_lightcone; i++)
		{
			j |= sim.out_lightcone[i];
		}

		if (j & MASK_PHI)
		{
			COUT << COLORTEXT_CYAN << " writing phi lightcone data" << COLORTEXT_RESET << endl;

			for (i = 0; i < LIGHTCONE_THICKNESS; i++)
			{
				if (LIGHTCONE_THICKNESS > 1)
					sprintf(filename, "_phi_%d.h5", i);
				else
					sprintf(filename, "_phi.h5");

				lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_PHI_OFFSET+i]->saveHDF5(h5filename + sim.basename_lightcone + filename);

				delete lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_PHI_OFFSET+i];
			}

			delete lclat[LIGHTCONE_PHI_OFFSET];
		}

		if (j & MASK_CHI)
		{
			COUT << COLORTEXT_CYAN << " writing chi lightcone data" << COLORTEXT_RESET << endl;

			for (i = 0; i < LIGHTCONE_THICKNESS; i++)
			{
				if (LIGHTCONE_THICKNESS > 1)
					sprintf(filename, "_chi_%d.h5", i);
				else
					sprintf(filename, "_chi.h5");

				lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_CHI_OFFSET+i]->saveHDF5(h5filename + sim.basename_lightcone + filename);

				delete lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_CHI_OFFSET+i];
			}

			delete lclat[LIGHTCONE_CHI_OFFSET];
		}

		if (j & MASK_B)
		{
			COUT << COLORTEXT_CYAN << " writing B lightcone data" << COLORTEXT_RESET << endl;

			for (i = 0; i < 3*LIGHTCONE_THICKNESS; i++)
			{
				if (LIGHTCONE_THICKNESS > 1)
					sprintf(filename, "_B%d_%d.h5", (i%3)+1, i/3);
				else
					sprintf(filename, "_B%d.h5", i+1);

				lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_B_OFFSET+i]->saveHDF5(h5filename + sim.basename_lightcone + filename);

				delete lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_B_OFFSET+i];
			}

			delete lclat[LIGHTCONE_B_OFFSET];
		}

		if (j & MASK_HIJ)
		{
			COUT << COLORTEXT_CYAN << " writing hij lightcone data" << COLORTEXT_RESET << endl;

			for (i = 0; i < 5*LIGHTCONE_THICKNESS; i++)
			{
				if (LIGHTCONE_THICKNESS > 1)
					sprintf(filename, "_h%d%d_%d.h5", ((i%5) < 3 ? 1 : 2), (i%5) + ((i%5) < 3 ? 1 : -1), i/5);
				else
					sprintf(filename, "_h%d%d.h5", (i < 3 ? 1 : 2), i + (i < 3 ? 1 : -1));

				lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+i]->saveHDF5(h5filename + sim.basename_lightcone + filename);

				delete lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+i];
			}

			delete lclat[LIGHTCONE_HIJ_OFFSET];
		}
	}

#ifdef HAVE_CLASS
#ifndef COSIRA_HACK
	if (sim.radiation_flag > 0)
#endif // COSIRA_HACK
		freeCLASSstructures(class_background, class_perturbs, class_spectra);
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
