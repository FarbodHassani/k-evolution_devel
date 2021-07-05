//////////////////////////
// output.hpp
//////////////////////////
//
// Output of snapshots, light cones and spectra
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris)
//
// Last modified: April 2019
//
//////////////////////////

#ifndef OUTPUT_HEADER
#define OUTPUT_HEADER

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

using namespace std;


//////////////////////////
// writeSnapshots
//////////////////////////
// Description:
//   output of snapshots
//
// Arguments:
//   sim            simulation metadata structure
//   cosmo          cosmological parameter structure
//   fourpiG        4 pi G (in code units)
//   a              scale factor
//   snapcount      snapshot index
//   h5filename     base name for HDF5 output file
//   pcls_cdm       pointer to particle handler for CDM
//   pcls_b         pointer to particle handler for baryons
//   pcls_ncdm      array of particle handlers for
//                  non-cold DM (may be set to NULL)
//   phi            pointer to allocated field
//   pi_k           pointer to allocated field
//   zeta         pointer to allocated field
//   chi            pointer to allocated field
//   Bi             pointer to allocated field
//   source         pointer to allocated field
//   Sij            pointer to allocated field
//   T00_Kess       pointer to allocated field
//   T0i_Kess       pointer to allocated field
//   Tij_Kess       pointer to allocated field
//   scalarFT       pointer to allocated field
//   BiFT           pointer to allocated field
//   SijFT          pointer to allocated field
//   plan_phi       pointer to FFT planner
//   plan_chi       pointer to FFT planner
//   plan_Bi        pointer to FFT planner
//   plan_source    pointer to FFT planner
//   plan_Sij       pointer to FFT planner
//   Bi_check       pointer to allocated field
//   BiFT_check     pointer to allocated field
//   plan_Bi_check  pointer to FFT planner
//   vi             pointer to allocated field
//
// Returns:
//
//////////////////////////

void writeSnapshots(metadata & sim, cosmology & cosmo, const double fourpiG, const double a, const double dtau_old, const int done_hij, const int snapcount, string h5filename,
#ifdef HAVE_CLASS_BG
gsl_spline * H_spline, gsl_interp_accel * acc,
#endif
Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_cdm, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_b, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_ncdm, Field<Real> * phi, Field<Real> * pi_k, Field<Real> * zeta, Field<Real> * chi, Field<Real> * Bi, Field<Real> * T00_Kess, Field<Real> * T0i_Kess, Field<Real> * Tij_Kess, Field<Real> * source, Field<Real> * Sij, Field<Cplx> * scalarFT, Field<Cplx> * BiFT, Field<Cplx> * SijFT, PlanFFT<Cplx> * plan_phi, PlanFFT<Cplx> * plan_chi, PlanFFT<Cplx> * plan_Bi, PlanFFT<Cplx> * plan_source, PlanFFT<Cplx> * plan_Sij
#ifdef CHECK_B
, Field<Real> * Bi_check, Field<Cplx> * BiFT_check, PlanFFT<Cplx> * plan_Bi_check
#endif
#ifdef VELOCITY
, Field<Real> * vi
#endif
)
{
	char filename[2*PARAM_MAX_LENGTH+24];
	char buffer[64];
	int i;
	gadget2_header hdr;
	Site x(phi->lattice());
	Real divB, curlB, divh, traceh, normh;
	double dtau_pos = 0.;

	sprintf(filename, "%03d", snapcount);

#ifdef EXTERNAL_IO
	while (ioserver.openOstream()== OSTREAM_FAIL);

	if (sim.out_snapshot & MASK_PCLS)
	{
		pcls_cdm->saveHDF5_server_open(h5filename + filename + "_cdm");
		if (sim.baryon_flag)
			pcls_b->saveHDF5_server_open(h5filename + filename + "_b");
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			if (sim.numpcl[1+sim.baryon_flag+i] == 0) continue;
			sprintf(buffer, "_ncdm%d", i);
			pcls_ncdm[i].saveHDF5_server_open(h5filename + filename + buffer);
		}
	}

	if (sim.out_snapshot & MASK_T00)
		source->saveHDF5_server_open(h5filename + filename + "_T00");

#ifdef VELOCITY
	if (sim.out_snapshot & MASK_VEL)
		vi->saveHDF5_server_open(h5filename + filename + "_v");
#endif

	if (sim.out_snapshot & MASK_B)
		Bi->saveHDF5_server_open(h5filename + filename + "_B");

	if (sim.out_snapshot & MASK_PHI)
		phi->saveHDF5_server_open(h5filename + filename + "_phi");

	if (sim.out_snapshot & MASK_PI_K)
			pi_k->saveHDF5_server_open(h5filename + filename + "_pi_k");

	if (sim.out_snapshot & MASK_zeta)
			zeta->saveHDF5_server_open(h5filename + filename + "_zeta");
	//Kessence end


	if (sim.out_snapshot & MASK_CHI)
		chi->saveHDF5_server_open(h5filename + filename + "_chi");

	if (sim.out_snapshot & MASK_HIJ)
		Sij->saveHDF5_server_open(h5filename + filename + "_hij");

#ifdef CHECK_B
	if (sim.out_snapshot & MASK_B)
		Bi_check->saveHDF5_server_open(h5filename + filename + "_B_check");
#endif
#endif

	if (sim.out_snapshot & MASK_RBARE || sim.out_snapshot & MASK_POT)
	{
		projection_init(source);
		scalarProjectionCIC_project(pcls_cdm, source);
		if (sim.baryon_flag)
			scalarProjectionCIC_project(pcls_b, source);
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			if (sim.numpcl[1+sim.baryon_flag+i] == 0) continue;
			scalarProjectionCIC_project(pcls_ncdm+i, source);
		}
		scalarProjectionCIC_comm(source);
	}

	if (sim.out_snapshot & MASK_RBARE)
	{
		if (sim.downgrade_factor > 1)
			source->saveHDF5_coarseGrain3D(h5filename + filename + "_rhoN.h5", sim.downgrade_factor);
		else
			source->saveHDF5(h5filename + filename + "_rhoN.h5");
	}

	if (sim.out_snapshot & MASK_POT)
	{
		plan_source->execute(FFT_FORWARD);
		solveModifiedPoissonFT(*scalarFT, *scalarFT, fourpiG / a);
		plan_source->execute(FFT_BACKWARD);
		if (sim.downgrade_factor > 1)
			source->saveHDF5_coarseGrain3D(h5filename + filename + "_psiN.h5", sim.downgrade_factor);
		else
			source->saveHDF5(h5filename + filename + "_psiN.h5");
	}

	if (sim.out_snapshot & MASK_T00)
	{
		projection_init(source);
		if (sim.gr_flag > 0)
		{
			projection_T00_project(pcls_cdm, source, a, phi);
			if (sim.baryon_flag)
				projection_T00_project(pcls_b, source, a, phi);
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (sim.numpcl[1+sim.baryon_flag+i] == 0) continue;
				projection_T00_project(pcls_ncdm+i, source, a, phi);
			}
		}
		else
		{
			scalarProjectionCIC_project(pcls_cdm, source);
			if (sim.baryon_flag)
				scalarProjectionCIC_project(pcls_b, source);
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (sim.numpcl[1+sim.baryon_flag+i] == 0) continue;
				scalarProjectionCIC_project(pcls_ncdm+i, source);
			}
		}
		projection_T00_comm(source);
#ifdef EXTERNAL_IO
		source->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if (sim.downgrade_factor > 1)
			source->saveHDF5_coarseGrain3D(h5filename + filename + "_T00.h5", sim.downgrade_factor);
		else
			source->saveHDF5(h5filename + filename + "_T00.h5");
#endif
	}

#ifdef VELOCITY
	if (sim.out_snapshot & MASK_VEL)
	{
#ifdef EXTERNAL_IO
		vi->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if (sim.downgrade_factor > 1)
			vi->saveHDF5_coarseGrain3D(h5filename + filename + "_v.h5", sim.downgrade_factor);
		else
			vi->saveHDF5(h5filename + filename + "_v.h5");
#endif
	}
#endif

	if (sim.out_snapshot & MASK_B)
	{
		if (sim.gr_flag == 0)
		{
			plan_Bi->execute(FFT_BACKWARD);
		}
		for (x.first(); x.test(); x.next())
		{
			(*Bi)(x,0) /= a * a * sim.numpts;
			(*Bi)(x,1) /= a * a * sim.numpts;
			(*Bi)(x,2) /= a * a * sim.numpts;
		}
		Bi->updateHalo();

		computeVectorDiagnostics(*Bi, divB, curlB);
		COUT << " B diagnostics: max |divB| = " << divB << ", max |curlB| = " << curlB << endl;

#ifdef EXTERNAL_IO
		Bi->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if (sim.downgrade_factor > 1)
			Bi->saveHDF5_coarseGrain3D(h5filename + filename + "_B.h5", sim.downgrade_factor);
		else
			Bi->saveHDF5(h5filename + filename + "_B.h5");
#endif

		if (sim.gr_flag > 0)
		{
			plan_Bi->execute(FFT_BACKWARD);
			Bi->updateHalo();
		}
	}

	if (sim.out_snapshot & MASK_PHI)
#ifdef EXTERNAL_IO
		phi->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if (sim.downgrade_factor > 1)
			phi->saveHDF5_coarseGrain3D(h5filename + filename + "_phi.h5", sim.downgrade_factor);
		else
			phi->saveHDF5(h5filename + filename + "_phi.h5");
#endif


	if (sim.out_snapshot & MASK_PI_K)
#ifdef EXTERNAL_IO
		pi_k->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if (sim.downgrade_factor > 1)
			pi_k->saveHDF5_coarseGrain3D(h5filename + filename + "_pi_k.h5", sim.downgrade_factor);
		else
			pi_k->saveHDF5(h5filename + filename + "_pi_k.h5");
#endif

	if (sim.out_snapshot & MASK_CHI)
#ifdef EXTERNAL_IO
		chi->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if (sim.downgrade_factor > 1)
			chi->saveHDF5_coarseGrain3D(h5filename + filename + "_chi.h5", sim.downgrade_factor);
		else
			chi->saveHDF5(h5filename + filename + "_chi.h5");
#endif

	if (sim.out_snapshot & MASK_HIJ)
	{
		if (done_hij == 0)
		{
			projectFTtensor(*SijFT, *SijFT);
			plan_Sij->execute(FFT_BACKWARD);
			Sij->updateHalo();
		}

		computeTensorDiagnostics(*Sij, divh, traceh, normh);
		COUT << " GW diagnostics: max |divh| = " << divh << ", max |traceh| = " << traceh << ", max |h| = " << normh << endl;

#ifdef EXTERNAL_IO
		Sij->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if (sim.downgrade_factor > 1)
			Sij->saveHDF5_coarseGrain3D(h5filename + filename + "_hij.h5", sim.downgrade_factor);
		else
			Sij->saveHDF5(h5filename + filename + "_hij.h5");
#endif
	}

	if (sim.out_snapshot & MASK_TIJ)
	{
		projection_init(Sij);
		projection_Tij_project(pcls_cdm, Sij, a, phi);
		if (sim.baryon_flag)
			projection_Tij_project(pcls_b, Sij, a, phi);
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			if (sim.numpcl[1+sim.baryon_flag+i] == 0) continue;
			projection_Tij_project(pcls_ncdm+i, Sij, a, phi);
		}
		projection_Tij_comm(Sij);

		if (sim.downgrade_factor > 1)
			Sij->saveHDF5_coarseGrain3D(h5filename + filename + "_Tij.h5", sim.downgrade_factor);
		else
			Sij->saveHDF5(h5filename + filename + "_Tij.h5");
	}

	if (sim.out_snapshot & MASK_P)
	{
		projection_init(Bi);
		projection_T0i_project(pcls_cdm, Bi, phi);
		if (sim.baryon_flag)
			projection_T0i_project(pcls_b, Bi, phi);
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			if (sim.numpcl[1+sim.baryon_flag+i] == 0) continue;
			projection_T0i_project(pcls_ncdm+i, Bi, phi);
		}
		projection_T0i_comm(Bi);
		if (sim.downgrade_factor > 1)
			Bi->saveHDF5_coarseGrain3D(h5filename + filename + "_p.h5", sim.downgrade_factor);
		else
			Bi->saveHDF5(h5filename + filename + "_p.h5");
		if (sim.gr_flag > 0)
		{
			plan_Bi->execute(FFT_BACKWARD);
			Bi->updateHalo();
		}
	}

#ifdef CHECK_B
	if (sim.out_snapshot & MASK_B)
	{
		if (sim.vector_flag == VECTOR_PARABOLIC)
		{
			projection_init(Bi_check);
			projection_T0i_project(pcls_cdm, Bi_check, phi);
			if (sim.baryon_flag)
				projection_T0i_project(pcls_b, Bi_check, phi);
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (sim.numpcl[1+sim.baryon_flag+i] == 0) continue;
				projection_T0i_project(pcls_ncdm+i, Bi_check, phi);
			}
			projection_T0i_comm(Bi_check);
			plan_Bi_check->execute(FFT_FORWARD);
			projectFTvector(*BiFT_check, *BiFT_check, fourpiG / (double) sim.numpts / (double) sim.numpts);
		}
		plan_Bi_check->execute(FFT_BACKWARD);

		for (x.first(); x.test(); x.next())
		{
			(*Bi_check)(x,0) /= a * a * sim.numpts;
			(*Bi_check)(x,1) /= a * a * sim.numpts;
			(*Bi_check)(x,2) /= a * a * sim.numpts;
		}
#ifdef EXTERNAL_IO
		Bi_check->saveHDF5_server_write(NUMBER_OF_IO_FILES);
#else
		if (sim.downgrade_factor > 1)
			Bi_check->saveHDF5_coarseGrain3D(h5filename + filename + "_B_check.h5", sim.downgrade_factor);
		else
			Bi_check->saveHDF5(h5filename + filename + "_B_check.h5");
#endif
	}
#endif
//Kessence
if (sim.out_snapshot & MASK_T_KESS)
{
	// COUT<<"output kesssence energy"<<endl;
	#ifdef EXTERNAL_IO
			T00_Kess->saveHDF5_server_write(NUMBER_OF_IO_FILES);
			// T0i_Kess->saveHDF5_server_write(NUMBER_OF_IO_FILES);
			// Tij_Kess->saveHDF5_server_write(NUMBER_OF_IO_FILES);
	#else
			if (sim.downgrade_factor > 1)
			{
				T00_Kess->saveHDF5_coarseGrain3D(h5filename + filename + "_T00_Kess.h5", sim.downgrade_factor);
				// T0i_Kess->saveHDF5_coarseGrain3D(h5filename + filename + "_T0i_Kess.h5", sim.downgrade_factor);
				// Tij_Kess->saveHDF5_coarseGrain3D(h5filename + filename + "_Tij_Kess.h5", sim.downgrade_factor);
			}
			else
			{
				T00_Kess->saveHDF5(h5filename + filename + "_T00_Kess.h5");
				// T0i_Kess->saveHDF5(h5filename + filename + "_T0i_Kess.h5");
				// Tij_Kess->saveHDF5(h5filename + filename + "_Tij_Kess.h5");
			}
	#endif
}
//Kessence end
	if (sim.out_snapshot & MASK_GADGET)
	{
		if (sim.out_snapshot & MASK_MULTI)
			hdr.num_files = parallel.grid_size()[1];
		else
			hdr.num_files = 1;
		hdr.Omega0 = cosmo.Omega_m;
		hdr.OmegaLambda = cosmo.Omega_Lambda;
		hdr.HubbleParam = cosmo.h;
		hdr.BoxSize = sim.boxsize / GADGET_LENGTH_CONVERSION;
		hdr.flag_sfr = 0;
		hdr.flag_cooling = 0;
		hdr.flag_feedback = 0;
		hdr.flag_age = 0;
		hdr.flag_metals = 0;
		for (i = 0; i < 256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8 - 2 * 4 - 6 * 4; i++)
			hdr.fill[i] = 0;
		for (i = 0; i < 6; i++)
		{
			hdr.npart[i] = 0;
			hdr.npartTotal[i] = 0;
			hdr.npartTotalHW[i] = 0;
			hdr.mass[i] = 0.;
		}

#ifdef EXACT_OUTPUT_REDSHIFTS
		hdr.time = 1. / (sim.z_snapshot[snapcount] + 1.);
		hdr.redshift = sim.z_snapshot[snapcount];
		dtau_pos = (hdr.time - a) / a / Hconf(a, fourpiG,//TODO_EB
		#ifdef HAVE_CLASS_BG
			H_spline, acc
		#else
			cosmo
		#endif
		);
#else
		hdr.time = a;
		hdr.redshift = (1./a) - 1.;
#endif

		if (sim.tracer_factor[0] > 0)
		{
			hdr.npart[1] = (uint32_t) (((sim.numpcl[0] % sim.tracer_factor[0]) ? (1 + (sim.numpcl[0] / sim.tracer_factor[0])) : (sim.numpcl[0] / sim.tracer_factor[0])) % (1ll << 32));
			hdr.npartTotal[1] = hdr.npart[1];
			hdr.npartTotalHW[1] = (uint32_t) (((sim.numpcl[0] % sim.tracer_factor[0]) ? (1 + (sim.numpcl[0] / sim.tracer_factor[0])) : (sim.numpcl[0] / sim.tracer_factor[0])) / (1ll << 32));
			if (sim.baryon_flag)
				hdr.mass[1] = (double) sim.tracer_factor[0] * C_RHO_CRIT * cosmo.Omega_cdm * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[0] / GADGET_MASS_CONVERSION;
			else
				hdr.mass[1] = (double) sim.tracer_factor[0] * C_RHO_CRIT * (cosmo.Omega_cdm + cosmo.Omega_b) * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[0] / GADGET_MASS_CONVERSION;

			if (hdr.npartTotalHW[1] != 0 && hdr.num_files == 1)
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": number of particles (" << (sim.numpcl[0] / sim.tracer_factor[0]) << ") in Gadget2 file exceeds limit (4294967295). Try using multi-Gadget2 output format." << endl;
			}
			else
				pcls_cdm->saveGadget2(h5filename + filename + "_cdm", hdr, sim.tracer_factor[0], dtau_pos, dtau_pos + 0.5 * dtau_old, phi);
		}

		if (sim.baryon_flag && sim.tracer_factor[1] > 0)
		{
			hdr.npart[1] = (uint32_t) (((sim.numpcl[1] % sim.tracer_factor[1]) ? (1 + (sim.numpcl[1] / sim.tracer_factor[1])) : (sim.numpcl[1] / sim.tracer_factor[1])) % (1ll << 32));
			hdr.npartTotal[1] = hdr.npart[1];
			hdr.npartTotalHW[1] = (uint32_t) (((sim.numpcl[1] % sim.tracer_factor[1]) ? (1 + (sim.numpcl[1] / sim.tracer_factor[1])) : (sim.numpcl[1] / sim.tracer_factor[1])) / (1ll << 32));
			hdr.mass[1] = (double) sim.tracer_factor[1] * C_RHO_CRIT * cosmo.Omega_b * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[1] / GADGET_MASS_CONVERSION;
			if (hdr.npartTotalHW[1] != 0 && hdr.num_files == 1)
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": number of particles (" << (sim.numpcl[1] / sim.tracer_factor[1]) << ") in Gadget2 file exceeds limit (4294967295). Try using multi-Gadget2 output format." << endl;
			}
			else
				pcls_b->saveGadget2(h5filename + filename + "_b", hdr, sim.tracer_factor[1], dtau_pos, dtau_pos + 0.5 * dtau_old, phi);
		}

		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			if (sim.numpcl[1+sim.baryon_flag+i] == 0 || sim.tracer_factor[i+1+sim.baryon_flag] == 0) continue;
			sprintf(buffer, "_ncdm%d", i);
			hdr.npart[1] = (uint32_t) (((sim.numpcl[i+1+sim.baryon_flag] % sim.tracer_factor[i+1+sim.baryon_flag]) ? (1 + (sim.numpcl[i+1+sim.baryon_flag] / sim.tracer_factor[i+1+sim.baryon_flag])) : (sim.numpcl[i+1+sim.baryon_flag] / sim.tracer_factor[i+1+sim.baryon_flag])) % (1ll << 32));
			hdr.npartTotal[1] = hdr.npart[1];
			hdr.npartTotalHW[1] = (uint32_t) (((sim.numpcl[i+1+sim.baryon_flag] % sim.tracer_factor[i+1+sim.baryon_flag]) ? (1 + (sim.numpcl[i+1+sim.baryon_flag] / sim.tracer_factor[i+1+sim.baryon_flag])) : (sim.numpcl[i+1+sim.baryon_flag] / sim.tracer_factor[i+1+sim.baryon_flag])) / (1ll << 32));
			hdr.mass[1] = (double) sim.tracer_factor[i+1+sim.baryon_flag] * C_RHO_CRIT * cosmo.Omega_ncdm[i] * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[i+1+sim.baryon_flag] / GADGET_MASS_CONVERSION;
			if (hdr.npartTotalHW[1] != 0 && hdr.num_files == 1)
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": number of particles (" << (sim.numpcl[i+1+sim.baryon_flag] / sim.tracer_factor[i+1+sim.baryon_flag]) << ") in Gadget2 file exceeds limit (4294967295). Try using multi-Gadget2 output format." << endl;
			}
			else
				pcls_ncdm[i].saveGadget2(h5filename + filename + buffer, hdr, sim.tracer_factor[i+1+sim.baryon_flag], dtau_pos, dtau_pos + 0.5 * dtau_old, phi);
		}
	}

	if (sim.out_snapshot & MASK_PCLS)
	{
#ifdef EXTERNAL_IO
		pcls_cdm->saveHDF5_server_write();
		if (sim.baryon_flag)
			pcls_b->saveHDF5_server_write();
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			if (sim.numpcl[1+sim.baryon_flag+i] == 0) continue;
			pcls_ncdm[i].saveHDF5_server_write();
		}
#else
		pcls_cdm->saveHDF5(h5filename + filename + "_cdm", 1);
		if (sim.baryon_flag)
			pcls_b->saveHDF5(h5filename + filename + "_b", 1);
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			if (sim.numpcl[1+sim.baryon_flag+i] == 0) continue;
			sprintf(buffer, "_ncdm%d", i);
			pcls_ncdm[i].saveHDF5(h5filename + filename + buffer, 1);
		}
#endif
	}

#ifdef EXTERNAL_IO
	ioserver.closeOstream();
#endif
}


//////////////////////////
// writeLightcones
//////////////////////////
// Description:
//   output of light cones
//
// Arguments:
//   sim            simulation metadata structure
//   cosmo          cosmological parameter structure
//   fourpiG        4 pi G (in code units)
//   a              scale factor
//   tau            conformal time
//   dtau           conformal time step
//   dtau_old       conformal time step of previous cycle
//   maxvel         maximum cdm velocity
//   cycle          current simulation cycle
//   h5filename     base name for HDF5 output file
//   pcls_cdm       pointer to particle handler for CDM
//   pcls_b         pointer to particle handler for baryons
//   pcls_ncdm      array of particle handlers for
//                  non-cold DM (may be set to NULL)
//   phi            pointer to allocated field
//   chi            pointer to allocated field
//   Bi             pointer to allocated field
//   Sij            pointer to allocated field
//   BiFT           pointer to allocated field
//   SijFT          pointer to allocated field
//   plan_Bi        pointer to FFT planner
//   plan_Sij       pointer to FFT planner
//   done_hij       reference to tensor projection flag
//   IDbacklog      IDs of particles written in previous cycle
//
// Returns:
//
//////////////////////////

void writeLightcones(metadata & sim, cosmology & cosmo, const double fourpiG, const double a, const double tau, const double dtau, const double dtau_old, const double maxvel, const int cycle, string h5filename,
#ifdef HAVE_CLASS_BG
background & class_background, gsl_spline * H_spline, gsl_interp_accel * acc,
#endif
Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_cdm, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_b, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_ncdm, Field<Real> * phi, Field<Real> * chi, Field<Real> * Bi, Field<Real> * Sij, Field<Cplx> * BiFT, Field<Cplx> * SijFT, PlanFFT<Cplx> * plan_Bi, PlanFFT<Cplx> * plan_Sij, int & done_hij, set<long> * IDbacklog)
{
	int i, j, n, p;
	double d;
	double vertex[MAX_INTERSECTS][3];
	double domain[6];
	double pos[3];
	double s[2];
	char filename[2*PARAM_MAX_LENGTH+24];
	char buffer[268];
	FILE * outfile;
	gadget2_header hdr;
	set<long> IDprelog[MAX_PCL_SPECIES];
	long * IDcombuf;
	long * IDcombuf2;
	Site xsim;
	int done_B = 0;
#ifdef HAVE_HEALPIX
	int64_t pix, pix2, q;
	vector<int> pixbatch_id;
	vector<int> sender_proc;
	vector<int> pixbatch_size[3];
	vector<int> pixbatch_delim[3];
	int pixbatch_type;
	int commdir[2];
	Real * pixbuf[LIGHTCONE_MAX_FIELDS][9];
	Real * commbuf;
	int pixbuf_size[9];
	int pixbuf_reserve[9];
	int64_t bytes, bytes2, offset2;
	vector<MPI_Offset> offset;
	char ** outbuf = new char*[LIGHTCONE_MAX_FIELDS];
	healpix_header maphdr;
	double R[3][3];
	double w[3];
	double temp;
	int base_pos[3];
	int shell, shell_inner, shell_outer, shell_write;
	uint32_t blocksize;
	MPI_File mapfile;
	MPI_Status status;
	MPI_Datatype patch;
	int io_group_size;
	long pixcount = 0;

	for (j = 0; j < 9*LIGHTCONE_MAX_FIELDS; j++)
		pixbuf[j/9][j%9] = NULL;

	for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
		outbuf[j] = NULL;
#endif

	done_hij = 0;

	domain[0] = -0.5;
	domain[1] = phi->lattice().coordSkip()[1] - 0.5;
	domain[2] = phi->lattice().coordSkip()[0] - 0.5;
	for (i = 0; i < 3; i++)
		domain[i+3] = domain[i] + phi->lattice().sizeLocal(i) + 1.;

	for (i = 0; i < 6; i++)
		domain[i] /= (double) sim.numpts;

	for (i = 0; i < sim.num_lightcone; i++)
	{
		if (parallel.isRoot())
		{
			if (sim.num_lightcone > 1)
				sprintf(filename, "%s%s%d_info.dat", sim.output_path, sim.basename_lightcone, i);
			else
				sprintf(filename, "%s%s_info.dat", sim.output_path, sim.basename_lightcone);

			outfile = fopen(filename, "a");
			if (outfile == NULL)
			{
				cout << " error opening file for lightcone info!" << endl;
			}
			else if (cycle == 0)
			{
				if (sim.num_lightcone > 1)
					fprintf(outfile, "# information file for lightcone %d\n# geometric parameters:\n# vertex = (%f, %f, %f) Mpc/h\n# redshift = %f\n# distance = (%f - %f) Mpc/h\n# opening half-angle = %f degrees\n# direction = (%f, %f, %f)\n# cycle   tau/boxsize    a              pcl_inner        pcl_outer        metric_inner     metric_outer\n", i, sim.lightcone[i].vertex[0]*sim.boxsize, sim.lightcone[i].vertex[1]*sim.boxsize, sim.lightcone[i].vertex[2]*sim.boxsize, sim.lightcone[i].z, sim.lightcone[i].distance[0]*sim.boxsize, sim.lightcone[i].distance[1]*sim.boxsize, (sim.lightcone[i].opening > -1.) ? acos(sim.lightcone[i].opening) * 180. / M_PI : 180., sim.lightcone[i].direction[0], sim.lightcone[i].direction[1], sim.lightcone[i].direction[2]);
				else
					fprintf(outfile, "# information file for lightcone\n# geometric parameters:\n# vertex = (%f, %f, %f) Mpc/h\n# redshift = %f\n# distance = (%f - %f) Mpc/h\n# opening half-angle = %f degrees\n# direction = (%f, %f, %f)\n# cycle   tau/boxsize    a              pcl_inner        pcl_outer        metric_inner     metric_outer\n", sim.lightcone[i].vertex[0]*sim.boxsize, sim.lightcone[i].vertex[1]*sim.boxsize, sim.lightcone[i].vertex[2]*sim.boxsize, sim.lightcone[i].z, sim.lightcone[i].distance[0]*sim.boxsize, sim.lightcone[i].distance[1]*sim.boxsize, (sim.lightcone[i].opening > -1.) ? acos(sim.lightcone[i].opening) * 180. / M_PI : 180., sim.lightcone[i].direction[0], sim.lightcone[i].direction[1], sim.lightcone[i].direction[2]);
			}
		}

		d = particleHorizon(1. / (1. + sim.lightcone[i].z), fourpiG,
			#ifdef HAVE_CLASS_BG
			gsl_spline_eval(H_spline, 1., acc), class_background
			#else
			cosmo
			#endif
		);

		s[0] = d - tau - 0.5 * sim.covering[i] * dtau;
		s[1] = d - tau + 0.5 * sim.covering[i] * dtau_old;

#ifdef HAVE_HEALPIX
		shell_inner = (s[0] > sim.lightcone[i].distance[1]) ? ceil(s[0] * sim.numpts * sim.shellfactor[i]) : ceil(sim.lightcone[i].distance[1] * sim.numpts * sim.shellfactor[i]);
		shell_outer = (sim.lightcone[i].distance[0] > s[1]) ? floor(s[1] * sim.numpts * sim.shellfactor[i]) : (ceil(sim.lightcone[i].distance[0] * sim.numpts * sim.shellfactor[i])-1);
		if (shell_outer < shell_inner && s[1] > 0) shell_outer = shell_inner;

		maphdr.precision = sizeof(Real);
		maphdr.Ngrid = sim.numpts;
		maphdr.direction[0] = sim.lightcone[i].direction[0];
		maphdr.direction[1] = sim.lightcone[i].direction[1];
		maphdr.direction[2] = sim.lightcone[i].direction[2];
		maphdr.boxsize = sim.boxsize;
		memset((void *) maphdr.fill, 0, 256 - 5 * 4 - 5 * 8);

		xsim.initialize(phi->lattice());

		if (sim.lightcone[i].direction[0] == 0 && sim.lightcone[i].direction[1] == 0)
		{
			R[0][0] = sim.lightcone[i].direction[2];
			R[0][1] = 0;
			R[0][2] = 0;
			R[1][0] = 0;
			R[1][1] = 1;
			R[1][2] = 0;
			R[2][0] = 0;
			R[2][1] = 0;
			R[2][2] = sim.lightcone[i].direction[2];
		}
		else
		{
			temp = atan2(sim.lightcone[i].direction[1], sim.lightcone[i].direction[0]);
			R[0][0] = cos(temp) * sim.lightcone[i].direction[2];
			R[0][1] = -sin(temp);
			R[0][2] = sim.lightcone[i].direction[0];
			R[1][0] = sin(temp) * sim.lightcone[i].direction[2];
			R[1][1] = cos(temp);
			R[1][2] = sim.lightcone[i].direction[1];
			R[2][0] = -sqrt(1. - sim.lightcone[i].direction[2]*sim.lightcone[i].direction[2]);
			R[2][1] = 0;
			R[2][2] = sim.lightcone[i].direction[2];
		}
#endif

		if (sim.lightcone[i].distance[0] > s[0] && sim.lightcone[i].distance[1] <= s[1] && s[1] > 0.)
		{
			if (parallel.isRoot() && outfile != NULL)
			{
				fprintf(outfile, "%6d   %e   %e   %2.12f   %2.12f   %2.12f   %2.12f\n", cycle, tau, a, d - tau - 0.5 * dtau, d - tau + 0.5 * dtau_old, s[0], s[1]);
				fclose(outfile);

				if (sim.num_lightcone > 1)
					sprintf(filename, "%s%s%d_info.bin", sim.output_path, sim.basename_lightcone, i);
				else
					sprintf(filename, "%s%s_info.bin", sim.output_path, sim.basename_lightcone);

				outfile = fopen(filename, "a");
				if (outfile == NULL)
				{
					cout << " error opening file for lightcone info!" << endl;
				}
				else
				{
					((double *) buffer)[0] = tau;
					((double *) buffer)[1] = a;
					((double *) buffer)[2] = d - tau - 0.5 * dtau;
					((double *) buffer)[3] = d - tau + 0.5 * dtau_old;

					fwrite((const void *) &cycle, sizeof(int), 1, outfile);
					fwrite((const void *) buffer, sizeof(double), 4, outfile);
					fwrite((const void *) s, sizeof(double), 2, outfile);

					fclose(outfile);
				}
			}

#ifdef HAVE_HEALPIX
			bytes = 0;
			bytes2 = 0;

			for (j = 0; j < 9; j++)
				pixbuf_reserve[j] = PIXBUFFER;

			if (sim.out_lightcone[i] & MASK_PHI)
			{
				for (j = 0; j < 9; j++)
					pixbuf[LIGHTCONE_PHI_OFFSET][j] = (Real *) malloc(sizeof(Real) * PIXBUFFER);
			}

			if (sim.out_lightcone[i] & MASK_CHI)
			{
				for (j = 0; j < 9; j++)
					pixbuf[LIGHTCONE_CHI_OFFSET][j] = (Real *) malloc(sizeof(Real) * PIXBUFFER);
			}

			if (sim.out_lightcone[i] & MASK_B)
			{
				for (j = 0; j < 9; j++)
				{
					pixbuf[LIGHTCONE_B_OFFSET][j] = (Real *) malloc(sizeof(Real) * PIXBUFFER);
					pixbuf[LIGHTCONE_B_OFFSET+1][j] = (Real *) malloc(sizeof(Real) * PIXBUFFER);
					pixbuf[LIGHTCONE_B_OFFSET+2][j] = (Real *) malloc(sizeof(Real) * PIXBUFFER);
				}
			}

			if (sim.out_lightcone[i] & MASK_HIJ)
			{
				for (j = 0; j < 9; j++)
				{
					pixbuf[LIGHTCONE_HIJ_OFFSET][j] = (Real *) malloc(sizeof(Real) * PIXBUFFER);
					pixbuf[LIGHTCONE_HIJ_OFFSET+1][j] = (Real *) malloc(sizeof(Real) * PIXBUFFER);
					pixbuf[LIGHTCONE_HIJ_OFFSET+2][j] = (Real *) malloc(sizeof(Real) * PIXBUFFER);
					pixbuf[LIGHTCONE_HIJ_OFFSET+3][j] = (Real *) malloc(sizeof(Real) * PIXBUFFER);
					pixbuf[LIGHTCONE_HIJ_OFFSET+4][j] = (Real *) malloc(sizeof(Real) * PIXBUFFER);
				}
			}

			if (sim.gr_flag == 0 && sim.out_lightcone[i] & MASK_B && done_B == 0)
			{
				plan_Bi->execute(FFT_BACKWARD);
				Bi->updateHalo();
				done_B = 1;
			}

			if (sim.out_lightcone[i] & MASK_HIJ && done_hij == 0)
			{
				projectFTtensor(*SijFT, *SijFT);
				plan_Sij->execute(FFT_BACKWARD);
				Sij->updateHalo();
				done_hij = 1;
			}

			if ((shell_outer + 1 - shell_inner) > parallel.size())
			{
				shell_write = ((shell_outer + 1 - shell_inner) * parallel.rank() + parallel.size() - 1) / parallel.size();
				io_group_size = 0;
			}
			else
			{
				shell_write = ((shell_outer + 1 - shell_inner) * parallel.rank()) / parallel.size();
				io_group_size = (((shell_write+1) * parallel.size() + shell_outer - shell_inner) / (shell_outer + 1 - shell_inner)) - ((shell_write * parallel.size() + shell_outer - shell_inner) / (shell_outer + 1 - shell_inner));
			}

			for (shell = shell_inner; shell <= shell_outer; shell++)
			{
				maphdr.distance = (double) shell / (double) sim.numpts / sim.shellfactor[i];

				for (maphdr.Nside = sim.Nside[i][0]; maphdr.Nside < sim.Nside[i][1]; maphdr.Nside *= 2)
				{
					if (12. * maphdr.Nside * maphdr.Nside > sim.pixelfactor[i] * 4. * M_PI * maphdr.distance * maphdr.distance * sim.numpts * sim.numpts) break;
				}

				for (maphdr.Nside_ring = 2; 2.137937882409166 * sim.numpts * maphdr.distance / maphdr.Nside_ring > phi->lattice().sizeLocal(1) && 2.137937882409166 * sim.numpts * maphdr.distance / maphdr.Nside_ring > phi->lattice().sizeLocal(2) && maphdr.Nside_ring < maphdr.Nside; maphdr.Nside_ring *= 2);

				if (sim.lightcone[i].opening > 2./3.)
				{
					p = 1 + (int) floor(maphdr.Nside * sqrt(3. - 3. * sim.lightcone[i].opening));
					maphdr.Npix = 2 * p * (p+1);
				}
				else if (sim.lightcone[i].opening > -2./3.)
				{
					p = 1 + (int) floor(maphdr.Nside * (2. - 1.5 * sim.lightcone[i].opening));
					maphdr.Npix = 2 * maphdr.Nside * (maphdr.Nside+1) + (p-maphdr.Nside) * 4 * maphdr.Nside;
				}
				else if (sim.lightcone[i].opening > -1.)
				{
					p = (int) floor(maphdr.Nside * sqrt(3. + 3. * sim.lightcone[i].opening));
					maphdr.Npix = 12 * maphdr.Nside * maphdr.Nside - 2 * p * (p+1);
					p = 4 * maphdr.Nside - 1 - p;
				}
				else
				{
					maphdr.Npix = 12 * maphdr.Nside * maphdr.Nside;
					p = 4 * maphdr.Nside - 1;
				}

				pixbatch_size[0].push_back(maphdr.Nside / maphdr.Nside_ring);

				pixbatch_delim[1].push_back(p / pixbatch_size[0].back());
				pixbatch_delim[0].push_back((pixbatch_delim[1].back() > 0) ? pixbatch_delim[1].back()-1 : 0);
				pixbatch_delim[2].push_back(pixbatch_delim[1].back()+1);
				pixbatch_size[1].push_back((pixbatch_size[0].back() * (pixbatch_size[0].back()+1) + (2*pixbatch_size[0].back() - 1 - p%pixbatch_size[0].back()) * (p%pixbatch_size[0].back())) / 2);
				pixbatch_size[2].push_back(((p%pixbatch_size[0].back() + 1) * (p%pixbatch_size[0].back())) / 2);
				pixbatch_size[0].back() *= pixbatch_size[0].back();
				for (p = 0; p < 3; p++)
				{
					if (pixbatch_delim[p].back() <= maphdr.Nside_ring)
						pixbatch_delim[p].back() = 2 * pixbatch_delim[p].back() * (pixbatch_delim[p].back()+1);
					else if (pixbatch_delim[p].back() <= 3 * maphdr.Nside_ring)
						pixbatch_delim[p].back() = 2 * maphdr.Nside_ring * (maphdr.Nside_ring+1) + (pixbatch_delim[p].back()-maphdr.Nside_ring) * 4 * maphdr.Nside_ring;
					else if (pixbatch_delim[p].back() < 4 * maphdr.Nside_ring)
						pixbatch_delim[p].back() = 12 * maphdr.Nside_ring * maphdr.Nside_ring - 2 * (4 * maphdr.Nside_ring - 1 - pixbatch_delim[p].back()) * (4 * maphdr.Nside_ring - pixbatch_delim[p].back());
					else
						pixbatch_delim[p].back() = 12 * maphdr.Nside_ring * maphdr.Nside_ring;
				}

				if (pixbatch_size[1].back() == pixbatch_size[0].back())
					pixbatch_delim[0].back() = pixbatch_delim[1].back();

				for (j = 0; j < 9; j++)
					pixbuf_size[j] = 0;

				for (p = 0; p < pixbatch_delim[2].back(); p++)
				{
					pix2vec_ring64(maphdr.Nside_ring, p, w);

					base_pos[1] = (int) floor((maphdr.distance * (R[1][0] * w[0] + R[1][1] * w[1] + R[1][2] * w[2]) + sim.lightcone[i].vertex[1]) * sim.numpts) % sim.numpts;
					if (base_pos[1] < 0) base_pos[1] += sim.numpts;

					commdir[1] = phi->lattice().getRankDim1(base_pos[1]);
					j = commdir[1]*parallel.grid_size()[0];
					commdir[1] -= parallel.grid_rank()[1];

					if (commdir[1] < -1) commdir[1] += parallel.grid_size()[1];
					else if (commdir[1] > 1) commdir[1] -= parallel.grid_size()[1];

					base_pos[2] = (int) floor((maphdr.distance * (R[2][0] * w[0] + R[2][1] * w[1] + R[2][2] * w[2]) + sim.lightcone[i].vertex[2]) * sim.numpts) % sim.numpts;
					if (base_pos[2] < 0) base_pos[2] += sim.numpts;

					commdir[0] = phi->lattice().getRankDim0(base_pos[2]);
					j += commdir[0];
					commdir[0] -= parallel.grid_rank()[0];

					if (commdir[0] < -1) commdir[0] += parallel.grid_size()[0];
					else if (commdir[0] > 1) commdir[0] -= parallel.grid_size()[0];

					if ((io_group_size == 0 && parallel.rank() == ((shell - shell_inner) * parallel.size()) / (shell_outer + 1 - shell_inner)) || (io_group_size > 0 && shell - shell_inner == shell_write && ((pixbatch_delim[2].back() >= io_group_size && p / (pixbatch_delim[2].back() / io_group_size) < io_group_size && p / (pixbatch_delim[2].back() / io_group_size) == parallel.rank() - (shell_write * parallel.size() + shell_outer - shell_inner) / (shell_outer + 1 - shell_inner)) || (parallel.rank() - (shell_write * parallel.size() + shell_outer - shell_inner) / (shell_outer + 1 - shell_inner) == io_group_size - 1 && (pixbatch_delim[2].back() < io_group_size || p / (pixbatch_delim[2].back() / io_group_size) >= io_group_size))))) {
						sender_proc.push_back(j);
					}

					if (commdir[0] * commdir[0] > 1 || commdir[1] * commdir[1] > 1) continue;

					ring2nest64(maphdr.Nside_ring, p, &pix);
					pix *= pixbatch_size[0].back();

					if (p < pixbatch_delim[0].back()) pixbatch_type = 0;
					else if (p < pixbatch_delim[1].back()) pixbatch_type = 1;
					else pixbatch_type = 2;

					j = 3*commdir[0]+commdir[1]+4;

					if (pixbuf_size[j] + pixbatch_size[pixbatch_type].back() > pixbuf_reserve[j])
					{
						do
						{
							pixbuf_reserve[j] += PIXBUFFER;
						}
						while (pixbuf_size[j] + pixbatch_size[pixbatch_type].back() > pixbuf_reserve[j]);

						for (q = 0; q < LIGHTCONE_MAX_FIELDS; q++)
						{
							if (pixbuf[q][j] != NULL)
							{
								pixbuf[q][j] = (Real *) realloc((void *) pixbuf[q][j], sizeof(Real) * pixbuf_reserve[j]);
								if (pixbuf[q][j] == NULL)
								{
									cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": proc#" << parallel.rank() << " unable to allocate memory for pixelisation!" << endl;
									parallel.abortForce();
								}
							}
						}
					}

					for (q = 0; q < pixbatch_size[pixbatch_type].back(); pix++)
					{
						pixcount++;

						if (pixbatch_type)
						{
							nest2ring64(maphdr.Nside, pix, &pix2);
							if (pix2 >= maphdr.Npix) continue;
						}

						pix2vec_nest64(maphdr.Nside, pix, w);

						pos[0] = (maphdr.distance * (R[0][0] * w[0] + R[0][1] * w[1] + R[0][2] * w[2]) + sim.lightcone[i].vertex[0]) * sim.numpts;
						pos[1] = (maphdr.distance * (R[1][0] * w[0] + R[1][1] * w[1] + R[1][2] * w[2]) + sim.lightcone[i].vertex[1]) * sim.numpts;
						pos[2] = (maphdr.distance * (R[2][0] * w[0] + R[2][1] * w[1] + R[2][2] * w[2]) + sim.lightcone[i].vertex[2]) * sim.numpts;

						if (pos[0] >= 0)
						{
							w[0] = modf(pos[0], &temp);
							base_pos[0] = (int) temp % sim.numpts;
						}
						else
						{
							w[0] = 1. + modf(pos[0], &temp);
							base_pos[0] = sim.numpts - 1 - (((int) -temp) % sim.numpts);
						}
						if (pos[1] >= 0)
						{
							w[1] = modf(pos[1], &temp);
							base_pos[1] = (int) temp % sim.numpts;
						}
						else
						{
							w[1] = 1. + modf(pos[1], &temp);
							base_pos[1] = sim.numpts - 1 - (((int) -temp) % sim.numpts);
						}
						if (pos[2] >= 0)
						{
							w[2] = modf(pos[2], &temp);
							base_pos[2] = (int) temp % sim.numpts;
						}
						else
						{
							w[2] = 1. + modf(pos[2], &temp);
							base_pos[2] = sim.numpts - 1 - (((int) -temp) % sim.numpts);
						}

						if (xsim.setCoord(base_pos))
						{
							if (sim.out_lightcone[i] & MASK_PHI)
							{
								*(pixbuf[LIGHTCONE_PHI_OFFSET][j]+pixbuf_size[j]+q) = (1.-w[0]) * (1.-w[1]) * ((1.-w[2]) * (*phi)(xsim) + w[2] * (*phi)(xsim+2));
								*(pixbuf[LIGHTCONE_PHI_OFFSET][j]+pixbuf_size[j]+q) += w[0] * (1.-w[1]) * ((1.-w[2]) * (*phi)(xsim+0) + w[2] * (*phi)(xsim+0+2));
								*(pixbuf[LIGHTCONE_PHI_OFFSET][j]+pixbuf_size[j]+q) += w[0] * w[1] * ((1.-w[2]) * (*phi)(xsim+0+1) + w[2] * (*phi)(xsim+0+1+2));
								*(pixbuf[LIGHTCONE_PHI_OFFSET][j]+pixbuf_size[j]+q) += (1.-w[0]) * w[1] * ((1.-w[2]) * (*phi)(xsim+1) + w[2] * (*phi)(xsim+1+2));
							}
							if (sim.out_lightcone[i] & MASK_CHI)
							{
								*(pixbuf[LIGHTCONE_CHI_OFFSET][j]+pixbuf_size[j]+q) = (1.-w[0]) * (1.-w[1]) * ((1.-w[2]) * (*chi)(xsim) + w[2] * (*chi)(xsim+2));
								*(pixbuf[LIGHTCONE_CHI_OFFSET][j]+pixbuf_size[j]+q) += w[0] * (1.-w[1]) * ((1.-w[2]) * (*chi)(xsim+0) + w[2] * (*chi)(xsim+0+2));
								*(pixbuf[LIGHTCONE_CHI_OFFSET][j]+pixbuf_size[j]+q) += w[0] * w[1] * ((1.-w[2]) * (*chi)(xsim+0+1) + w[2] * (*chi)(xsim+0+1+2));
								*(pixbuf[LIGHTCONE_CHI_OFFSET][j]+pixbuf_size[j]+q) += (1.-w[0]) * w[1] * ((1.-w[2]) * (*chi)(xsim+1) + w[2] * (*chi)(xsim+1+2));
							}
							if (sim.out_lightcone[i] & MASK_B)
							{
#ifdef LIGHTCONE_INTERPOLATE
								*(pixbuf[LIGHTCONE_B_OFFSET][j]+pixbuf_size[j]+q) = (1.-w[2]) * (1.-w[1]) * ((1.-w[0]) * (*Bi)(xsim-0,0) + w[0] * (*Bi)(xsim+0,0) + (*Bi)(xsim,0));
								*(pixbuf[LIGHTCONE_B_OFFSET][j]+pixbuf_size[j]+q) += w[2] * (1.-w[1]) * ((1.-w[0]) * (*Bi)(xsim-0+2,0) + w[0] * (*Bi)(xsim+0+2,0) + (*Bi)(xsim+2,0));
								*(pixbuf[LIGHTCONE_B_OFFSET][j]+pixbuf_size[j]+q) += w[2] * w[1] * ((1.-w[0]) * (*Bi)(xsim-0+1+2,0) + w[0] * (*Bi)(xsim+0+1+2,0) + (*Bi)(xsim+1+2,0));
								*(pixbuf[LIGHTCONE_B_OFFSET][j]+pixbuf_size[j]+q) += (1.-w[2]) * w[1] * ((1.-w[0]) * (*Bi)(xsim-0+1,0) + w[0] * (*Bi)(xsim+0+1,0) + (*Bi)(xsim+1,0));

								*(pixbuf[LIGHTCONE_B_OFFSET+1][j]+pixbuf_size[j]+q) = (1.-w[2]) * (1.-w[0]) * ((1.-w[1]) * (*Bi)(xsim-1,1) + w[1] * (*Bi)(xsim+1,1) + (*Bi)(xsim,1));
								*(pixbuf[LIGHTCONE_B_OFFSET+1][j]+pixbuf_size[j]+q) += w[2] * (1.-w[0]) * ((1.-w[1]) * (*Bi)(xsim-1+2,1) + w[1] * (*Bi)(xsim+1+2,1) + (*Bi)(xsim+2,1));
								*(pixbuf[LIGHTCONE_B_OFFSET+1][j]+pixbuf_size[j]+q) += w[2] * w[0] * ((1.-w[1]) * (*Bi)(xsim+0-1+2,1) + w[1] * (*Bi)(xsim+0+1+2,1) + (*Bi)(xsim+0+2,1));
								*(pixbuf[LIGHTCONE_B_OFFSET+1][j]+pixbuf_size[j]+q) += (1.-w[2]) * w[0] * ((1.-w[1]) * (*Bi)(xsim+0-1,1) + w[1] * (*Bi)(xsim+0+1,1) + (*Bi)(xsim+0,1));

								*(pixbuf[LIGHTCONE_B_OFFSET+2][j]+pixbuf_size[j]+q) = (1.-w[0]) * (1.-w[1]) * ((1.-w[2]) * (*Bi)(xsim-2,2) + w[2] * (*Bi)(xsim+2,2) + (*Bi)(xsim,2));
								*(pixbuf[LIGHTCONE_B_OFFSET+2][j]+pixbuf_size[j]+q) += w[0] * (1.-w[1]) * ((1.-w[2]) * (*Bi)(xsim+0-2,2) + w[2] * (*Bi)(xsim+0+2,2) + (*Bi)(xsim+0,2));
								*(pixbuf[LIGHTCONE_B_OFFSET+2][j]+pixbuf_size[j]+q) += w[0] * w[1] * ((1.-w[2]) * (*Bi)(xsim+0+1-2,2) + w[2] * (*Bi)(xsim+0+1+2,2) + (*Bi)(xsim+0+1,2));
								*(pixbuf[LIGHTCONE_B_OFFSET+2][j]+pixbuf_size[j]+q) += (1.-w[0]) * w[1] * ((1.-w[2]) * (*Bi)(xsim+1-2,2) + w[2] * (*Bi)(xsim+1+2,2) + (*Bi)(xsim+1,2));

								*(pixbuf[LIGHTCONE_B_OFFSET][j]+pixbuf_size[j]+q) /= 2. * a * a * sim.numpts;
								*(pixbuf[LIGHTCONE_B_OFFSET+1][j]+pixbuf_size[j]+q) /= 2. * a * a * sim.numpts;
								*(pixbuf[LIGHTCONE_B_OFFSET+2][j]+pixbuf_size[j]+q) /= 2. * a * a * sim.numpts;
#else
								if (w[0] > 0.5)
								{
									*(pixbuf[LIGHTCONE_B_OFFSET][j]+pixbuf_size[j]+q) = (1.5-w[0]) * (1.-w[1]) * ((1.-w[2]) * (*Bi)(xsim,0) + w[2] * (*Bi)(xsim+2,0));
									*(pixbuf[LIGHTCONE_B_OFFSET][j]+pixbuf_size[j]+q) += (1.5-w[0]) * w[1] * ((1.-w[2]) * (*Bi)(xsim+1,0) + w[2] * (*Bi)(xsim+1+2,0));
									*(pixbuf[LIGHTCONE_B_OFFSET][j]+pixbuf_size[j]+q) += (w[0]-0.5) * w[1] * ((1.-w[2]) * (*Bi)(xsim+0+1,0) + w[2] * (*Bi)(xsim+0+1+2,0));
									*(pixbuf[LIGHTCONE_B_OFFSET][j]+pixbuf_size[j]+q) += (w[0]-0.5) * (1.-w[1]) * ((1.-w[2]) * (*Bi)(xsim+0,0) + w[2] * (*Bi)(xsim+0+2,0));
								}
								else
								{
									*(pixbuf[LIGHTCONE_B_OFFSET][j]+pixbuf_size[j]+q) = (0.5+w[0]) * (1.-w[1]) * ((1.-w[2]) * (*Bi)(xsim,0) + w[2] * (*Bi)(xsim+2,0));
									*(pixbuf[LIGHTCONE_B_OFFSET][j]+pixbuf_size[j]+q) += (0.5+w[0]) * w[1] * ((1.-w[2]) * (*Bi)(xsim+1,0) + w[2] * (*Bi)(xsim+1+2,0));
									*(pixbuf[LIGHTCONE_B_OFFSET][j]+pixbuf_size[j]+q) += (0.5-w[0]) * w[1] * ((1.-w[2]) * (*Bi)(xsim-0+1,0) + w[2] * (*Bi)(xsim-0+1+2,0));
									*(pixbuf[LIGHTCONE_B_OFFSET][j]+pixbuf_size[j]+q) += (0.5-w[0]) * (1.-w[1]) * ((1.-w[2]) * (*Bi)(xsim-0,0) + w[2] * (*Bi)(xsim-0+2,0));
								}
								if (w[1] > 0.5)
								{
									*(pixbuf[LIGHTCONE_B_OFFSET+1][j]+pixbuf_size[j]+q) = (1.-w[0]) * (1.5-w[1]) * ((1.-w[2]) * (*Bi)(xsim,1) + w[2] * (*Bi)(xsim+2,1));
									*(pixbuf[LIGHTCONE_B_OFFSET+1][j]+pixbuf_size[j]+q) += (1.-w[0]) * (w[1]-0.5) * ((1.-w[2]) * (*Bi)(xsim+1,1) + w[2] * (*Bi)(xsim+1+2,1));
									*(pixbuf[LIGHTCONE_B_OFFSET+1][j]+pixbuf_size[j]+q) += w[0] * (w[1]-0.5) * ((1.-w[2]) * (*Bi)(xsim+0+1,1) + w[2] * (*Bi)(xsim+0+1+2,1));
									*(pixbuf[LIGHTCONE_B_OFFSET+1][j]+pixbuf_size[j]+q) += w[0] * (1.5-w[1]) * ((1.-w[2]) * (*Bi)(xsim+0,1) + w[2] * (*Bi)(xsim+0+2,1));
								}
								else
								{
									*(pixbuf[LIGHTCONE_B_OFFSET+1][j]+pixbuf_size[j]+q) = (1.-w[0]) * (0.5+w[1]) * ((1.-w[2]) * (*Bi)(xsim,1) + w[2] * (*Bi)(xsim+2,1));
									*(pixbuf[LIGHTCONE_B_OFFSET+1][j]+pixbuf_size[j]+q) += (1.-w[0]) * (0.5-w[1]) * ((1.-w[2]) * (*Bi)(xsim-1,1) + w[2] * (*Bi)(xsim-1+2,1));
									*(pixbuf[LIGHTCONE_B_OFFSET+1][j]+pixbuf_size[j]+q) += w[0] * (0.5-w[1]) * ((1.-w[2]) * (*Bi)(xsim+0-1,1) + w[2] * (*Bi)(xsim+0-1+2,1));
									*(pixbuf[LIGHTCONE_B_OFFSET+1][j]+pixbuf_size[j]+q) += w[0] * (0.5+w[1]) * ((1.-w[2]) * (*Bi)(xsim+0,1) + w[2] * (*Bi)(xsim+0+2,1));
								}
								if (w[2] > 0.5)
								{
									*(pixbuf[LIGHTCONE_B_OFFSET+2][j]+pixbuf_size[j]+q) = (1.-w[0]) * (1.-w[1]) * ((1.5-w[2]) * (*Bi)(xsim,2) + (w[2]-0.5) * (*Bi)(xsim+2,2));
									*(pixbuf[LIGHTCONE_B_OFFSET+2][j]+pixbuf_size[j]+q) += (1.-w[0]) * w[1] * ((1.5-w[2]) * (*Bi)(xsim+1,2) + (w[2]-0.5) * (*Bi)(xsim+1+2,2));
									*(pixbuf[LIGHTCONE_B_OFFSET+2][j]+pixbuf_size[j]+q) += w[0] * w[1] * ((1.5-w[2]) * (*Bi)(xsim+0+1,2) + (w[2]-0.5) * (*Bi)(xsim+0+1+2,2));
									*(pixbuf[LIGHTCONE_B_OFFSET+2][j]+pixbuf_size[j]+q) += w[0] * (1.-w[1]) * ((1.5-w[2]) * (*Bi)(xsim+0,2) + (w[2]-0.5) * (*Bi)(xsim+0+2,2));
								}
								else
								{
									*(pixbuf[LIGHTCONE_B_OFFSET+2][j]+pixbuf_size[j]+q) = (1.-w[0]) * (1.-w[1]) * ((0.5+w[2]) * (*Bi)(xsim,2) + (0.5-w[2]) * (*Bi)(xsim-2,2));
									*(pixbuf[LIGHTCONE_B_OFFSET+2][j]+pixbuf_size[j]+q) += (1.-w[0]) * w[1] * ((0.5+w[2]) * (*Bi)(xsim+1,2) + (0.5-w[2]) * (*Bi)(xsim+1-2,2));
									*(pixbuf[LIGHTCONE_B_OFFSET+2][j]+pixbuf_size[j]+q) += w[0] * w[1] * ((0.5+w[2]) * (*Bi)(xsim+0+1,2) + (0.5-w[2]) * (*Bi)(xsim+0+1-2,2));
									*(pixbuf[LIGHTCONE_B_OFFSET+2][j]+pixbuf_size[j]+q) += w[0] * (1.-w[1]) * ((0.5+w[2]) * (*Bi)(xsim+0,2) + (0.5-w[2]) * (*Bi)(xsim+0-2,2));
								}
								*(pixbuf[LIGHTCONE_B_OFFSET][j]+pixbuf_size[j]+q) /= a * a * sim.numpts;
								*(pixbuf[LIGHTCONE_B_OFFSET+1][j]+pixbuf_size[j]+q) /= a * a * sim.numpts;
								*(pixbuf[LIGHTCONE_B_OFFSET+2][j]+pixbuf_size[j]+q) /= a * a * sim.numpts;
#endif
							}
							if (sim.out_lightcone[i] & MASK_HIJ)
							{
								*(pixbuf[LIGHTCONE_HIJ_OFFSET][j]+pixbuf_size[j]+q) = (1.-w[0]) * (1.-w[1]) * ((1.-w[2]) * (*Sij)(xsim,0,0) + w[2] * (*Sij)(xsim+2,0,0));
								*(pixbuf[LIGHTCONE_HIJ_OFFSET][j]+pixbuf_size[j]+q) += w[0] * (1.-w[1]) * ((1.-w[2]) * (*Sij)(xsim+0,0,0) + w[2] * (*Sij)(xsim+0+2,0,0));
								*(pixbuf[LIGHTCONE_HIJ_OFFSET][j]+pixbuf_size[j]+q) += w[0] * w[1] * ((1.-w[2]) * (*Sij)(xsim+0+1,0,0) + w[2] * (*Sij)(xsim+0+1+2,0,0));
								*(pixbuf[LIGHTCONE_HIJ_OFFSET][j]+pixbuf_size[j]+q) += (1.-w[0]) * w[1] * ((1.-w[2]) * (*Sij)(xsim+1,0,0) + w[2] * (*Sij)(xsim+1+2,0,0));

								*(pixbuf[LIGHTCONE_HIJ_OFFSET+1][j]+pixbuf_size[j]+q) = (1.-w[2]) * 0.25 * ((*Sij)(xsim,0,1) + (1.-w[0]) * ((*Sij)(xsim-0,0,1) + (1.-w[1]) * (*Sij)(xsim-0-1,0,1) + w[1] * (*Sij)(xsim-0+1,0,1)) + w[0] * ((*Sij)(xsim+0,0,1) + (1.-w[1]) * (*Sij)(xsim+0-1,0,1) + w[1] * (*Sij)(xsim+0+1,0,1)) + (1.-w[1]) * (*Sij)(xsim-1,0,1) + w[1] * (*Sij)(xsim+1,0,1));
								*(pixbuf[LIGHTCONE_HIJ_OFFSET+1][j]+pixbuf_size[j]+q) += w[2] * 0.25 * ((*Sij)(xsim+2,0,1) + (1.-w[0]) * ((*Sij)(xsim-0+2,0,1) + (1.-w[1]) * (*Sij)(xsim-0-1+2,0,1) + w[1] * (*Sij)(xsim-0+1+2,0,1)) + w[0] * ((*Sij)(xsim+0+2,0,1) + (1.-w[1]) * (*Sij)(xsim+0-1+2,0,1) + w[1] * (*Sij)(xsim+0+1+2,0,1)) + (1.-w[1]) * (*Sij)(xsim-1+2,0,1) + w[1] * (*Sij)(xsim+1+2,0,1));

								*(pixbuf[LIGHTCONE_HIJ_OFFSET+2][j]+pixbuf_size[j]+q) = (1.-w[1]) * 0.25 * ((*Sij)(xsim,0,2) + (1.-w[0]) * ((*Sij)(xsim-0,0,2) + (1.-w[2]) * (*Sij)(xsim-0-2,0,2) + w[2] * (*Sij)(xsim-0+2,0,2)) + w[0] * ((*Sij)(xsim+0,0,2) + (1.-w[2]) * (*Sij)(xsim+0-2,0,2) + w[2] * (*Sij)(xsim+0+2,0,2)) + (1.-w[2]) * (*Sij)(xsim-2,0,2) + w[2] * (*Sij)(xsim+2,0,2));
								*(pixbuf[LIGHTCONE_HIJ_OFFSET+2][j]+pixbuf_size[j]+q) += w[1] * 0.25 * ((*Sij)(xsim+1,0,2) + (1.-w[0]) * ((*Sij)(xsim-0+1,0,2) + (1.-w[2]) * (*Sij)(xsim-0+1-2,0,2) + w[2] * (*Sij)(xsim-0+1+2,0,2)) + w[0] * ((*Sij)(xsim+0+1,0,2) + (1.-w[2]) * (*Sij)(xsim+0+1-2,0,2) + w[2] * (*Sij)(xsim+0+1+2,0,2)) + (1.-w[2]) * (*Sij)(xsim+1-2,0,2) + w[2] * (*Sij)(xsim+1+2,0,2));


								*(pixbuf[LIGHTCONE_HIJ_OFFSET+3][j]+pixbuf_size[j]+q) = (1.-w[0]) * (1.-w[1]) * ((1.-w[2]) * (*Sij)(xsim,1,1) + w[2] * (*Sij)(xsim+2,1,1));
								*(pixbuf[LIGHTCONE_HIJ_OFFSET+3][j]+pixbuf_size[j]+q) += w[0] * (1.-w[1]) * ((1.-w[2]) * (*Sij)(xsim+0,1,1) + w[2] * (*Sij)(xsim+0+2,1,1));
								*(pixbuf[LIGHTCONE_HIJ_OFFSET+3][j]+pixbuf_size[j]+q) += w[0] * w[1] * ((1.-w[2]) * (*Sij)(xsim+0+1,1,1) + w[2] * (*Sij)(xsim+0+1+2,1,1));
								*(pixbuf[LIGHTCONE_HIJ_OFFSET+3][j]+pixbuf_size[j]+q) += (1.-w[0]) * w[1] * ((1.-w[2]) * (*Sij)(xsim+1,1,1) + w[2] * (*Sij)(xsim+1+2,1,1));

								*(pixbuf[LIGHTCONE_HIJ_OFFSET+4][j]+pixbuf_size[j]+q) = (1.-w[0]) * 0.25 * ((*Sij)(xsim,1,2) + (1.-w[1]) * ((*Sij)(xsim-1,1,2) + (1.-w[2]) * (*Sij)(xsim-1-2,1,2) + w[2] * (*Sij)(xsim-1+2,1,2)) + w[1] * ((*Sij)(xsim+1,1,2) + (1.-w[2]) * (*Sij)(xsim+1-2,1,2) + w[2] * (*Sij)(xsim+1+2,1,2)) + (1.-w[2]) * (*Sij)(xsim-2,1,2) + w[2] * (*Sij)(xsim+2,1,2));
								*(pixbuf[LIGHTCONE_HIJ_OFFSET+4][j]+pixbuf_size[j]+q) += w[0] * 0.25 * ((*Sij)(xsim+0,1,2) + (1.-w[1]) * ((*Sij)(xsim+0-1,1,2) + (1.-w[2]) * (*Sij)(xsim+0-1-2,1,2) + w[2] * (*Sij)(xsim+0-1+2,1,2)) + w[1] * ((*Sij)(xsim+0+1,1,2) + (1.-w[2]) * (*Sij)(xsim+0+1-2,1,2) + w[2] * (*Sij)(xsim+0+1+2,1,2)) + (1.-w[2]) * (*Sij)(xsim+0-2,1,2) + w[2] * (*Sij)(xsim+0+2,1,2));
							}
						}
						else
						{
							if (sim.out_lightcone[i] & MASK_PHI)
								*(pixbuf[LIGHTCONE_PHI_OFFSET][j]+pixbuf_size[j]+q) = 0;
							if (sim.out_lightcone[i] & MASK_CHI)
								*(pixbuf[LIGHTCONE_CHI_OFFSET][j]+pixbuf_size[j]+q) = 0;
							if (sim.out_lightcone[i] & MASK_B)
							{
								*(pixbuf[LIGHTCONE_B_OFFSET][j]+pixbuf_size[j]+q) = 0;
								*(pixbuf[LIGHTCONE_B_OFFSET+1][j]+pixbuf_size[j]+q) = 0;
								*(pixbuf[LIGHTCONE_B_OFFSET+2][j]+pixbuf_size[j]+q) = 0;
							}
							if (sim.out_lightcone[i] & MASK_HIJ)
							{
								*(pixbuf[LIGHTCONE_HIJ_OFFSET][j]+pixbuf_size[j]+q) = 0;
								*(pixbuf[LIGHTCONE_HIJ_OFFSET+1][j]+pixbuf_size[j]+q) = 0;
								*(pixbuf[LIGHTCONE_HIJ_OFFSET+2][j]+pixbuf_size[j]+q) = 0;
								*(pixbuf[LIGHTCONE_HIJ_OFFSET+3][j]+pixbuf_size[j]+q) = 0;
								*(pixbuf[LIGHTCONE_HIJ_OFFSET+4][j]+pixbuf_size[j]+q) = 0;
							}
						}

						q++;
					} // q-loop

					pixbuf_size[j] += pixbatch_size[pixbatch_type].back();

					if (j == 4)
					{
						pixbatch_id.push_back(p);
					}
				} // p-loop

				p = 0;
				for (j = 0; j < 3; j++)
				{
					if (pixbuf_size[3*j]+pixbuf_size[3*j+1]+pixbuf_size[3*j+2] > p) p = pixbuf_size[3*j]+pixbuf_size[3*j+1]+pixbuf_size[3*j+2];
				}

				if (p > 0)
				{
					commbuf = (Real *) malloc(sizeof(Real) * p);

					for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
					{
						if (pixbuf[j][0] != NULL)
						{
							if (parallel.grid_rank()[0] % 2 == 0)
							{
								if (pixbuf_size[0]+pixbuf_size[1]+pixbuf_size[2] > 0)
								{
									if (pixbuf_size[0] > 0)
										memcpy((void *) commbuf, (void *) pixbuf[j][0], pixbuf_size[0] * sizeof(Real));
									if (pixbuf_size[1] > 0)
										memcpy((void *) (commbuf+pixbuf_size[0]), (void *) pixbuf[j][1], pixbuf_size[1] * sizeof(Real));
									if (pixbuf_size[2] > 0)
										memcpy((void *) (commbuf+pixbuf_size[0]+pixbuf_size[1]), (void *) pixbuf[j][2], pixbuf_size[2] * sizeof(Real));
									parallel.send_dim0<Real>(commbuf, pixbuf_size[0]+pixbuf_size[1]+pixbuf_size[2], (parallel.grid_size()[0]+parallel.grid_rank()[0]-1) % parallel.grid_size()[0]);
								}
								if (pixbuf_size[3]+pixbuf_size[4]+pixbuf_size[5] > 0)
								{
									parallel.receive_dim0<Real>(commbuf, pixbuf_size[3]+pixbuf_size[4]+pixbuf_size[5], (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
									for (q = 0; q < pixbuf_size[3]; q++)
										*(pixbuf[j][3]+q) += commbuf[q];
									for (q = 0; q < pixbuf_size[4]; q++)
										*(pixbuf[j][4]+q) += commbuf[pixbuf_size[3]+q];
									for (q = 0; q < pixbuf_size[5]; q++)
										*(pixbuf[j][5]+q) += commbuf[pixbuf_size[3]+pixbuf_size[4]+q];
								}
							}
							else
							{
								if (pixbuf_size[3]+pixbuf_size[4]+pixbuf_size[5] > 0 && parallel.grid_size()[0] > 2)
								{
									parallel.receive_dim0<Real>(commbuf, pixbuf_size[3]+pixbuf_size[4]+pixbuf_size[5], (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
									for (q = 0; q < pixbuf_size[3]; q++)
										*(pixbuf[j][3]+q) += commbuf[q];
									for (q = 0; q < pixbuf_size[4]; q++)
										*(pixbuf[j][4]+q) += commbuf[pixbuf_size[3]+q];
									for (q = 0; q < pixbuf_size[5]; q++)
										*(pixbuf[j][5]+q) += commbuf[pixbuf_size[3]+pixbuf_size[4]+q];
								}
								if (pixbuf_size[0]+pixbuf_size[1]+pixbuf_size[2] > 0)
								{
									if (pixbuf_size[0] > 0)
										memcpy((void *) commbuf, (void *) pixbuf[j][0], pixbuf_size[0] * sizeof(Real));
									if (pixbuf_size[1] > 0)
										memcpy((void *) (commbuf+pixbuf_size[0]), (void *) pixbuf[j][1], pixbuf_size[1] * sizeof(Real));
									if (pixbuf_size[2] > 0)
										memcpy((void *) (commbuf+pixbuf_size[0]+pixbuf_size[1]), (void *) pixbuf[j][2], pixbuf_size[2] * sizeof(Real));
									parallel.send_dim0<Real>(commbuf, pixbuf_size[0]+pixbuf_size[1]+pixbuf_size[2], (parallel.grid_size()[0]+parallel.grid_rank()[0]-1) % parallel.grid_size()[0]);
								}
							}

							if (parallel.grid_rank()[0] % 2 == 0)
							{
								if (pixbuf_size[6]+pixbuf_size[7]+pixbuf_size[8] > 0)
								{
									if (pixbuf_size[6] > 0)
										memcpy((void *) commbuf, (void *) pixbuf[j][6], pixbuf_size[6] * sizeof(Real));
									if (pixbuf_size[7] > 0)
										memcpy((void *) (commbuf+pixbuf_size[6]), (void *) pixbuf[j][7], pixbuf_size[7] * sizeof(Real));
									if (pixbuf_size[8] > 0)
										memcpy((void *) (commbuf+pixbuf_size[6]+pixbuf_size[7]), (void *) pixbuf[j][8], pixbuf_size[8] * sizeof(Real));
									parallel.send_dim0<Real>(commbuf, pixbuf_size[6]+pixbuf_size[7]+pixbuf_size[8], (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
								}
								if (pixbuf_size[3]+pixbuf_size[4]+pixbuf_size[5] > 0 && parallel.grid_size()[0] > 2)
								{
									parallel.receive_dim0<Real>(commbuf, pixbuf_size[3]+pixbuf_size[4]+pixbuf_size[5], (parallel.grid_size()[0]+parallel.grid_rank()[0]-1) % parallel.grid_size()[0]);
									for (q = 0; q < pixbuf_size[3]; q++)
										*(pixbuf[j][3]+q) += commbuf[q];
									for (q = 0; q < pixbuf_size[4]; q++)
										*(pixbuf[j][4]+q) += commbuf[pixbuf_size[3]+q];
									for (q = 0; q < pixbuf_size[5]; q++)
										*(pixbuf[j][5]+q) += commbuf[pixbuf_size[3]+pixbuf_size[4]+q];
								}
							}
							else
							{
								if (pixbuf_size[3]+pixbuf_size[4]+pixbuf_size[5] > 0)
								{
									parallel.receive_dim0<Real>(commbuf, pixbuf_size[3]+pixbuf_size[4]+pixbuf_size[5], (parallel.grid_size()[0]+parallel.grid_rank()[0]-1) % parallel.grid_size()[0]);
									for (q = 0; q < pixbuf_size[3]; q++)
										*(pixbuf[j][3]+q) += commbuf[q];
									for (q = 0; q < pixbuf_size[4]; q++)
										*(pixbuf[j][4]+q) += commbuf[pixbuf_size[3]+q];
									for (q = 0; q < pixbuf_size[5]; q++)
										*(pixbuf[j][5]+q) += commbuf[pixbuf_size[3]+pixbuf_size[4]+q];
								}
								if (pixbuf_size[6]+pixbuf_size[7]+pixbuf_size[8] > 0)
								{
									if (pixbuf_size[6] > 0)
										memcpy((void *) commbuf, (void *) pixbuf[j][6], pixbuf_size[6] * sizeof(Real));
									if (pixbuf_size[7] > 0)
										memcpy((void *) (commbuf+pixbuf_size[6]), (void *) pixbuf[j][7], pixbuf_size[7] * sizeof(Real));
									if (pixbuf_size[8] > 0)
										memcpy((void *) (commbuf+pixbuf_size[6]+pixbuf_size[7]), (void *) pixbuf[j][8], pixbuf_size[8] * sizeof(Real));
									parallel.send_dim0<Real>(commbuf, pixbuf_size[6]+pixbuf_size[7]+pixbuf_size[8], (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
								}
							}

							if (parallel.grid_rank()[1] % 2 == 0)
							{
								if (pixbuf_size[3] > 0)
									parallel.send_dim1<Real>(pixbuf[j][3], pixbuf_size[3], (parallel.grid_size()[1]+parallel.grid_rank()[1]-1) % parallel.grid_size()[1]);
								if (pixbuf_size[4] > 0)
								{
									parallel.receive_dim1<Real>(commbuf, pixbuf_size[4], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
									for (q = 0; q < pixbuf_size[4]; q++)
										*(pixbuf[j][4]+q) += commbuf[q];
								}
							}
							else
							{
								if (pixbuf_size[4] > 0 && parallel.grid_size()[1] > 2)
								{
									parallel.receive_dim1<Real>(commbuf, pixbuf_size[4], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
									for (q = 0; q < pixbuf_size[4]; q++)
										*(pixbuf[j][4]+q) += commbuf[q];
								}
								if (pixbuf_size[3] > 0)
									parallel.send_dim1<Real>(pixbuf[j][3], pixbuf_size[3], (parallel.grid_size()[1]+parallel.grid_rank()[1]-1) % parallel.grid_size()[1]);
							}

							if (parallel.grid_rank()[1] % 2 == 0)
							{
								if (pixbuf_size[5] > 0)
									parallel.send_dim1<Real>(pixbuf[j][5], pixbuf_size[5], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
								if (pixbuf_size[4] > 0 && parallel.grid_size()[1] > 2)
								{
									parallel.receive_dim1<Real>(commbuf, pixbuf_size[4], (parallel.grid_size()[1]+parallel.grid_rank()[1]-1) % parallel.grid_size()[1]);
									for (q = 0; q < pixbuf_size[4]; q++)
										*(pixbuf[j][4]+q) += commbuf[q];
								}
							}
							else
							{
								if (pixbuf_size[4] > 0)
								{
									parallel.receive_dim1<Real>(commbuf, pixbuf_size[4], (parallel.grid_size()[1]+parallel.grid_rank()[1]-1) % parallel.grid_size()[1]);
									for (q = 0; q < pixbuf_size[4]; q++)
										*(pixbuf[j][4]+q) += commbuf[q];
								}
								if (pixbuf_size[5] > 0)
									parallel.send_dim1<Real>(pixbuf[j][5], pixbuf_size[5], (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
							}
						}
					}

					free(commbuf);
				}

				if (io_group_size == 0 && parallel.rank() == ((shell - shell_inner) * parallel.size() / (shell_outer + 1 - shell_inner)))
				{
					for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
					{
						if (pixbuf[j][0] != NULL)
						{
							if (bytes2 == 0)
							{
								outbuf[j] = (char *) malloc(maphdr.Npix * maphdr.precision + 272);

								if (outbuf[j] == NULL)
								{
									cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": proc#" << parallel.rank() << " unable to allocate " << maphdr.Npix * maphdr.precision + 272 << " bytes of memory for pixelisation!" << endl;
									parallel.abortForce();
								}
							}
							else
							{
								outbuf[j] = (char *) realloc((void *) outbuf[j], bytes2 + maphdr.Npix * maphdr.precision + 272);

								if (outbuf[j] == NULL)
								{
									cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": proc#" << parallel.rank() << " unable to reallocate " << bytes2 + maphdr.Npix * maphdr.precision + 272 << " bytes of memory (" << maphdr.Npix * maphdr.precision + 272 << " additional bytes) for pixelisation!" << endl;
									parallel.abortForce();
								}
							}

							blocksize = 256;
							memcpy((void *) (outbuf[j] + bytes2), (void *) &blocksize, 4);
							memcpy((void *) (outbuf[j] + bytes2 + 4), (void *) &maphdr, 256);
							memcpy((void *) (outbuf[j] + bytes2 + 260), (void *) &blocksize, 4);
							blocksize = maphdr.precision * maphdr.Npix;
							memcpy((void *) (outbuf[j] + bytes2 + 264), (void *) &blocksize, 4);
							memcpy((void *) (outbuf[j] + bytes2 + 268 + blocksize), (void *) &blocksize, 4);
						}
					}
					offset2 = bytes2 + 268;
					bytes2 += maphdr.Npix * maphdr.precision + 272;
					p = 0;
					q = pixbatch_delim[2].back();
				}
				else if (io_group_size > 0 && shell - shell_inner == shell_write)
				{
					q = pixbatch_delim[2].back() / io_group_size;
					p = parallel.rank() - (shell_write * parallel.size() + shell_outer - shell_inner) / (shell_outer + 1 - shell_inner);

					for (j = 0; p * q >= pixbatch_delim[j].back(); j++);

					if ((p+1) * q >= pixbatch_delim[j].back() && j < 2)
					{
						bytes2 = (pixbatch_delim[j].back() - p * q) * pixbatch_size[j].back() * maphdr.precision;
						if ((p+1) * q >= pixbatch_delim[j+1].back() && j < 1)
						{
							bytes2 += (pixbatch_delim[j+1].back() - pixbatch_delim[j].back()) * pixbatch_size[j+1].back() * maphdr.precision;
							bytes2 += ((p+1) * q - pixbatch_delim[j+1].back()) * pixbatch_size[j+2].back() * maphdr.precision;
						}
						else
							bytes2 += ((p+1) * q - pixbatch_delim[j].back()) * pixbatch_size[j+1].back() * maphdr.precision;
					}
					else
						bytes2 = q * pixbatch_size[j].back() * maphdr.precision;

					if (p == 0)
					{
						bytes2 += 268;
						offset2 = 268;
					}
					else
						offset2 = 0;

					if (p == io_group_size-1)
					{
						bytes2 += 4;
						q = pixbatch_delim[2].back() % io_group_size;
						if (pixbatch_delim[2].back()-pixbatch_delim[1].back() < q)
						{
							bytes2 += (pixbatch_delim[2].back()-pixbatch_delim[1].back()) * pixbatch_size[2].back() * maphdr.precision;
							if (pixbatch_delim[2].back()-pixbatch_delim[0].back() < q)
							{
								bytes2 += (pixbatch_delim[1].back()-pixbatch_delim[0].back()) * pixbatch_size[1].back() * maphdr.precision;
								bytes2 += (q-pixbatch_delim[2].back()+pixbatch_delim[0].back()) * pixbatch_size[0].back() * maphdr.precision;
							}
							else
								bytes2 += (q-pixbatch_delim[2].back()+pixbatch_delim[1].back()) * pixbatch_size[1].back() * maphdr.precision;
						}
						else
							bytes2 += q * pixbatch_size[2].back() * maphdr.precision;

						q += pixbatch_delim[2].back() / io_group_size;
					}

					for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
					{
						if (pixbuf[j][0] != NULL)
						{
							if (bytes2 > 0)
							{
								outbuf[j] = (char *) malloc(bytes2);

								if (outbuf[j] == NULL)
								{
									cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": proc#" << parallel.rank() << " unable to allocate " << bytes2 << " bytes of memory for pixelisation!" << endl;
									parallel.abortForce();
								}
							}

							if (p == 0)
							{
								blocksize = 256;
								memcpy((void *) outbuf[j], (void *) &blocksize, 4);
								memcpy((void *) (outbuf[j] + 4), (void *) &maphdr, 256);
								memcpy((void *) (outbuf[j] + 260), (void *) &blocksize, 4);
								blocksize = maphdr.precision * maphdr.Npix;
								memcpy((void *) (outbuf[j] + 264), (void *) &blocksize, 4);
							}

							if (p == io_group_size-1)
							{
								blocksize = maphdr.precision * maphdr.Npix;
								memcpy((void *) (outbuf[j] + bytes2 - 4), (void *) &blocksize, 4);
							}
						}
					}

					p *= pixbatch_delim[2].back() / io_group_size;
				}

				pix = 0;
				pix2 = 0;

				if ((io_group_size == 0 && parallel.rank() == ((shell - shell_inner) * parallel.size()) / (shell_outer + 1 - shell_inner)) || (io_group_size > 0 && shell - shell_inner == shell_write))
				{
					if (q != sender_proc.size())
					{
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": proc#" << parallel.rank() << " pixel batch count mismatch! expecting " << q << " but sender list contains " << sender_proc.size() << " entries!" << endl;
						exit(-99);
					}

					for (int64_t p2 = p; p2 < p+q; p2 += n)
					{
						while (pix < pixbatch_id.size() && pixbatch_id[pix] < p2)
						{
							for (pixbatch_type = 0; pixbatch_delim[pixbatch_type].back() <= pixbatch_id[pix]; pixbatch_type++);
							if (io_group_size > 0 && pixbatch_delim[2].back() >= io_group_size && pixbatch_id[pix] / (pixbatch_delim[2].back() / io_group_size) < io_group_size)
							{
								for (n = 1; pix+n < pixbatch_id.size() && pixbatch_id[pix+n] == pixbatch_id[pix+n-1]+1 && pixbatch_id[pix+n] < pixbatch_delim[pixbatch_type].back() && (pixbatch_id[pix+n] / (pixbatch_delim[2].back() / io_group_size) == pixbatch_id[pix] / (pixbatch_delim[2].back() / io_group_size) || pixbatch_id[pix] / (pixbatch_delim[2].back() / io_group_size) == io_group_size-1); n++);
								for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
								{
									if (pixbuf[j][4] != NULL && pixbatch_size[pixbatch_type].back() > 0)
										parallel.send<Real>(pixbuf[j][4]+pix2, n*pixbatch_size[pixbatch_type].back(), (pixbatch_id[pix] / (pixbatch_delim[2].back() / io_group_size)) + ((shell - shell_inner) * parallel.size() + shell_outer - shell_inner) / (shell_outer + 1 - shell_inner));
								}
							}
							else
							{
								for (n = 1; pix+n < pixbatch_id.size() && pixbatch_id[pix+n] == pixbatch_id[pix+n-1]+1 && pixbatch_id[pix+n] < pixbatch_delim[pixbatch_type].back(); n++);
								for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
								{
									if (pixbuf[j][4] != NULL && pixbatch_size[pixbatch_type].back() > 0)
										parallel.send<Real>(pixbuf[j][4]+pix2, n*pixbatch_size[pixbatch_type].back(), (io_group_size ? io_group_size - 1 : 0) + ((shell - shell_inner) * parallel.size() + (io_group_size ? shell_outer - shell_inner : 0)) / (shell_outer + 1 - shell_inner));
								}
							}
							pix += n;
							pix2 += n*pixbatch_size[pixbatch_type].back();
						}

						for (pixbatch_type = 0; pixbatch_delim[pixbatch_type].back() <= p2; pixbatch_type++);

						for (n = 1; p2+n < p+q && sender_proc[p2+n-p] == sender_proc[p2-p] && p2+n < pixbatch_delim[pixbatch_type].back(); n++);

						if (sender_proc[p2-p] == parallel.rank())
						{
							if (pix+n-1 >= pixbatch_id.size())
							{
								cerr << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": proc#" << parallel.rank() << " pixel batch index mismatch! expecting " << p2 << " but ID list contains not enough elements!" << endl;
								exit(-99);
							}
							else if (pixbatch_id[pix] != p2)
							{
								cerr << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": proc#" << parallel.rank() << " pixel batch index mismatch! expecting " << p2 << " but ID list says " << pixbatch_id[pix] << "!" << endl;
								exit(-99);
							}
							for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
							{
								if (pixbuf[j][4] != NULL && pixbatch_size[pixbatch_type].back() > 0)
									memcpy((void *) (outbuf[j]+offset2), (void *) (pixbuf[j][4]+pix2), n*pixbatch_size[pixbatch_type].back()*maphdr.precision);
							}
							pix += n;
							pix2 += n*pixbatch_size[pixbatch_type].back();
						}
						else
						{
							for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
							{
								if (outbuf[j] != NULL && pixbatch_size[pixbatch_type].back() > 0)
									parallel.receive<Real>((Real *) (outbuf[j]+offset2), n*pixbatch_size[pixbatch_type].back(), sender_proc[p2-p]);
							}
						}

						offset2 += n*pixbatch_size[pixbatch_type].back()*maphdr.precision;
					}

					if (io_group_size > 0)
					{
						if (p > 0)
						{
							if (p >= pixbatch_delim[0].back())
							{
								offset2 = 268 + pixbatch_delim[0].back() * pixbatch_size[0].back() * maphdr.precision;
								if (p >= pixbatch_delim[1].back())
								{
									offset2 += (pixbatch_delim[1].back()-pixbatch_delim[0].back()) * pixbatch_size[1].back() * maphdr.precision;
									offset2 += (p - pixbatch_delim[1].back()) * pixbatch_size[2].back() * maphdr.precision;
								}
								else offset2 += (p - pixbatch_delim[0].back()) * pixbatch_size[1].back() * maphdr.precision;
							}
							else
								offset2 = 268 + p * pixbatch_size[0].back() * maphdr.precision;
						}
						else if (parallel.rank() == (shell_write * parallel.size() + shell_outer - shell_inner) / (shell_outer + 1 - shell_inner)) offset2 = 0;
						else offset2 = 268;
					}
				}

				p = ((((shell + 1 - shell_inner) * parallel.size() + shell_outer - shell_inner) / (shell_outer + 1 - shell_inner)) - (((shell - shell_inner) * parallel.size() + shell_outer - shell_inner) / (shell_outer + 1 - shell_inner)));

				while (pix < pixbatch_id.size())
				{
					for (pixbatch_type = 0; pixbatch_delim[pixbatch_type].back() <= pixbatch_id[pix]; pixbatch_type++);

					if (p > 0 && pixbatch_delim[2].back() >= p && pixbatch_id[pix] / (pixbatch_delim[2].back() / p) < p)
					{
						for (n = 1; pix+n < pixbatch_id.size() && pixbatch_id[pix+n] == pixbatch_id[pix+n-1]+1 && pixbatch_id[pix+n] < pixbatch_delim[pixbatch_type].back() && (pixbatch_id[pix+n] / (pixbatch_delim[2].back() / p) == pixbatch_id[pix] / (pixbatch_delim[2].back() / p) || pixbatch_id[pix] / (pixbatch_delim[2].back() / p) == p-1); n++);
						for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
						{
							if (pixbuf[j][4] != NULL && pixbatch_size[pixbatch_type].back() > 0)
								parallel.send<Real>(pixbuf[j][4]+pix2, n*pixbatch_size[pixbatch_type].back(), (pixbatch_id[pix] / (pixbatch_delim[2].back() / p)) + ((shell - shell_inner) * parallel.size() + (io_group_size ? shell_outer - shell_inner : 0)) / (shell_outer + 1 - shell_inner));
						}
					}
					else
					{
						for (n = 1; pix+n < pixbatch_id.size() && pixbatch_id[pix+n] == pixbatch_id[pix+n-1]+1 && pixbatch_id[pix+n] < pixbatch_delim[pixbatch_type].back(); n++);
						for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
						{
							if (pixbuf[j][4] != NULL && pixbatch_size[pixbatch_type].back() > 0)
								parallel.send<Real>(pixbuf[j][4]+pix2, n*pixbatch_size[pixbatch_type].back(), (p ? p - 1 : 0) + ((shell - shell_inner) * parallel.size() + (p ? shell_outer - shell_inner : 0)) / (shell_outer + 1 - shell_inner));
						}
					}

					pix += n;
					pix2 += n*pixbatch_size[pixbatch_type].back();
				}

				offset.push_back(bytes);
				bytes += maphdr.Npix * maphdr.precision + 272;

				pixbatch_id.clear();
				sender_proc.clear();
			} // shell-loop

			if (io_group_size == 0)
				offset2 = 0;

			for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
			{
				if (pixbuf[j][0] == NULL || shell_outer < shell_inner) continue;

				if (sim.num_lightcone > 1)
				{
					if (j == LIGHTCONE_PHI_OFFSET)
						sprintf(filename, "%s%s%d_%04d_phi.map", sim.output_path, sim.basename_lightcone, i, cycle);
					else if (j == LIGHTCONE_CHI_OFFSET)
						sprintf(filename, "%s%s%d_%04d_chi.map", sim.output_path, sim.basename_lightcone, i, cycle);
					else if (j >= LIGHTCONE_B_OFFSET && j < LIGHTCONE_B_OFFSET+3)
						sprintf(filename, "%s%s%d_%04d_B%d.map", sim.output_path, sim.basename_lightcone, i, cycle, j+1-LIGHTCONE_B_OFFSET);
					else if (j >= LIGHTCONE_HIJ_OFFSET && j < LIGHTCONE_HIJ_OFFSET+5)
						sprintf(filename, "%s%s%d_%04d_h%d%d.map", sim.output_path, sim.basename_lightcone, i, cycle, (j - LIGHTCONE_HIJ_OFFSET < 3 ? 1 : 2), j - LIGHTCONE_HIJ_OFFSET + (j - LIGHTCONE_HIJ_OFFSET < 3 ? 1 : -1));
				}
				else
				{
					if (j == LIGHTCONE_PHI_OFFSET)
						sprintf(filename, "%s%s_%04d_phi.map", sim.output_path, sim.basename_lightcone, cycle);
					else if (j == LIGHTCONE_CHI_OFFSET)
						sprintf(filename, "%s%s_%04d_chi.map", sim.output_path, sim.basename_lightcone, cycle);
					else if (j >= LIGHTCONE_B_OFFSET && j < LIGHTCONE_B_OFFSET+3)
						sprintf(filename, "%s%s_%04d_B%d.map", sim.output_path, sim.basename_lightcone, cycle, j+1-LIGHTCONE_B_OFFSET);
					else if (j >= LIGHTCONE_HIJ_OFFSET && j < LIGHTCONE_HIJ_OFFSET+5)
						sprintf(filename, "%s%s_%04d_h%d%d.map", sim.output_path, sim.basename_lightcone, cycle, (j - LIGHTCONE_HIJ_OFFSET < 3 ? 1 : 2), j - LIGHTCONE_HIJ_OFFSET + (j - LIGHTCONE_HIJ_OFFSET < 3 ? 1 : -1));
				}

				MPI_File_open(parallel.lat_world_comm(), filename, MPI_MODE_WRONLY | MPI_MODE_CREATE,  MPI_INFO_NULL, &mapfile);
				MPI_File_set_size(mapfile, (MPI_Offset) bytes);
				MPI_File_write_at_all(mapfile, (MPI_Offset) offset[shell_write] + offset2, (void *) outbuf[j], bytes2, MPI_BYTE, &status);
				MPI_File_close(&mapfile);
			}

			for (j = 0; j < 3; j++)
			{
				pixbatch_size[j].clear();
				pixbatch_delim[j].clear();
			}

			offset.clear();

			for (j = 0; j < 9*LIGHTCONE_MAX_FIELDS; j++)
			{
				if (pixbuf[j/9][j%9] != NULL)
				{
					free(pixbuf[j/9][j%9]);
					pixbuf[j/9][j%9] = NULL;
				}
			}

			for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
			{
				if (outbuf[j] != NULL)
				{
					free(outbuf[j]);
					outbuf[j] = NULL;
				}
			}
#endif // HAVE_HEALPIX
		}
		else if (parallel.isRoot() && outfile != NULL)
			fclose(outfile);

		if (sim.out_lightcone[i] & MASK_GADGET && sim.lightcone[i].distance[0] > d - tau + 0.5 * dtau_old && sim.lightcone[i].distance[1] <= d - tau + 0.5 * dtau_old && d - tau + 0.5 * dtau_old > 0.)
		{
			n = findIntersectingLightcones(sim.lightcone[i], d - tau + (0.5 + LIGHTCONE_IDCHECK_ZONE) * dtau_old, d - tau - 0.5 * dtau, domain, vertex);

			hdr.num_files = 1;
			hdr.Omega0 = cosmo.Omega_m;
			hdr.OmegaLambda = cosmo.Omega_Lambda;
			hdr.HubbleParam = cosmo.h;
			hdr.BoxSize = sim.boxsize / GADGET_LENGTH_CONVERSION;
			hdr.flag_sfr = 0;
			hdr.flag_cooling = 0;
			hdr.flag_feedback = 0;
			hdr.flag_age = 0;
			hdr.flag_metals = 0;
			for (p = 0; p < 256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8 - 2 * 4 - 6 * 4; p++)
				hdr.fill[p] = 0;
			for (p = 0; p < 6; p++)
			{
				hdr.npart[p] = 0;
				hdr.npartTotal[p] = 0;
				hdr.npartTotalHW[p] = 0;
				hdr.mass[p] = 0.;
			}

			hdr.time = a;
			hdr.redshift = (1./a) - 1.;

			if (sim.baryon_flag)
				hdr.mass[1] = (double) sim.tracer_factor[0] * C_RHO_CRIT * cosmo.Omega_cdm * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[0] / GADGET_MASS_CONVERSION;
			else
				hdr.mass[1] = (double) sim.tracer_factor[0] * C_RHO_CRIT * (cosmo.Omega_cdm + cosmo.Omega_b) * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[0] / GADGET_MASS_CONVERSION;

			if (sim.num_lightcone > 1)
				sprintf(filename, "%d_%04d", i, cycle);
			else
				sprintf(filename, "_%04d", cycle);

			if (sim.tracer_factor[0] > 0)
				pcls_cdm->saveGadget2(h5filename + filename + "_cdm", hdr, sim.lightcone[i], d - tau, dtau, dtau_old, a * Hconf(a, fourpiG,//TODO_EB
				#ifdef HAVE_CLASS_BG
					H_spline, acc
				#else
					cosmo
				#endif
				), vertex, n, IDbacklog[0], IDprelog[0], phi, sim.tracer_factor[0]);

			if (sim.baryon_flag && sim.tracer_factor[1] > 0)
			{
				hdr.mass[1] = (double) sim.tracer_factor[1] * C_RHO_CRIT * cosmo.Omega_b * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[1] / GADGET_MASS_CONVERSION;
        pcls_b->saveGadget2(h5filename + filename + "_b", hdr, sim.lightcone[i], d - tau, dtau, dtau_old, a * Hconf(a, fourpiG,//TODO_EB
				#ifdef HAVE_CLASS_BG
					H_spline, acc
				#else
					cosmo
				#endif
				), vertex, n, IDbacklog[1], IDprelog[1], phi, sim.tracer_factor[1]);
			}

			for (p = 0; p < cosmo.num_ncdm; p++)
			{
				if (sim.numpcl[1+sim.baryon_flag+p] == 0 || sim.tracer_factor[p+1+sim.baryon_flag] == 0) continue;
				sprintf(buffer, "_ncdm%d", p);
				hdr.mass[1] = (double) sim.tracer_factor[p+1+sim.baryon_flag] * C_RHO_CRIT * cosmo.Omega_ncdm[p] * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[p+1+sim.baryon_flag] / GADGET_MASS_CONVERSION;
        pcls_ncdm[p].saveGadget2(h5filename + filename + buffer, hdr, sim.lightcone[i], d - tau, dtau, dtau_old, a * Hconf(a, fourpiG,//TODO_EB
				#ifdef HAVE_CLASS_BG
					H_spline, acc
				#else
					cosmo
				#endif
				), vertex, n, IDbacklog[p+1+sim.baryon_flag], IDprelog[p+1+sim.baryon_flag], phi, sim.tracer_factor[p+1+sim.baryon_flag]);
			}
		}
	}

#ifdef HAVE_HEALPIX
	delete[] outbuf;
#endif

	for (p = 0; p <= cosmo.num_ncdm + sim.baryon_flag; p++)
	{
		IDbacklog[p] = IDprelog[p];
		IDprelog[p].clear();

		n = IDbacklog[p].size();
		// dim 0 send/rec
		if (parallel.grid_rank()[0] % 2 == 0)
		{
			parallel.send_dim0<int>(n, (parallel.grid_size()[0]+parallel.grid_rank()[0]-1) % parallel.grid_size()[0]);
			parallel.receive_dim0<int>(i, (parallel.grid_size()[0]+parallel.grid_rank()[0]-1) % parallel.grid_size()[0]);
			parallel.send_dim0<int>(n, (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
			parallel.receive_dim0<int>(j, (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
		}
		else
		{
			parallel.receive_dim0<int>(i, (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
			parallel.send_dim0<int>(n, (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
			parallel.receive_dim0<int>(j, (parallel.grid_size()[0]+parallel.grid_rank()[0]-1) % parallel.grid_size()[0]);
			parallel.send_dim0<int>(n, (parallel.grid_size()[0]+parallel.grid_rank()[0]-1) % parallel.grid_size()[0]);
		}

		if (n+i+j > 0)
		{
			IDcombuf = (long *) malloc((n+i+j) * sizeof(long));

			n = 0;
			for (std::set<long>::iterator it = IDbacklog[p].begin(); it != IDbacklog[p].end(); it++)
				IDcombuf[n++] = *it;

			if (parallel.grid_rank()[0] % 2 == 0)
			{
				if (n > 0)
					parallel.send_dim0<long>(IDcombuf, n, (parallel.grid_size()[0]+parallel.grid_rank()[0]-1) % parallel.grid_size()[0]);
				if (i > 0)
					parallel.receive_dim0<long>(IDcombuf+n, i, (parallel.grid_size()[0]+parallel.grid_rank()[0]-1) % parallel.grid_size()[0]);
				if (n > 0)
					parallel.send_dim0<long>(IDcombuf, n, (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
				if (j > 0)
					parallel.receive_dim0<long>(IDcombuf+n+i, j, (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
			}
			else
			{
				if (i > 0)
					parallel.receive_dim0<long>(IDcombuf+n, i, (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
				if (n > 0)
					parallel.send_dim0<long>(IDcombuf, n, (parallel.grid_rank()[0]+1) % parallel.grid_size()[0]);
				if (j > 0)
					parallel.receive_dim0<long>(IDcombuf+n+i, j, (parallel.grid_size()[0]+parallel.grid_rank()[0]-1) % parallel.grid_size()[0]);
				if (n > 0)
					parallel.send_dim0<long>(IDcombuf, n, (parallel.grid_size()[0]+parallel.grid_rank()[0]-1) % parallel.grid_size()[0]);
			}

			n += i + j;

			for (i = IDbacklog[p].size(); i < n; i++)
				IDbacklog[p].insert(IDcombuf[i]);
		}

		// dim 1 send/rec
		if (parallel.grid_rank()[1] % 2 == 0)
		{
			parallel.send_dim1<int>(n, (parallel.grid_size()[1]+parallel.grid_rank()[1]-1) % parallel.grid_size()[1]);
			parallel.receive_dim1<int>(i, (parallel.grid_size()[1]+parallel.grid_rank()[1]-1) % parallel.grid_size()[1]);
			parallel.send_dim1<int>(n, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
			parallel.receive_dim1<int>(j, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);

			if (n > 0)
				parallel.send_dim1<long>(IDcombuf, n, (parallel.grid_size()[1]+parallel.grid_rank()[1]-1) % parallel.grid_size()[1]);

			if (i > 0)
			{
				IDcombuf2 = (long *) malloc(i * sizeof(long));
				parallel.receive_dim1<long>(IDcombuf2, i, (parallel.grid_size()[1]+parallel.grid_rank()[1]-1) % parallel.grid_size()[1]);
				while (i > 0)
					IDbacklog[p].insert(IDcombuf2[--i]);
				free(IDcombuf2);
			}

			if (n > 0)
			{
				parallel.send_dim1<long>(IDcombuf, n, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
				free(IDcombuf);
			}

			if (j > 0)
			{
				IDcombuf2 = (long *) malloc(j * sizeof(long));
				parallel.receive_dim1<long>(IDcombuf2, j, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
				while (j > 0)
					IDbacklog[p].insert(IDcombuf2[--j]);
				free(IDcombuf2);
			}
		}
		else
		{
			parallel.receive_dim1<int>(i, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
			parallel.send_dim1<int>(n, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
			parallel.receive_dim1<int>(j, (parallel.grid_size()[1]+parallel.grid_rank()[1]-1) % parallel.grid_size()[1]);
			parallel.send_dim1<int>(n, (parallel.grid_size()[1]+parallel.grid_rank()[1]-1) % parallel.grid_size()[1]);

			if (i > 0)
			{
				IDcombuf2 = (long *) malloc(i * sizeof(long));
				parallel.receive_dim1<long>(IDcombuf2, i, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);
				while (i > 0)
					IDbacklog[p].insert(IDcombuf2[--i]);
				free(IDcombuf2);
			}

			if (n > 0)
				parallel.send_dim1<long>(IDcombuf, n, (parallel.grid_rank()[1]+1) % parallel.grid_size()[1]);

			if (j > 0)
			{
				IDcombuf2 = (long *) malloc(j * sizeof(long));
				parallel.receive_dim1<long>(IDcombuf2, j, (parallel.grid_size()[1]+parallel.grid_rank()[1]-1) % parallel.grid_size()[1]);
				while (j > 0)
					IDbacklog[p].insert(IDcombuf2[--j]);
				free(IDcombuf2);
			}

			if (n > 0)
			{
				parallel.send_dim1<long>(IDcombuf, n, (parallel.grid_size()[1]+parallel.grid_rank()[1]-1) % parallel.grid_size()[1]);
				free(IDcombuf);
			}
		}
	}
}



//////////////////////////
// writeSpectra
//////////////////////////
// Description:
//   output of spectra
//
// Arguments:
//   sim            simulation metadata structure
//   cosmo          cosmological parameter structure
//   fourpiG        4 pi G (in code units)
//   a              scale factor
//   pkcount        spectrum output index
//   pcls_cdm       pointer to particle handler for CDM
//   pcls_b         pointer to particle handler for baryons
//   pcls_ncdm      array of particle handlers for
//                  non-cold DM (may be set to NULL)
//   phi            pointer to allocated field
//   chi            pointer to allocated field
//   Bi             pointer to allocated field
//   source         pointer to allocated field
//   Sij            pointer to allocated field
//   scalarFT       pointer to allocated field
//   BiFT           pointer to allocated field
//   SijFT          pointer to allocated field
//   plan_phi       pointer to FFT planner
//   plan_chi       pointer to FFT planner
//   plan_Bi        pointer to FFT planner
//   plan_source    pointer to FFT planner
//   plan_Sij       pointer to FFT planner
//   Bi_check       pointer to allocated field (or NULL)
//   BiFT_check     pointer to allocated field (or NULL)
//   plan_Bi_check  pointer to FFT planner (or NULL)
//
// Returns:
//
//////////////////////////

void writeSpectra(metadata & sim, cosmology & cosmo, const double fourpiG, const double a, const int pkcount,
#ifdef HAVE_CLASS
background & class_background, perturbs & class_perturbs, icsettings & ic,
#endif
Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_cdm, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_b, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_ncdm, Field<Real> * phi, Field<Real> * pi_k ,Field<Real> * zeta, Field<Real> * chi, Field<Real> * Bi,  Field<Real> * T00_Kess, Field<Real> * T0i_Kess, Field<Real> * Tij_Kess, Field<Real> * source, Field<Real> * Sij, Field<Cplx> * scalarFT, Field<Cplx> * scalarFT_pi, Field<Cplx> * scalarFT_zeta, Field<Cplx> * BiFT, Field<Cplx> * T00_KessFT, Field<Cplx> * T0i_KessFT, Field<Cplx> * Tij_KessFT, Field<Cplx> * SijFT, PlanFFT<Cplx> * plan_phi, PlanFFT<Cplx> * plan_pi_k , PlanFFT<Cplx> * plan_zeta, PlanFFT<Cplx> * plan_chi, PlanFFT<Cplx> * plan_Bi, PlanFFT<Cplx> * plan_T00_Kess, PlanFFT<Cplx> * plan_T0i_Kess, PlanFFT<Cplx> * plan_Tij_Kess, PlanFFT<Cplx> * plan_source, PlanFFT<Cplx> * plan_Sij
#ifdef CHECK_B
, Field<Real> * Bi_check, Field<Cplx> * BiFT_check, PlanFFT<Cplx> * plan_Bi_check
#endif
#ifdef VELOCITY
, Field<Real> * vi, Field<Cplx> * viFT, PlanFFT<Cplx> * plan_vi
#endif
)
{
	char filename[2*PARAM_MAX_LENGTH+24];
	char buffer[64];
	int i, j;
	Site x(phi->lattice());
	rKSite kFT(scalarFT->lattice());
	long numpts3d = (long) sim.numpts * (long) sim.numpts * (long) sim.numpts;
	Cplx tempk;
	double Omega_ncdm;

	Real * kbin;
	Real * power;
	Real * kscatter;
	Real * pscatter;
	int * occupation;

	kbin = (Real *) malloc(sim.numbins * sizeof(Real));
	power = (Real *) malloc(sim.numbins * sizeof(Real));
	kscatter = (Real *) malloc(sim.numbins * sizeof(Real));
	pscatter = (Real *) malloc(sim.numbins * sizeof(Real));
	occupation = (int *) malloc(sim.numbins * sizeof(int));

	if (sim.out_pk & MASK_RBARE || sim.out_pk & MASK_DBARE || sim.out_pk & MASK_POT || ((sim.out_pk & MASK_T00 || sim.out_pk & MASK_DELTA) && sim.gr_flag == 0))
	{
		projection_init(source);
#ifdef HAVE_CLASS
		if ((sim.radiation_flag > 0 || sim.fluid_flag > 0) && sim.gr_flag == 0)
		{
			projection_T00_project(class_background, class_perturbs, *source, *scalarFT, plan_source, sim, ic, cosmo, fourpiG, a);
			if (sim.out_pk & MASK_DELTA)
			{
				Omega_ncdm = 0;
				for (i = 0; i < cosmo.num_ncdm; i++)
				{
					if (a < 1. / (sim.z_switch_deltancdm[i] + 1.) && cosmo.Omega_ncdm[i] > 0)
						Omega_ncdm += bg_ncdm(a, cosmo, i);
				}
				// plan_source->execute(FFT_FORWARD);
				// extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				// sprintf(filename, "%s%s%03d_deltaclass.dat", sim.output_path, sim.basename_pk, pkcount);
				// writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * ((a < 1. / (sim.z_switch_deltarad + 1.) ? sim.radiation_flag : 0) * cosmo.Omega_rad / a + sim.fluid_flag * cosmo.Omega_fld / pow(a, 3. * cosmo.w0_fld) + Omega_ncdm) * ((a < 1. / (sim.z_switch_deltarad + 1.) ? sim.radiation_flag : 0) * cosmo.Omega_rad / a + sim.fluid_flag * cosmo.Omega_fld / pow(a, 3. * cosmo.w0_fld) + Omega_ncdm), filename, "power spectrum of delta for linear fields (CLASS)", a, sim.z_pk[pkcount]);
			}
		}
#endif
		scalarProjectionCIC_project(pcls_cdm, source);
		if (sim.baryon_flag)
			scalarProjectionCIC_project(pcls_b, source);
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			if (sim.numpcl[1+sim.baryon_flag+i] == 0) continue;
			scalarProjectionCIC_project(pcls_ncdm+i, source);
		}
		scalarProjectionCIC_comm(source);
		plan_source->execute(FFT_FORWARD);

		if (sim.out_pk & MASK_RBARE || sim.out_pk & MASK_DBARE || ((sim.out_pk & MASK_T00 || sim.out_pk & MASK_DELTA) && sim.gr_flag == 0))
			extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);

		if (sim.out_pk & MASK_RBARE)
		{
			sprintf(filename, "%s%s%03d_rhoN.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of rho_N", a, sim.z_pk[pkcount]);
		}

		if (sim.out_pk & MASK_DBARE)
		{
			// For deltaN we have k which is in unit of boxsize must divided by boxsize to give h/Mpc unit.
			//
			sprintf(filename, "%s%s%03d_deltaN.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_m * cosmo.Omega_m, filename, "power spectrum of delta_N", a, sim.z_pk[pkcount]);
		}

		if (sim.out_pk & MASK_T00 && sim.gr_flag == 0)
		{
			sprintf(filename, "%s%s%03d_T00.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of T00", a, sim.z_pk[pkcount]);
		}

		if (sim.out_pk & MASK_DELTA && sim.gr_flag == 0)
		{
			sprintf(filename, "%s%s%03d_delta.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_m * cosmo.Omega_m, filename, "power spectrum of delta", a, sim.z_pk[pkcount]);
		}

		if (sim.out_pk & MASK_POT)
		{
			solveModifiedPoissonFT(*scalarFT, *scalarFT, fourpiG / a);
			extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
			sprintf(filename, "%s%s%03d_psiN.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of psi_N", a, sim.z_pk[pkcount]);
		}

		if ((cosmo.num_ncdm > 0 || sim.baryon_flag) && (sim.out_pk & MASK_DBARE || (sim.out_pk & MASK_DELTA && sim.gr_flag == 0)))
		{
			projection_init(source);
			scalarProjectionCIC_project(pcls_cdm, source);
			scalarProjectionCIC_comm(source);
			plan_source->execute(FFT_FORWARD);
			extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
			sprintf(filename, "%s%s%03d_cdm.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * (sim.baryon_flag ? (cosmo.Omega_cdm * cosmo.Omega_cdm) : ((cosmo.Omega_cdm + cosmo.Omega_b) * (cosmo.Omega_cdm + cosmo.Omega_b))), filename, "power spectrum of delta_N for cdm", a, sim.z_pk[pkcount]);
			if (sim.baryon_flag)
			{
				// store k-space information for cross-spectra using SijFT as temporary array
				if (sim.out_pk & MASK_XSPEC)
				{
					for (kFT.first(); kFT.test(); kFT.next())
						(*SijFT)(kFT, 0) = (*scalarFT)(kFT);
				}
				projection_init(source);
				scalarProjectionCIC_project(pcls_b, source);
				scalarProjectionCIC_comm(source);
				plan_source->execute(FFT_FORWARD);
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				sprintf(filename, "%s%s%03d_b.dat", sim.output_path, sim.basename_pk, pkcount);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_b * cosmo.Omega_b, filename, "power spectrum of delta_N for baryons", a, sim.z_pk[pkcount]);
				if (sim.out_pk & MASK_XSPEC)
				{
					extractCrossSpectrum(*scalarFT, *SijFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
					sprintf(filename, "%s%s%03d_cdmxb.dat", sim.output_path, sim.basename_pk, pkcount);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_cdm * cosmo.Omega_b, filename, "cross power spectrum of delta_N for cdm x baryons", a, sim.z_pk[pkcount]);
				}
			}
			Omega_ncdm = 0.;
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				projection_init(source);
				if (sim.numpcl[1+sim.baryon_flag+i] > 0)
					scalarProjectionCIC_project(pcls_ncdm+i, source);
				scalarProjectionCIC_comm(source);
				plan_source->execute(FFT_FORWARD);
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				sprintf(filename, "%s%s%03d_ncdm%d.dat", sim.output_path, sim.basename_pk, pkcount, i);
				sprintf(buffer, "power spectrum of delta_N for ncdm %d", i);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_ncdm[i] * cosmo.Omega_ncdm[i], filename, buffer, a, sim.z_pk[pkcount]);
				Omega_ncdm += cosmo.Omega_ncdm[i];
				// store k-space information for cross-spectra using SijFT as temporary array
				if (cosmo.num_ncdm > 1 && i < 6)
				{
					for (kFT.first(); kFT.test(); kFT.next())
						(*SijFT)(kFT, i) = (*scalarFT)(kFT);
				}
			}
			if (cosmo.num_ncdm > 1 && cosmo.num_ncdm <= 7)
			{
				for (kFT.first(); kFT.test(); kFT.next())
				{
					for (i = 0; i < cosmo.num_ncdm-1; i++)
						(*scalarFT)(kFT) += (*SijFT)(kFT, i);
				}
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				sprintf(filename, "%s%s%03d_ncdm.dat", sim.output_path, sim.basename_pk, pkcount);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * Omega_ncdm * Omega_ncdm, filename, "power spectrum of delta_N for total ncdm", a, sim.z_pk[pkcount]);
			}
			if (cosmo.num_ncdm > 1)
			{
				for (i = 0; i < cosmo.num_ncdm-1 && i < 5; i++)
				{
					for (j = i+1; j < cosmo.num_ncdm && j < 6; j++)
					{
						if (sim.out_pk & MASK_XSPEC || (i == 0 && j == 1) || (i == 2 && j == 3) || (i == 4 && j == 5))
						{
							extractCrossSpectrum(*SijFT, *SijFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR, i, j);
							sprintf(filename, "%s%s%03d_ncdm%dx%d.dat", sim.output_path, sim.basename_pk, pkcount, i, j);
							sprintf(buffer, "cross power spectrum of delta_N for ncdm %d x %d", i, j);
							writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_ncdm[i] * cosmo.Omega_ncdm[j], filename, buffer, a, sim.z_pk[pkcount]);
						}
					}
				}
			}
		}
	}

	if (sim.out_pk & MASK_PHI)
	{
		plan_phi->execute(FFT_FORWARD);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s%03d_phi.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of phi", a, sim.z_pk[pkcount]);
	}


	 //KESSENCE PART

   // Note that according to definition and since pi is dimensionfull, so in the output we write \pi * H_conf which is dimensionless and can be compared to class and hiclass
	  if (sim.out_pk & MASK_PI_K)
		{
			plan_pi_k->execute(FFT_FORWARD);
			extractPowerSpectrum(*scalarFT_pi , kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
			sprintf(filename, "%s%s%03d_pi_k.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of pi_k  (dimensionless)", a, sim.z_pk[pkcount]);
		}

	     if (sim.out_pk & MASK_zeta)
		{
			plan_zeta->execute(FFT_FORWARD);
			extractPowerSpectrum(*scalarFT_zeta, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
			sprintf(filename, "%s%s%03d_zeta.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of zeta (dimensionless)", a, sim.z_pk[pkcount]);
		}

	   if (sim.out_pk & MASK_Delta_KESS)
		{
			// P (\delta)= deltarho_kess^2/ Omega_kess *a^(-3(1+w)) ) Omega_kess *a^(-3(1+w)) ) since in the defnition we have a^3 T00
			// We already included a^(-3) in the denominator, so we only need take the rest into account.
			plan_T00_Kess->execute(FFT_FORWARD);
			extractPowerSpectrum(*T00_KessFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
			sprintf(filename, "%s%s%03d_delta_kess.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI* cosmo.Omega_kessence * cosmo.Omega_kessence * pow(a, -3.* cosmo.w_kessence) * pow(a, -3.* cosmo.w_kessence), filename, "power spectrum of delta_kessence", a, sim.z_pk[pkcount]);
		 }

     // Phi_prime is dimensionful so we divide  to Hconf to make dimensionless!
     // if (sim.out_pk & MASK_PHI_PRIME)
     //   {
     //     phi_prime_plan->execute(FFT_FORWARD);
     //     extractPowerSpectrum(*phi_prime_scalarFT , kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
     //     sprintf(filename, "%s%s%03d_phi_prime.dat", sim.output_path, sim.basename_pk, pkcount);
     //     writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * (Hconf(a,fourpiG,cosmo) * Hconf(a,fourpiG,cosmo)), filename, "power spectrum of phi_prime/H (dimensionless)", a);
     //   }
		 //TODO_EB:here I did not change Hconf (in this function it appears only commented)

	   //KESSENCE END

	if (sim.out_pk & MASK_CHI)
	{
		plan_chi->execute(FFT_FORWARD);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s%03d_chi.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of chi", a, sim.z_pk[pkcount]);
	}

	if (sim.out_pk & MASK_HIJ)
	{
		projection_init(Sij);
		projection_Tij_project(pcls_cdm, Sij, a, phi);
		if (sim.baryon_flag)
			projection_Tij_project(pcls_b, Sij, a, phi);
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			if (sim.numpcl[1+sim.baryon_flag+i] == 0) continue;
			projection_Tij_project(pcls_ncdm+i, Sij, a, phi);
		}
		projection_Tij_comm(Sij);

		prepareFTsource<Real>(*phi, *Sij, *Sij, 2. * fourpiG / (double) sim.numpts / (double) sim.numpts / a);
		plan_Sij->execute(FFT_FORWARD);
		projectFTtensor(*SijFT, *SijFT);

		extractPowerSpectrum(*SijFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s%03d_hij.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, 2. * M_PI * M_PI, filename, "power spectrum of hij", a, sim.z_pk[pkcount]);
	}

	if ((sim.out_pk & MASK_T00 || sim.out_pk & MASK_DELTA) && sim.gr_flag > 0)
	{
		projection_init(source);
#ifdef HAVE_CLASS
		if (sim.radiation_flag > 0 || sim.fluid_flag > 0)
		{
			projection_T00_project(class_background, class_perturbs, *source, *scalarFT, plan_source, sim, ic, cosmo, fourpiG, a);
			if (sim.out_pk & MASK_DELTA)
			{
				Omega_ncdm = 0;
				for (i = 0; i < cosmo.num_ncdm; i++)
				{
					if (a < 1. / (sim.z_switch_deltancdm[i] + 1.) && cosmo.Omega_ncdm[i] > 0)
						Omega_ncdm += bg_ncdm(a, cosmo, i);
				}
				// plan_source->execute(FFT_FORWARD);
				// extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				// sprintf(filename, "%s%s%03d_deltaclass.dat", sim.output_path, sim.basename_pk, pkcount);
				// writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * ((a < 1. / (sim.z_switch_deltarad + 1.) ? sim.radiation_flag : 0) * cosmo.Omega_rad / a + sim.fluid_flag * cosmo.Omega_fld / pow(a, 3. * cosmo.w0_fld) + Omega_ncdm) * ((a < 1. / (sim.z_switch_deltarad + 1.) ? sim.radiation_flag : 0) * cosmo.Omega_rad / a + sim.fluid_flag * cosmo.Omega_fld / pow(a, 3. * cosmo.w0_fld) + Omega_ncdm), filename, "power spectrum of delta for linear fields (CLASS)", a, sim.z_pk[pkcount]);
			}
		}
#endif
		projection_T00_project(pcls_cdm, source, a, phi);
		if (sim.baryon_flag)
			projection_T00_project(pcls_b, source, a, phi);
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			if (sim.numpcl[1+sim.baryon_flag+i] == 0) continue;
			projection_T00_project(pcls_ncdm+i, source, a, phi);
		}
		projection_T00_comm(source);

		plan_source->execute(FFT_FORWARD);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);

		if (sim.out_pk & MASK_T00)
		{
			sprintf(filename, "%s%s%03d_T00.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of T00", a, sim.z_pk[pkcount]);
		}

		if (sim.out_pk & MASK_DELTA)
		{
			sprintf(filename, "%s%s%03d_delta.dat", sim.output_path, sim.basename_pk, pkcount);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * (cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) * (cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)), filename, "power spectrum of delta", a, sim.z_pk[pkcount]);
		}

    //Kessence Cross Power delta_kess * delta_m
    // cout<<"MASK_DELTAKESS_DELTA: "<<MASK_DELTAKESS_DELTA<<"MASK_Delta_KESS: "<<MASK_Delta_KESS<<endl;
    if ( sim.gr_flag > 0 && sim.out_pk & MASK_Delta_KESS && sim.out_pk & MASK_DELTAKESS_DELTA  && sim.out_pk & MASK_Delta_KESS)
    {
       // P (\deltam \delta_kess)= deltarho_kess * \delta_m / Omega_kess *a^(-3(1+w)) ) Omega_m *a^-3 since in the defnition we have a^3 T00
       // We already included a^(-3) in the denominator, so we only need take the rest into account.
       // Which are just a^{-3w} and Omega_kess and Omega_m
      extractCrossSpectrum(*scalarFT, *T00_KessFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
      sprintf(filename, "%s%s%03d_deltakess_deltam.dat", sim.output_path, sim.basename_pk, pkcount);
      writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * (cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) * cosmo.Omega_kessence *  pow(a, -3.* cosmo.w_kessence)  , filename, "cross power spectrum of delta_m and  delta_kess", a, sim.z_pk[pkcount]);
    }

		if (cosmo.num_ncdm > 0 || sim.baryon_flag)
		{
			projection_init(source);
			projection_T00_project(pcls_cdm, source, a, phi);
			projection_T00_comm(source);
			plan_source->execute(FFT_FORWARD);
			extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
			if (sim.out_pk & MASK_T00)
			{
				sprintf(filename, "%s%s%03d_T00cdm.dat", sim.output_path, sim.basename_pk, pkcount);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of T00 for cdm", a, sim.z_pk[pkcount]);
			}


			if (sim.out_pk & MASK_DELTA)
			{
				sprintf(filename, "%s%s%03d_deltacdm.dat", sim.output_path, sim.basename_pk, pkcount);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * (sim.baryon_flag ? (cosmo.Omega_cdm * cosmo.Omega_cdm) : ((cosmo.Omega_cdm + cosmo.Omega_b) * (cosmo.Omega_cdm + cosmo.Omega_b))), filename, "power spectrum of delta for cdm", a, sim.z_pk[pkcount]);
			}
			if (sim.baryon_flag)
			{
				// store k-space information for cross-spectra using SijFT as temporary array
				if (sim.out_pk & MASK_XSPEC)
				{
					for (kFT.first(); kFT.test(); kFT.next())
						(*SijFT)(kFT, 0) = (*scalarFT)(kFT);
				}
				projection_init(source);
				projection_T00_project(pcls_b, source, a, phi);
				projection_T00_comm(source);
				plan_source->execute(FFT_FORWARD);
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				if (sim.out_pk & MASK_T00)
				{
					sprintf(filename, "%s%s%03d_T00b.dat", sim.output_path, sim.basename_pk, pkcount);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of T00 for baryons", a, sim.z_pk[pkcount]);
				}
				if (sim.out_pk & MASK_DELTA)
				{
					sprintf(filename, "%s%s%03d_deltab.dat", sim.output_path, sim.basename_pk, pkcount);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_b * cosmo.Omega_b, filename, "power spectrum of delta for baryons", a, sim.z_pk[pkcount]);
				}
				if (sim.out_pk & MASK_XSPEC)
				{
					extractCrossSpectrum(*scalarFT, *SijFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
					sprintf(filename, "%s%s%03d_deltacdmxb.dat", sim.output_path, sim.basename_pk, pkcount);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * cosmo.Omega_b * cosmo.Omega_cdm, filename, "cross power spectrum of delta for cdm x baryons", a, sim.z_pk[pkcount]);
				}
			}
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				projection_init(source);
				if (sim.numpcl[1+sim.baryon_flag+i] > 0)
					projection_T00_project(pcls_ncdm+i, source, a, phi);
				projection_T00_comm(source);
				plan_source->execute(FFT_FORWARD);
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				if (sim.out_pk & MASK_T00)
				{
					sprintf(filename, "%s%s%03d_T00ncdm%d.dat", sim.output_path, sim.basename_pk, pkcount, i);
					sprintf(buffer, "power spectrum of T00 for ncdm %d", i);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, buffer, a, sim.z_pk[pkcount]);
				}
				if (sim.out_pk & MASK_DELTA)
				{
					sprintf(filename, "%s%s%03d_deltancdm%d.dat", sim.output_path, sim.basename_pk, pkcount, i);
					sprintf(buffer, "power spectrum of delta for ncdm %d", i);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * bg_ncdm(a, cosmo, i) * bg_ncdm(a, cosmo, i), filename, buffer, a, sim.z_pk[pkcount]);
				}
				// store k-space information for cross-spectra using SijFT as temporary array
				if (cosmo.num_ncdm > 1 && i < 6)
				{
					for (kFT.first(); kFT.test(); kFT.next())
						(*SijFT)(kFT, i) = (*scalarFT)(kFT);
				}
			}
			if (cosmo.num_ncdm > 1 && cosmo.num_ncdm <= 7)
			{
				for (kFT.first(); kFT.test(); kFT.next())
				{
					for (i = 0; i < cosmo.num_ncdm-1; i++)
						(*scalarFT)(kFT) += (*SijFT)(kFT, i);
				}
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				if (sim.out_pk & MASK_T00)
				{
					sprintf(filename, "%s%s%03d_T00ncdm.dat", sim.output_path, sim.basename_pk, pkcount);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, "power spectrum of T00 for total ncdm", a, sim.z_pk[pkcount]);
				}
				if (sim.out_pk & MASK_DELTA)
				{
					sprintf(filename, "%s%s%03d_deltancdm.dat", sim.output_path, sim.basename_pk, pkcount);
					writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * bg_ncdm(a, cosmo) * bg_ncdm(a, cosmo), filename, "power spectrum of delta for total ncdm", a, sim.z_pk[pkcount]);
				}
			}
			if (cosmo.num_ncdm > 1)
			{
				for (i = 0; i < cosmo.num_ncdm-1 && i < 5; i++)
				{
					for (j = i+1; j < cosmo.num_ncdm && j < 6; j++)
					{
						if (sim.out_pk & MASK_XSPEC || (i == 0 && j == 1) || (i == 2 && j == 3) || (i == 4 && j == 5))
						{
							extractCrossSpectrum(*SijFT, *SijFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR, i, j);
							if (sim.out_pk & MASK_T00)
							{
								sprintf(filename, "%s%s%03d_T00ncdm%dx%d.dat", sim.output_path, sim.basename_pk, pkcount, i, j);
								sprintf(buffer, "cross power spectrum of T00 for ncdm %d x %d", i, j);
								writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * pow(a, 6.0), filename, buffer, a, sim.z_pk[pkcount]);
							}
							if (sim.out_pk & MASK_DELTA)
							{
								sprintf(filename, "%s%s%03d_deltancdm%dx%d.dat", sim.output_path, sim.basename_pk, pkcount, i, j);
								sprintf(buffer, "cross power spectrum of delta for ncdm %d x %d", i, j);
								writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * bg_ncdm(a, cosmo, i) * bg_ncdm(a, cosmo, j), filename, buffer, a, sim.z_pk[pkcount]);
							}
						}
					}
				}
			}
		}
	}

	if (sim.out_pk & MASK_B)
	{
		extractPowerSpectrum(*BiFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s%03d_B.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, a * a * a * a * sim.numpts * sim.numpts * 2. * M_PI * M_PI, filename, "power spectrum of B", a, sim.z_pk[pkcount]);

#ifdef CHECK_B
		if (sim.vector_flag == VECTOR_PARABOLIC)
		{
			projection_init(Bi_check);
			projection_T0i_project(pcls_cdm, Bi_check, phi);
			if (sim.baryon_flag)
				projection_T0i_project(pcls_b, Bi_check, phi);
			for (i = 0; i < cosmo.num_ncdm; i++)
			{
				if (sim.numpcl[1+sim.baryon_flag+i] == 0) continue;
				projection_T0i_project(pcls_ncdm+i, Bi_check, phi);
			}
			projection_T0i_comm(Bi_check);
			plan_Bi_check->execute(FFT_FORWARD);
			projectFTvector(*BiFT_check, *BiFT_check, fourpiG / (double) sim.numpts / (double) sim.numpts);
		}
		extractPowerSpectrum(*BiFT_check, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s%03d_B_check.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, a * a * a * a * sim.numpts * sim.numpts * 2. * M_PI * M_PI, filename, "power spectrum of B", a, sim.z_pk[pkcount]);
#endif
	}

#ifdef VELOCITY
	if (sim.out_pk & MASK_VEL)
	{
		plan_vi->execute(FFT_FORWARD);
		extractPowerSpectrum(*viFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s%03d_v.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of velocity", a, sim.z_pk[pkcount]);

		projectFTtheta(*scalarFT, *viFT);
		extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s%03d_theta.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * sim.boxsize * sim.boxsize / cosmo.h / cosmo.h, filename, "power spectrum of theta (div v)", a, sim.z_pk[pkcount]);

		projectFTomega(*viFT);
		extractPowerSpectrum(*viFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
		sprintf(filename, "%s%s%03d_omega.dat", sim.output_path, sim.basename_pk, pkcount);
		writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * sim.boxsize * sim.boxsize / cosmo.h / cosmo.h, filename, "power spectrum of omega (curl v)", a, sim.z_pk[pkcount]);
	}
#endif

	free(kbin);
	free(power);
	free(kscatter);
	free(pscatter);
	free(occupation);
}

#endif


void writeSpectra_phi_prime(metadata & sim, cosmology & cosmo, const double fourpiG, const double a, const int pkcount,
	#ifdef HAVE_CLASS_BG
	gsl_spline * H_spline, gsl_interp_accel * acc,
	#endif
	Field<Real> * phi_prime ,Field<Cplx> * phi_prime_scalarFT , PlanFFT<Cplx> * phi_prime_plan)
{
char filename[2*PARAM_MAX_LENGTH+24];
char buffer[64];
int i, j;
Site x(phi_prime->lattice());
rKSite kFT(phi_prime_scalarFT->lattice());
  long numpts3d = (long) sim.numpts * (long) sim.numpts * (long) sim.numpts;
Cplx tempk;

Real * kbin;
Real * power;
Real * kscatter;
Real * pscatter;
int * occupation;

kbin = (Real *) malloc(sim.numbins * sizeof(Real));
power = (Real *) malloc(sim.numbins * sizeof(Real));
kscatter = (Real *) malloc(sim.numbins * sizeof(Real));
pscatter = (Real *) malloc(sim.numbins * sizeof(Real));
occupation = (int *) malloc(sim.numbins * sizeof(int));

double Hc = Hconf(a, fourpiG,//TODO_EB
	#ifdef HAVE_CLASS_BG
		H_spline, acc
	#else
		cosmo
	#endif
	);

// Phi_prime is dimensionful so we divide  to Hconf to make dimensionless!
if (sim.out_pk & MASK_PHI_PRIME)
  {
    phi_prime_plan->execute(FFT_FORWARD);
    extractPowerSpectrum(*phi_prime_scalarFT , kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
    sprintf(filename, "%s%s%03d_phi_prime.dat", sim.output_path, sim.basename_pk, pkcount);
    writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * (Hc * Hc), filename, "power spectrum of phi_prime/H (dimensionless)", a, sim.z_pk[pkcount]);
  }

free(kbin);
free(power);
free(kscatter);
free(pscatter);
free(occupation);
}


// Modified Poisson equation terms- TESTS
void writeSpectra_PoissonTerms(metadata & sim, cosmology & cosmo, const double fourpiG, const double a, const int pkcount,
	#ifdef HAVE_CLASS_BG
	gsl_spline * H_spline, gsl_interp_accel * acc,
	#endif
	Field<Real> * short_wave ,Field<Cplx> * short_wave_scalarFT , PlanFFT<Cplx> * short_wave_plan)
{
char filename[2*PARAM_MAX_LENGTH+24];
char buffer[64];
int i, j;
Site x(short_wave->lattice());
rKSite kFT(short_wave_scalarFT->lattice());
long numpts3d = (long) sim.numpts * (long) sim.numpts * (long) sim.numpts;
Cplx tempk;

Real * kbin;
Real * power;
Real * kscatter;
Real * pscatter;
int * occupation;

kbin = (Real *) malloc(sim.numbins * sizeof(Real));
power = (Real *) malloc(sim.numbins * sizeof(Real));
kscatter = (Real *) malloc(sim.numbins * sizeof(Real));
pscatter = (Real *) malloc(sim.numbins * sizeof(Real));
occupation = (int *) malloc(sim.numbins * sizeof(int));

double Hc = Hconf(a, fourpiG,//TODO_EB
	#ifdef HAVE_CLASS_BG
		H_spline, acc
	#else
		cosmo
	#endif
	);

    // It's dimensionless
    // relativistic_term_plan->execute(FFT_FORWARD);
    // extractPowerSpectrum(*relativistic_term_scalarFT , kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
    // sprintf(filename, "%s%s%03d_RE_phi_prime_over_dtau.dat", sim.output_path, sim.basename_pk, pkcount);
    // writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of relativistic term in Poisson equation; phi_prime/dtau (dimensionless)", a, sim.z_pk[pkcount]);


    // We divide by H_conf^2 to make it dimentionless!
    short_wave_plan->execute(FFT_FORWARD);
    extractPowerSpectrum(*short_wave_scalarFT , kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
    sprintf(filename, "%s%s%03d_short_wave.dat", sim.output_path, sim.basename_pk, pkcount);
    writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * (Hc * Hc,  Hc,  Hc) , filename, "power spectrum of short wave correction (dimensionless)", a, sim.z_pk[pkcount]);


    // stress_tensor_plan->execute(FFT_FORWARD);
    // extractPowerSpectrum(*stress_tensor_scalarFT , kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
    // sprintf(filename, "%s%s%03d_stress_tensor.dat", sim.output_path, sim.basename_pk, pkcount);
    // writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI , filename, "power spectrum of stress tensor part to the Poisson equation (dimensionless)", a, sim.z_pk[pkcount]);

free(kbin);
free(power);
free(kscatter);
free(pscatter);
free(occupation);
}
