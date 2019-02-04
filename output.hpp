//////////////////////////
// output.hpp
//////////////////////////
//
// Output of snapshots and spectra
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: November 2018
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
//   pcls_cdm       pointer to (uninitialized) particle handler for CDM
//   pcls_b         pointer to (uninitialized) particle handler for baryons
//   pcls_ncdm      array of (uninitialized) particle handlers for
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
void writeSnapshots(metadata & sim, cosmology & cosmo, const double fourpiG, const double a, const int done_hij, const int snapcount, string h5filename,
  #ifdef HAVE_CLASS
background & class_background, perturbs & class_perturbs, spectra & class_spectra, icsettings & ic,
#endif
Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_cdm, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_b, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_ncdm, Field<Real> * phi, Field<Real> * chi, Field<Real> * Bi, Field<Real> * source, Field<Real> * Sij, Field<Cplx> * scalarFT, Field<Cplx> * BiFT, Field<Cplx> * SijFT, PlanFFT<Cplx> * plan_phi, PlanFFT<Cplx> * plan_chi, PlanFFT<Cplx> * plan_Bi, PlanFFT<Cplx> * plan_source, PlanFFT<Cplx> * plan_Sij, Field<Real> * Bi_check = NULL, Field<Cplx> * BiFT_check = NULL, PlanFFT<Cplx> * plan_Bi_check = NULL)
{
	char filename[2*PARAM_MAX_LENGTH+24];
	char buffer[64];
	int i;
	gadget2_header hdr;
	Site x(phi->lattice());
	Real divB, curlB, divh, traceh, normh;

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

	if (sim.out_snapshot & MASK_B)
		Bi->saveHDF5_server_open(h5filename + filename + "_B");

	if (sim.out_snapshot & MASK_PHI)
		phi->saveHDF5_server_open(h5filename + filename + "_phi");

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

    #ifdef HAVE_CLASS
      if (sim.radiation_flag > 0 || sim.fluid_flag > 0)
      {
        projection_T00_project(class_background, class_perturbs, class_spectra, *source, *scalarFT, plan_source, sim, ic, cosmo, fourpiG, a);
           projection_T00_comm(source);
        #ifdef EXTERNAL_IO
           source->saveHDF5_server_write(NUMBER_OF_IO_FILES);
        #else
           if (sim.downgrade_factor > 1)
             source->saveHDF5_coarseGrain3D(h5filename + filename + "_T00_delta_kessence.h5", sim.downgrade_factor);
           else
             source->saveHDF5(h5filename + filename + "_T00_delta_kessence.h5");
        #endif
        }
    #endif
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

	if (sim.out_snapshot & MASK_GADGET)
	{
		if (sim.out_snapshot & MASK_MULTI)
			hdr.num_files = parallel.grid_size()[1];
		else
			hdr.num_files = 1;
		hdr.Omega0 = cosmo.Omega_m;
		hdr.OmegaLambda = cosmo.Omega_Lambda;
    hdr.OmegaLambda = 1.0 - cosmo.Omega_m;
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

		hdr.time = a;
		hdr.redshift = (1./a) - 1.;

		hdr.npart[1] = (uint32_t) ((sim.numpcl[0] / sim.tracer_factor[0]) % (1l << 32));
		hdr.npartTotal[1] = hdr.npart[1];
		hdr.npartTotalHW[1] = (uint32_t) ((sim.numpcl[0] / sim.tracer_factor[0]) >> 32);
		if (sim.baryon_flag)
			hdr.mass[1] = (double) sim.tracer_factor[0] * C_RHO_CRIT * cosmo.Omega_cdm * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[0] / GADGET_MASS_CONVERSION;
		else
			hdr.mass[1] = (double) sim.tracer_factor[0] * C_RHO_CRIT * (cosmo.Omega_cdm + cosmo.Omega_b) * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[0] / GADGET_MASS_CONVERSION;

		if (hdr.npartTotalHW[1] != 0 && hdr.num_files == 1)
		{
			COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": number of particles (" << (sim.numpcl[0] / sim.tracer_factor[0]) << ") in Gadget2 file exceeds limit (4294967295). Try using multi-Gadget2 output format." << endl;
		}
		else
			pcls_cdm->saveGadget2(h5filename + filename + "_cdm", hdr, sim.tracer_factor[0]);

		if (sim.baryon_flag)
		{
			hdr.npart[1] = (uint32_t) ((sim.numpcl[1] / sim.tracer_factor[1]) % (1l << 32));
			hdr.npartTotal[1] = hdr.npart[1];
			hdr.npartTotalHW[1] = (uint32_t) ((sim.numpcl[1] / sim.tracer_factor[1]) >> 32);
			hdr.mass[1] = (double) sim.tracer_factor[1] * C_RHO_CRIT * cosmo.Omega_b * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[1] / GADGET_MASS_CONVERSION;
			if (hdr.npartTotalHW[1] != 0 && hdr.num_files == 1)
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": number of particles (" << (sim.numpcl[1] / sim.tracer_factor[1]) << ") in Gadget2 file exceeds limit (4294967295). Try using multi-Gadget2 output format." << endl;
			}
			else
				pcls_b->saveGadget2(h5filename + filename + "_b", hdr, sim.tracer_factor[1]);
		}
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			if (sim.numpcl[1+sim.baryon_flag+i] == 0) continue;
			sprintf(buffer, "_ncdm%d", i);
			hdr.npart[1] = (uint32_t) ((sim.numpcl[i+1+sim.baryon_flag] / sim.tracer_factor[i+1+sim.baryon_flag]) % (1l << 32));
			hdr.npartTotal[1] = hdr.npart[1];
			hdr.npartTotalHW[1] = (uint32_t) ((sim.numpcl[i+1+sim.baryon_flag] / sim.tracer_factor[i+1+sim.baryon_flag]) >> 32);
			hdr.mass[1] = (double) sim.tracer_factor[i+1+sim.baryon_flag] * C_RHO_CRIT * cosmo.Omega_ncdm[i] * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[i+1+sim.baryon_flag] / GADGET_MASS_CONVERSION;
			if (hdr.npartTotalHW[1] != 0 && hdr.num_files == 1)
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": number of particles (" << (sim.numpcl[i+1+sim.baryon_flag] / sim.tracer_factor[i+1+sim.baryon_flag]) << ") in Gadget2 file exceeds limit (4294967295). Try using multi-Gadget2 output format." << endl;
			}
			else
				pcls_ncdm[i].saveGadget2(h5filename + filename + buffer, hdr, sim.tracer_factor[i+1+sim.baryon_flag]);
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


#if LIGHTCONE_THICKNESS > 3
#error This implementation of writeLightcones does not support LIGHTCONE_THICKNESS > 3
#endif

void writeLightcones(metadata & sim, cosmology & cosmo, const double fourpiG, const double a, const double tau, const double dtau, const double dtau_old, const double dtau_older, const double maxvel, const int cycle, string h5filename, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_cdm, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_b, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_ncdm, Field<Real> * phi, Field<Real> * chi, Field<Real> * Bi, Field<Real> * Sij, Field<Cplx> * BiFT, Field<Cplx> * SijFT, PlanFFT<Cplx> * plan_Bi, PlanFFT<Cplx> * plan_Sij, Field<Real> * lcbuffer[LIGHTCONE_MAX_FIELDS*LIGHTCONE_THICKNESS], Lattice * lclat[LIGHTCONE_MAX_FIELDS], int & done_hij)
{
	int i, j, n, p;
	double d;
	double vertex[MAX_INTERSECTS][3];
	double domain[6];
	double pos[3];
	double s[4];
	char filename[2*PARAM_MAX_LENGTH+24];
	char buffer[268];
	FILE * outfile;
	gadget2_header hdr;
	Site xsim;
	Site xlc;
	int done_B = 0;
#ifdef HAVE_HEALPIX
	Real * pixel[LIGHTCONE_MAX_FIELDS];
	char * outbuf;
	healpix_header map;
	double R[3][3];
	double w[3];
	double temp;
	int base_pos[3];
	int shell, shell_inner, shell_outer;
	int proc[LIGHTCONE_MAX_FIELDS];
	int proc_start[LIGHTCONE_MAX_FIELDS];
	int64_t bytes[LIGHTCONE_MAX_FIELDS];
	int64_t batch;
	uint32_t blocksize;
	MPI_File mapfile;
	MPI_Status status;
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
					fprintf(outfile, "# information file for lightcone %d\n# geometric parameters:\n# vertex = (%f, %f, %f) Mpc/h\n# redshift = %f\n# distance = (%f - %f) Mpc/h\n# opening half-angle = %f degrees\n# direction = (%f, %f, %f)\n# cycle   tau/boxsize    a              pcl_inner        pcl_outer        metric_inner", i, sim.lightcone[i].vertex[0]*sim.boxsize, sim.lightcone[i].vertex[1]*sim.boxsize, sim.lightcone[i].vertex[2]*sim.boxsize, sim.lightcone[i].z, sim.lightcone[i].distance[0]*sim.boxsize, sim.lightcone[i].distance[1]*sim.boxsize, (sim.lightcone[i].opening > -1.) ? acos(sim.lightcone[i].opening) * 180. / M_PI : 180., sim.lightcone[i].direction[0], sim.lightcone[i].direction[1], sim.lightcone[i].direction[2]);
				else
					fprintf(outfile, "# information file for lightcone\n# geometric parameters:\n# vertex = (%f, %f, %f) Mpc/h\n# redshift = %f\n# distance = (%f - %f) Mpc/h\n# opening half-angle = %f degrees\n# direction = (%f, %f, %f)\n# cycle   tau/boxsize    a              pcl_inner        pcl_outer        metric_inner", sim.lightcone[i].vertex[0]*sim.boxsize, sim.lightcone[i].vertex[1]*sim.boxsize, sim.lightcone[i].vertex[2]*sim.boxsize, sim.lightcone[i].z, sim.lightcone[i].distance[0]*sim.boxsize, sim.lightcone[i].distance[1]*sim.boxsize, (sim.lightcone[i].opening > -1.) ? acos(sim.lightcone[i].opening) * 180. / M_PI : 180., sim.lightcone[i].direction[0], sim.lightcone[i].direction[1], sim.lightcone[i].direction[2]);

#ifndef HAVE_HEALPIX
				for (j = 1; j < LIGHTCONE_THICKNESS; j++)
					fprintf(outfile, "     metric_(%d/%d)", j-1, j);
#endif

				fprintf(outfile, "     metric_outer\n");
			}
		}

		d = particleHorizon(1. / (1. + sim.lightcone[i].z), fourpiG, cosmo);

#ifdef HAVE_HEALPIX
		s[0] = d - tau - 0.5 * sim.covering[i] * dtau;
		s[1] = d - tau + 0.5 * sim.covering[i] * dtau_old;

		shell_inner = (s[0] > 0) ? ceil(s[0] * sim.numpts * sim.shellfactor[i]) : 0;
		shell_outer = (sim.lightcone[i].distance[0] > s[1]) ? floor(s[1] * sim.numpts * sim.shellfactor[i]) : floor(sim.lightcone[i].distance[0] * sim.numpts * sim.shellfactor[i]);
		if (shell_outer < shell_inner && s[1] > 0) shell_outer = shell_inner;

		map.precision = sizeof(Real);
		map.Ngrid = sim.numpts;
		map.direction[0] = sim.lightcone[i].direction[0];
		map.direction[1] = sim.lightcone[i].direction[1];
		map.direction[2] = sim.lightcone[i].direction[2];
		map.boxsize = sim.boxsize;

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

		if (sim.lightcone[i].distance[0] > s[0] && sim.lightcone[i].distance[1] <= s[1] && s[1] > 0.)
#else
#if LIGHTCONE_THICKNESS == 2
		s[0] = d - tau - dtau;
		s[1] = d - tau;
		s[2] = d - tau + dtau_old;
#elif LIGHTCONE_THICKNESS == 3
		s[0] = d - tau - 2. * dtau;
		s[1] = d - tau - 0.5 * dtau;
		s[2] = d - tau + 0.5 * dtau_old;
		s[3] = d - tau + dtau_old + 0.5 * dtau_older;
#else
		s[0] = d - tau - 0.5 * dtau;
		s[1] = d - tau + 0.5 * dtau_old;
#endif

		if (sim.lightcone[i].distance[0] > s[1] && sim.lightcone[i].distance[1] <= s[LIGHTCONE_THICKNESS] && s[LIGHTCONE_THICKNESS] > 0.)
#endif
		{
			if (parallel.isRoot() && outfile != NULL)
			{
				fprintf(outfile, "%6d   %e   %e   %2.12f   %2.12f   %2.12f", cycle, tau, a, d - tau - 0.5 * dtau, d - tau + 0.5 * dtau_old, s[0]);
#ifndef HAVE_HEALPIX
				for (j = 1; j < LIGHTCONE_THICKNESS; j++)
					fprintf(outfile, "   %2.12f", s[j]);
				fprintf(outfile, "   %2.12f\n", s[LIGHTCONE_THICKNESS]);
#else
				fprintf(outfile, "   %2.12f\n", s[1]);
#endif
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
					fwrite((const void *) s, sizeof(double), LIGHTCONE_THICKNESS+1, outfile);

					fclose(outfile);
				}
			}

#ifdef HAVE_HEALPIX
			outbuf = (char *) malloc(sizeof(Real) * PIXBUFFER);

			for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
				proc[j] = -1;

			j = (sim.out_lightcone[i] & MASK_PHI) ? 1 : 0;
			if (sim.out_lightcone[i] & MASK_CHI) j += 1;
			if (sim.out_lightcone[i] & MASK_B) j += 3;
			if (sim.out_lightcone[i] & MASK_HIJ) j += 5;

			if (j > parallel.size())
			{
				COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": number of tasks must be larger than number of field components in light cone output! Some data will not be written." << endl;
				j = parallel.size();
			}

			p = 0;
			if (sim.out_lightcone[i] & MASK_PHI)
			{
				pixel[LIGHTCONE_PHI_OFFSET] = (Real *) malloc(sizeof(Real) * PIXBUFFER);
				proc[LIGHTCONE_PHI_OFFSET] = p;
				p += parallel.size() / j;
			}
			if (sim.out_lightcone[i] & MASK_CHI)
			{
				pixel[LIGHTCONE_CHI_OFFSET] = (Real *) malloc(sizeof(Real) * PIXBUFFER);
				proc[LIGHTCONE_CHI_OFFSET] = p;
				p += parallel.size() / j;
			}
			if (sim.out_lightcone[i] & MASK_B)
			{
				pixel[LIGHTCONE_B_OFFSET] = (Real *) malloc(sizeof(Real) * PIXBUFFER);
				proc[LIGHTCONE_B_OFFSET] = p;
				p += parallel.size() / j;
				pixel[LIGHTCONE_B_OFFSET+1] = (Real *) malloc(sizeof(Real) * PIXBUFFER);
				proc[LIGHTCONE_B_OFFSET+1] = p;
				p += parallel.size() / j;
				pixel[LIGHTCONE_B_OFFSET+2] = (Real *) malloc(sizeof(Real) * PIXBUFFER);
				proc[LIGHTCONE_B_OFFSET+2] = p;
				p += parallel.size() / j;
			}
			if (sim.out_lightcone[i] & MASK_HIJ)
			{
				pixel[LIGHTCONE_HIJ_OFFSET] = (Real *) malloc(sizeof(Real) * PIXBUFFER);
				proc[LIGHTCONE_HIJ_OFFSET] = p;
				p += parallel.size() / j;
				pixel[LIGHTCONE_HIJ_OFFSET+1] = (Real *) malloc(sizeof(Real) * PIXBUFFER);
				proc[LIGHTCONE_HIJ_OFFSET+1] = p;
				p += parallel.size() / j;
				pixel[LIGHTCONE_HIJ_OFFSET+2] = (Real *) malloc(sizeof(Real) * PIXBUFFER);
				proc[LIGHTCONE_HIJ_OFFSET+2] = p;
				p += parallel.size() / j;
				pixel[LIGHTCONE_HIJ_OFFSET+3] = (Real *) malloc(sizeof(Real) * PIXBUFFER);
				proc[LIGHTCONE_HIJ_OFFSET+3] = p;
				p += parallel.size() / j;
				pixel[LIGHTCONE_HIJ_OFFSET+4] = (Real *) malloc(sizeof(Real) * PIXBUFFER);
				proc[LIGHTCONE_HIJ_OFFSET+4] = p;
				p += parallel.size() / j;
			}

			for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
				proc_start[j] = proc[j];

			for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
				bytes[j] = 0;

			batch = 0;

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

			for (shell = shell_inner; shell <= shell_outer; shell++)
			{
				map.distance = (double) shell / (double) sim.numpts / sim.shellfactor[i];

				for (map.Nside = sim.Nside[i][0]; map.Nside < sim.Nside[i][1]; map.Nside *= 2)
				{
					if (sim.pixelfactor[i] * 12. * map.Nside * map.Nside > 4. * M_PI * map.distance * map.distance * sim.numpts * sim.numpts) break;
				}

				map.Npix = 4;
				for (p = 2; p < 4 * map.Nside && map.Npix < (1. - sim.lightcone[i].opening) * 6 * map.Nside * map.Nside; p++)
				{
					if (p <= map.Nside) map.Npix += 4 * p;
					else if (p <= 3 * map.Nside - 1) map.Npix += 4 * map.Nside;
					else map.Npix += 4 * (4 * map.Nside - p);
				}

				blocksize = 256;
				memcpy((void *) buffer, (void *) &blocksize, sizeof(uint32_t));
				memcpy((void *) (buffer+4), (void *) &map, sizeof(healpix_header));
				memcpy((void *) (buffer+260), (void *) &blocksize, sizeof(uint32_t));
				blocksize = map.precision * map.Npix;
				memcpy((void *) (buffer+264), (void *) &blocksize, sizeof(uint32_t));

				for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
				{
					if (parallel.rank() == proc[j])
					{
						if ((PIXBUFFER * sizeof(Real)) - batch < 268)
						{
							memcpy((void *) (outbuf+batch), (void *) buffer, (sizeof(Real) * PIXBUFFER) - batch);
							batch = PIXBUFFER * sizeof(Real);
						}
						else
						{
							memcpy((void *) (outbuf+batch), (void *) buffer, 268);
							batch += 268;
						}
					}

					if (proc_start[j] >= 0 && proc_start[j] < parallel.size())
						bytes[j] += 268;

					if (bytes[j] / (PIXBUFFER * sizeof(Real)) > proc[j] - proc_start[j])
					{
						proc[j]++;
						if (parallel.rank() == proc[j])
						{
							batch = bytes[j] % (PIXBUFFER * sizeof(Real));
							memcpy((void *) outbuf, (void *) (buffer+268-batch), batch);
						}
					}
				}

				for (p = 0; p < map.Npix; p++)
				{
					pix2vec_ring64(map.Nside, p, w);

					pos[0] = (map.distance * (R[0][0] * w[0] + R[0][1] * w[1] + R[0][2] * w[2]) + sim.lightcone[i].vertex[0]) * sim.numpts;
					pos[1] = (map.distance * (R[1][0] * w[0] + R[1][1] * w[1] + R[1][2] * w[2]) + sim.lightcone[i].vertex[1]) * sim.numpts;
					pos[2] = (map.distance * (R[2][0] * w[0] + R[2][1] * w[1] + R[2][2] * w[2]) + sim.lightcone[i].vertex[2]) * sim.numpts;

					if (pos[0] >= 0)
					{
						w[0] = modf(pos[0], &temp);
						base_pos[0] = (int) temp % sim.numpts;
					}
					else
					{
						w[0] = 1. + modf(pos[0], &temp);
						base_pos[0] = sim.numpts - (1 - (int) temp) % sim.numpts;
					}
					if (pos[1] >= 0)
					{
						w[1] = modf(pos[1], &temp);
						base_pos[1] = (int) temp % sim.numpts;
					}
					else
					{
						w[1] = 1. + modf(pos[1], &temp);
						base_pos[1] = sim.numpts - (1 - (int) temp) % sim.numpts;
					}
					if (pos[2] >= 0)
					{
						w[2] = modf(pos[2], &temp);
						base_pos[2] = (int) temp % sim.numpts;
					}
					else
					{
						w[2] = 1. + modf(pos[2], &temp);
						base_pos[2] = sim.numpts - (1 - (int) temp) % sim.numpts;
					}


					if (xsim.setCoord(base_pos))
					{
						if (sim.out_lightcone[i] & MASK_PHI)
						{
							pixel[LIGHTCONE_PHI_OFFSET][p] = (1.-w[0]) * (1.-w[1]) * ((1.-w[2]) * (*phi)(xsim) + w[2] * (*phi)(xsim+2));
							pixel[LIGHTCONE_PHI_OFFSET][p] += w[0] * (1.-w[1]) * ((1.-w[2]) * (*phi)(xsim+0) + w[2] * (*phi)(xsim+0+2));
							pixel[LIGHTCONE_PHI_OFFSET][p] += w[0] * w[1] * ((1.-w[2]) * (*phi)(xsim+0+1) + w[2] * (*phi)(xsim+0+1+2));
							pixel[LIGHTCONE_PHI_OFFSET][p] += (1.-w[0]) * w[1] * ((1.-w[2]) * (*phi)(xsim+1) + w[2] * (*phi)(xsim+1+2));
						}
						if (sim.out_lightcone[i] & MASK_CHI)
						{
							pixel[LIGHTCONE_CHI_OFFSET][p] = (1.-w[0]) * (1.-w[1]) * ((1.-w[2]) * (*chi)(xsim) + w[2] * (*chi)(xsim+2));
							pixel[LIGHTCONE_CHI_OFFSET][p] += w[0] * (1.-w[1]) * ((1.-w[2]) * (*chi)(xsim+0) + w[2] * (*chi)(xsim+0+2));
							pixel[LIGHTCONE_CHI_OFFSET][p] += w[0] * w[1] * ((1.-w[2]) * (*chi)(xsim+0+1) + w[2] * (*chi)(xsim+0+1+2));
							pixel[LIGHTCONE_CHI_OFFSET][p] += (1.-w[0]) * w[1] * ((1.-w[2]) * (*chi)(xsim+1) + w[2] * (*chi)(xsim+1+2));
						}
						if (sim.out_lightcone[i] & MASK_B)
						{
#ifdef LIGHTCONE_INTERPOLATE
							pixel[LIGHTCONE_B_OFFSET][p] = (1.-w[2]) * (1.-w[1]) * ((1.-w[0]) * (*Bi)(xsim-0,0) + w[0] * (*Bi)(xsim+0,0) + (*Bi)(xsim,0));
							pixel[LIGHTCONE_B_OFFSET][p] += w[2] * (1.-w[1]) * ((1.-w[0]) * (*Bi)(xsim-0+2,0) + w[0] * (*Bi)(xsim+0+2,0) + (*Bi)(xsim+2,0));
							pixel[LIGHTCONE_B_OFFSET][p] += w[2] * w[1] * ((1.-w[0]) * (*Bi)(xsim-0+1+2,0) + w[0] * (*Bi)(xsim+0+1+2,0) + (*Bi)(xsim+1+2,0));
							pixel[LIGHTCONE_B_OFFSET][p] += (1.-w[2]) * w[1] * ((1.-w[0]) * (*Bi)(xsim-0+1,0) + w[0] * (*Bi)(xsim+0+1,0) + (*Bi)(xsim+1,0));

							pixel[LIGHTCONE_B_OFFSET+1][p] = (1.-w[2]) * (1.-w[0]) * ((1.-w[1]) * (*Bi)(xsim-1,1) + w[1] * (*Bi)(xsim+1,1) + (*Bi)(xsim,1));
							pixel[LIGHTCONE_B_OFFSET+1][p] += w[2] * (1.-w[0]) * ((1.-w[1]) * (*Bi)(xsim-1+2,1) + w[1] * (*Bi)(xsim+1+2,1) + (*Bi)(xsim+2,1));
							pixel[LIGHTCONE_B_OFFSET+1][p] += w[2] * w[0] * ((1.-w[1]) * (*Bi)(xsim+0-1+2,1) + w[1] * (*Bi)(xsim+0+1+2,1) + (*Bi)(xsim+0+2,1));
							pixel[LIGHTCONE_B_OFFSET+1][p] += (1.-w[2]) * w[0] * ((1.-w[1]) * (*Bi)(xsim+0-1,1) + w[1] * (*Bi)(xsim+0+1,1) + (*Bi)(xsim+0,1));

							pixel[LIGHTCONE_B_OFFSET+2][p] = (1.-w[0]) * (1.-w[1]) * ((1.-w[2]) * (*Bi)(xsim-2,2) + w[2] * (*Bi)(xsim+2,2) + (*Bi)(xsim,2));
							pixel[LIGHTCONE_B_OFFSET+2][p] += w[0] * (1.-w[1]) * ((1.-w[2]) * (*Bi)(xsim+0-2,2) + w[2] * (*Bi)(xsim+0+2,2) + (*Bi)(xsim+0,2));
							pixel[LIGHTCONE_B_OFFSET+2][p] += w[0] * w[1] * ((1.-w[2]) * (*Bi)(xsim+0+1-2,2) + w[2] * (*Bi)(xsim+0+1+2,2) + (*Bi)(xsim+0+1,2));
							pixel[LIGHTCONE_B_OFFSET+2][p] += (1.-w[0]) * w[1] * ((1.-w[2]) * (*Bi)(xsim+1-2,2) + w[2] * (*Bi)(xsim+1+2,2) + (*Bi)(xsim+1,2));

							pixel[LIGHTCONE_B_OFFSET][p] /= 2. * a * a * sim.numpts;
							pixel[LIGHTCONE_B_OFFSET+1][p] /= 2. * a * a * sim.numpts;
							pixel[LIGHTCONE_B_OFFSET+2][p] /= 2. * a * a * sim.numpts;
#else
							if (w[0] > 0.5)
							{
								pixel[LIGHTCONE_B_OFFSET][p] = (1.5-w[0]) * (1.-w[1]) * ((1.-w[2]) * (*Bi)(xsim,0) + w[2] * (*Bi)(xsim+2,0));
								pixel[LIGHTCONE_B_OFFSET][p] += (1.5-w[0]) * w[1] * ((1.-w[2]) * (*Bi)(xsim+1,0) + w[2] * (*Bi)(xsim+1+2,0));
								pixel[LIGHTCONE_B_OFFSET][p] += (w[0]-0.5) * w[1] * ((1.-w[2]) * (*Bi)(xsim+0+1,0) + w[2] * (*Bi)(xsim+0+1+2,0));
								pixel[LIGHTCONE_B_OFFSET][p] += (w[0]-0.5) * (1.-w[1]) * ((1.-w[2]) * (*Bi)(xsim+0,0) + w[2] * (*Bi)(xsim+0+2,0));
							}
							else
							{
								pixel[LIGHTCONE_B_OFFSET][p] = (0.5+w[0]) * (1.-w[1]) * ((1.-w[2]) * (*Bi)(xsim,0) + w[2] * (*Bi)(xsim+2,0));
								pixel[LIGHTCONE_B_OFFSET][p] += (0.5+w[0]) * w[1] * ((1.-w[2]) * (*Bi)(xsim+1,0) + w[2] * (*Bi)(xsim+1+2,0));
								pixel[LIGHTCONE_B_OFFSET][p] += (0.5-w[0]) * w[1] * ((1.-w[2]) * (*Bi)(xsim-0+1,0) + w[2] * (*Bi)(xsim-0+1+2,0));
								pixel[LIGHTCONE_B_OFFSET][p] += (0.5-w[0]) * (1.-w[1]) * ((1.-w[2]) * (*Bi)(xsim-0,0) + w[2] * (*Bi)(xsim-0+2,0));
							}
							if (w[1] > 0.5)
							{
								pixel[LIGHTCONE_B_OFFSET+1][p] = (1.-w[0]) * (1.5-w[1]) * ((1.-w[2]) * (*Bi)(xsim,1) + w[2] * (*Bi)(xsim+2,1));
								pixel[LIGHTCONE_B_OFFSET+1][p] += (1.-w[0]) * (w[1]-0.5) * ((1.-w[2]) * (*Bi)(xsim+1,1) + w[2] * (*Bi)(xsim+1+2,1));
								pixel[LIGHTCONE_B_OFFSET+1][p] += w[0] * (w[1]-0.5) * ((1.-w[2]) * (*Bi)(xsim+0+1,1) + w[2] * (*Bi)(xsim+0+1+2,1));
								pixel[LIGHTCONE_B_OFFSET+1][p] += w[0] * (1.5-w[1]) * ((1.-w[2]) * (*Bi)(xsim+0,1) + w[2] * (*Bi)(xsim+0+2,1));
							}
							else
							{
								pixel[LIGHTCONE_B_OFFSET+1][p] = (1.-w[0]) * (0.5+w[1]) * ((1.-w[2]) * (*Bi)(xsim,1) + w[2] * (*Bi)(xsim+2,1));
								pixel[LIGHTCONE_B_OFFSET+1][p] += (1.-w[0]) * (0.5-w[1]) * ((1.-w[2]) * (*Bi)(xsim-1,1) + w[2] * (*Bi)(xsim-1+2,1));
								pixel[LIGHTCONE_B_OFFSET+1][p] += w[0] * (0.5-w[1]) * ((1.-w[2]) * (*Bi)(xsim+0-1,1) + w[2] * (*Bi)(xsim+0-1+2,1));
								pixel[LIGHTCONE_B_OFFSET+1][p] += w[0] * (0.5+w[1]) * ((1.-w[2]) * (*Bi)(xsim+0,1) + w[2] * (*Bi)(xsim+0+2,1));
							}
							if (w[2] > 0.5)
							{
								pixel[LIGHTCONE_B_OFFSET+2][p] = (1.-w[0]) * (1.-w[1]) * ((1.5-w[2]) * (*Bi)(xsim,2) + (w[2]-0.5) * (*Bi)(xsim+2,2));
								pixel[LIGHTCONE_B_OFFSET+2][p] += (1.-w[0]) * w[1] * ((1.5-w[2]) * (*Bi)(xsim+1,2) + (w[2]-0.5) * (*Bi)(xsim+1+2,2));
								pixel[LIGHTCONE_B_OFFSET+2][p] += w[0] * w[1] * ((1.5-w[2]) * (*Bi)(xsim+0+1,2) + (w[2]-0.5) * (*Bi)(xsim+0+1+2,2));
								pixel[LIGHTCONE_B_OFFSET+2][p] += w[0] * (1.-w[1]) * ((1.5-w[2]) * (*Bi)(xsim+0,2) + (w[2]-0.5) * (*Bi)(xsim+0+2,2));
							}
							else
							{
								pixel[LIGHTCONE_B_OFFSET+2][p] = (1.-w[0]) * (1.-w[1]) * ((0.5+w[2]) * (*Bi)(xsim,2) + (0.5-w[2]) * (*Bi)(xsim-2,2));
								pixel[LIGHTCONE_B_OFFSET+2][p] += (1.-w[0]) * w[1] * ((0.5+w[2]) * (*Bi)(xsim+1,2) + (0.5-w[2]) * (*Bi)(xsim+1-2,2));
								pixel[LIGHTCONE_B_OFFSET+2][p] += w[0] * w[1] * ((0.5+w[2]) * (*Bi)(xsim+0+1,2) + (0.5-w[2]) * (*Bi)(xsim+0+1-2,2));
								pixel[LIGHTCONE_B_OFFSET+2][p] += w[0] * (1.-w[1]) * ((0.5+w[2]) * (*Bi)(xsim+0,2) + (0.5-w[2]) * (*Bi)(xsim+0-2,2));
							}
							pixel[LIGHTCONE_B_OFFSET][p] /= a * a * sim.numpts;
							pixel[LIGHTCONE_B_OFFSET+1][p] /= a * a * sim.numpts;
							pixel[LIGHTCONE_B_OFFSET+2][p] /= a * a * sim.numpts;
#endif
						}
						if (sim.out_lightcone[i] & MASK_HIJ)
						{
							pixel[LIGHTCONE_HIJ_OFFSET][p] = (1.-w[0]) * (1.-w[1]) * ((1.-w[2]) * (*Sij)(xsim,0,0) + w[2] * (*Sij)(xsim+2,0,0));
							pixel[LIGHTCONE_HIJ_OFFSET][p] += w[0] * (1.-w[1]) * ((1.-w[2]) * (*Sij)(xsim+0,0,0) + w[2] * (*Sij)(xsim+0+2,0,0));
							pixel[LIGHTCONE_HIJ_OFFSET][p] += w[0] * w[1] * ((1.-w[2]) * (*Sij)(xsim+0+1,0,0) + w[2] * (*Sij)(xsim+0+1+2,0,0));
							pixel[LIGHTCONE_HIJ_OFFSET][p] += (1.-w[0]) * w[1] * ((1.-w[2]) * (*Sij)(xsim+1,0,0) + w[2] * (*Sij)(xsim+1+2,0,0));

							pixel[LIGHTCONE_HIJ_OFFSET+1][p] = (1.-w[2]) * 0.25 * ((*Sij)(xsim,0,1) + (1.-w[0]) * ((*Sij)(xsim-0,0,1) + (1.-w[1]) * (*Sij)(xsim-0-1,0,1) + w[1] * (*Sij)(xsim-0+1,0,1)) + w[0] * ((*Sij)(xsim+0,0,1) + (1.-w[1]) * (*Sij)(xsim+0-1,0,1) + w[1] * (*Sij)(xsim+0+1,0,1)) + (1.-w[1]) * (*Sij)(xsim-1,0,1) + w[1] * (*Sij)(xsim+1,0,1));
							pixel[LIGHTCONE_HIJ_OFFSET+1][p] += w[2] * 0.25 * ((*Sij)(xsim+2,0,1) + (1.-w[0]) * ((*Sij)(xsim-0+2,0,1) + (1.-w[1]) * (*Sij)(xsim-0-1+2,0,1) + w[1] * (*Sij)(xsim-0+1+2,0,1)) + w[0] * ((*Sij)(xsim+0+2,0,1) + (1.-w[1]) * (*Sij)(xsim+0-1+2,0,1) + w[1] * (*Sij)(xsim+0+1+2,0,1)) + (1.-w[1]) * (*Sij)(xsim-1+2,0,1) + w[1] * (*Sij)(xsim+1+2,0,1));

							pixel[LIGHTCONE_HIJ_OFFSET+2][p] = (1.-w[1]) * 0.25 * ((*Sij)(xsim,0,2) + (1.-w[0]) * ((*Sij)(xsim-0,0,2) + (1.-w[2]) * (*Sij)(xsim-0-2,0,2) + w[2] * (*Sij)(xsim-0+2,0,2)) + w[0] * ((*Sij)(xsim+0,0,2) + (1.-w[2]) * (*Sij)(xsim+0-2,0,2) + w[2] * (*Sij)(xsim+0+2,0,2)) + (1.-w[2]) * (*Sij)(xsim-2,0,2) + w[2] * (*Sij)(xsim+2,0,2));
							pixel[LIGHTCONE_HIJ_OFFSET+2][p] += w[1] * 0.25 * ((*Sij)(xsim+1,0,2) + (1.-w[0]) * ((*Sij)(xsim-0+1,0,2) + (1.-w[2]) * (*Sij)(xsim-0+1-2,0,2) + w[2] * (*Sij)(xsim-0+1+2,0,2)) + w[0] * ((*Sij)(xsim+0+1,0,2) + (1.-w[2]) * (*Sij)(xsim+0+1-2,0,2) + w[2] * (*Sij)(xsim+0+1+2,0,2)) + (1.-w[2]) * (*Sij)(xsim+1-2,0,2) + w[2] * (*Sij)(xsim+1+2,0,2));


							pixel[LIGHTCONE_HIJ_OFFSET+3][p] = (1.-w[0]) * (1.-w[1]) * ((1.-w[2]) * (*Sij)(xsim,1,1) + w[2] * (*Sij)(xsim+2,1,1));
							pixel[LIGHTCONE_HIJ_OFFSET+3][p] += w[0] * (1.-w[1]) * ((1.-w[2]) * (*Sij)(xsim+0,1,1) + w[2] * (*Sij)(xsim+0+2,1,1));
							pixel[LIGHTCONE_HIJ_OFFSET+3][p] += w[0] * w[1] * ((1.-w[2]) * (*Sij)(xsim+0+1,1,1) + w[2] * (*Sij)(xsim+0+1+2,1,1));
							pixel[LIGHTCONE_HIJ_OFFSET+3][p] += (1.-w[0]) * w[1] * ((1.-w[2]) * (*Sij)(xsim+1,1,1) + w[2] * (*Sij)(xsim+1+2,1,1));

							pixel[LIGHTCONE_HIJ_OFFSET+4][p] = (1.-w[0]) * 0.25 * ((*Sij)(xsim,1,2) + (1.-w[1]) * ((*Sij)(xsim-1,1,2) + (1.-w[2]) * (*Sij)(xsim-1-2,1,2) + w[2] * (*Sij)(xsim-1+2,1,2)) + w[1] * ((*Sij)(xsim+1,1,2) + (1.-w[2]) * (*Sij)(xsim+1-2,1,2) + w[2] * (*Sij)(xsim+1+2,1,2)) + (1.-w[2]) * (*Sij)(xsim-2,1,2) + w[2] * (*Sij)(xsim+2,1,2));
							pixel[LIGHTCONE_HIJ_OFFSET+4][p] += w[0] * 0.25 * ((*Sij)(xsim+0,1,2) + (1.-w[1]) * ((*Sij)(xsim+0-1,1,2) + (1.-w[2]) * (*Sij)(xsim+0-1-2,1,2) + w[2] * (*Sij)(xsim+0-1+2,1,2)) + w[1] * ((*Sij)(xsim+0+1,1,2) + (1.-w[2]) * (*Sij)(xsim+0+1-2,1,2) + w[2] * (*Sij)(xsim+0+1+2,1,2)) + (1.-w[2]) * (*Sij)(xsim+0-2,1,2) + w[2] * (*Sij)(xsim+0+2,1,2));
						}
					}
					else
					{
						if (sim.out_lightcone[i] & MASK_PHI)
							pixel[LIGHTCONE_PHI_OFFSET][p] = 0;
						if (sim.out_lightcone[i] & MASK_CHI)
							pixel[LIGHTCONE_CHI_OFFSET][p] = 0;
						if (sim.out_lightcone[i] & MASK_B)
						{
							pixel[LIGHTCONE_B_OFFSET][p] = 0;
							pixel[LIGHTCONE_B_OFFSET+1][p] = 0;
							pixel[LIGHTCONE_B_OFFSET+2][p] = 0;
						}
						if (sim.out_lightcone[i] & MASK_HIJ)
						{
							pixel[LIGHTCONE_HIJ_OFFSET][p] = 0;
							pixel[LIGHTCONE_HIJ_OFFSET+1][p] = 0;
							pixel[LIGHTCONE_HIJ_OFFSET+2][p] = 0;
							pixel[LIGHTCONE_HIJ_OFFSET+3][p] = 0;
							pixel[LIGHTCONE_HIJ_OFFSET+4][p] = 0;
						}
					}

					if ((p+1) % PIXBUFFER == 0)
					{
						for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
						{
							if (parallel.rank() == proc[j])
							{
#ifdef SINGLE
								MPI_Reduce(MPI_IN_PLACE, pixel[j], PIXBUFFER, MPI_FLOAT, MPI_SUM, proc[j], parallel.lat_world_comm());
#else
								MPI_Reduce(MPI_IN_PLACE, pixel[j], PIXBUFFER, MPI_DOUBLE, MPI_SUM, proc[j], parallel.lat_world_comm());
#endif
								memcpy((void *) (outbuf+batch), (void *) pixel[j], (PIXBUFFER * sizeof(Real)) - batch);
								if (batch > 0)
									parallel.send<char>(((char *) pixel[j]) + (PIXBUFFER * sizeof(Real)) - batch, batch, proc[j]+1);
								batch = PIXBUFFER * sizeof(Real);
							}
							else if (proc[j] >= parallel.size() && proc_start[j] < parallel.size())
							{
								COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": insufficient PIXBUFFER for parallel Healpix output!" << endl;
								parallel.abortForce();
							}
							else if (proc[j] >= 0 && proc[j] < parallel.size())
							{
#ifdef SINGLE
								MPI_Reduce(pixel[j], NULL, PIXBUFFER, MPI_FLOAT, MPI_SUM, proc[j], parallel.lat_world_comm());
#else
								MPI_Reduce(pixel[j], NULL, PIXBUFFER, MPI_DOUBLE, MPI_SUM, proc[j], parallel.lat_world_comm());
#endif
							}
							if (proc[j] >= 0)
							{
								proc[j]++;
								bytes[j] += PIXBUFFER * sizeof(Real);
								if (parallel.rank() == proc[j] && bytes[j] % (PIXBUFFER * sizeof(Real)) > 0)
								{
									if (batch > 0)
									{
										cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": task " << parallel.rank() << " insufficient PIXBUFFER for parallel Healpix output!" << endl;
										exit(-200);
									}
									batch = bytes[j] % (PIXBUFFER * sizeof(Real));
									parallel.receive<char>(outbuf, batch, proc[j]-1);
								}
							}
						}
					}
				}

				if (map.Npix % PIXBUFFER > 0)
				{
					for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
					{
						if (parallel.rank() == proc[j])
						{
#ifdef SINGLE
							MPI_Reduce(MPI_IN_PLACE, pixel[j], map.Npix % PIXBUFFER, MPI_FLOAT, MPI_SUM, proc[j], parallel.lat_world_comm());
#else
							MPI_Reduce(MPI_IN_PLACE, pixel[j], map.Npix % PIXBUFFER, MPI_DOUBLE, MPI_SUM, proc[j], parallel.lat_world_comm());
#endif
							if ((map.Npix % PIXBUFFER) * sizeof(Real) <= (PIXBUFFER * sizeof(Real)) - batch)
							{
								memcpy((void *) (outbuf+batch), (void *) pixel[j], (map.Npix % PIXBUFFER) * sizeof(Real));
								batch += (map.Npix % PIXBUFFER) * sizeof(Real);
							}
							else
							{
								memcpy((void *) (outbuf+batch), (void *) pixel[j], (PIXBUFFER * sizeof(Real)) - batch);
								parallel.send<char>(((char *) pixel[j]) + (PIXBUFFER * sizeof(Real)) - batch, ((map.Npix % PIXBUFFER) - PIXBUFFER) * sizeof(Real) + batch, proc[j]+1);
								batch = PIXBUFFER * sizeof(Real);
							}
						}
						else if (proc[j] >= parallel.size() && proc_start[j] < parallel.size())
						{
							COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": insufficient PIXBUFFER for parallel Healpix output!" << endl;
							parallel.abortForce();
						}
						else if (proc[j] >= 0 && proc[j] < parallel.size())
						{
#ifdef SINGLE
							MPI_Reduce(pixel[j], NULL, map.Npix % PIXBUFFER, MPI_FLOAT, MPI_SUM, proc[j], parallel.lat_world_comm());
#else
							MPI_Reduce(pixel[j], NULL, map.Npix % PIXBUFFER, MPI_DOUBLE, MPI_SUM, proc[j], parallel.lat_world_comm());
#endif
						}

						if (proc[j] >= 0)
						{
							bytes[j] += (map.Npix % PIXBUFFER) * sizeof(Real);

							if (bytes[j] / (PIXBUFFER * sizeof(Real)) > proc[j] - proc_start[j])
							{
								proc[j]++;
								if (parallel.rank() == proc[j])
								{
									if (batch > 0)
									{
										cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": task " << parallel.rank() << " insufficient PIXBUFFER for parallel Healpix output!" << endl;
										exit(-200);
									}
									batch = bytes[j] % (PIXBUFFER * sizeof(Real));
									if (batch > 0)
										parallel.receive<char>(outbuf, batch, proc[j]-1);
								}
							}
						}
					}
				}

				for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
				{
					if (parallel.rank() == proc[j])
					{
						if ((PIXBUFFER * sizeof(Real)) - batch < 4)
						{
							memcpy((void *) (outbuf+batch), (void *) (buffer+264), (sizeof(Real) * PIXBUFFER) - batch);
							batch = PIXBUFFER * sizeof(Real);
						}
						else
						{
							memcpy((void *) (outbuf+batch), (void *) (buffer+264), 4);
							batch += 4;
						}
					}

					if (proc[j] >= 0)
					{
						bytes[j] += 4;

						if (bytes[j] / (PIXBUFFER * sizeof(Real)) > proc[j] - proc_start[j])
						{
							proc[j]++;
							if (parallel.rank() == proc[j])
							{
								if (batch > 0)
								{
									cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": task " << parallel.rank() << " insufficient PIXBUFFER for parallel Healpix output!" << endl;
									exit(-200);
								}
								batch = bytes[j] % (PIXBUFFER * sizeof(Real));
								memcpy((void *) outbuf, (void *) (buffer+268-batch), batch);
							}
						}
					}
				}
			}

			for (j = 0; j < LIGHTCONE_MAX_FIELDS; j++)
			{
				if (proc_start[j] >= 0) free(pixel[j]);

				if (bytes[j] == 0) continue;

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
				MPI_File_set_size(mapfile, (MPI_Offset) bytes[j]);

				if (parallel.rank() >= proc_start[j] && parallel.rank() <= proc[j] && batch > 0)
					MPI_File_write_at(mapfile, (MPI_Offset) ((int64_t) (parallel.rank() - proc_start[j]) * PIXBUFFER * sizeof(Real)), outbuf, batch, MPI_BYTE, &status);

				MPI_File_close(&mapfile);
			}

			free(outbuf);
#else // !HAVE_HEALPIX
			for (j = 0; j < LIGHTCONE_THICKNESS; j++)
			{
				if (sim.lightcone[i].distance[0] <= s[j+1] || sim.lightcone[i].distance[1] > s[j+1] || s[j+1] <= 0) continue;

				n = findIntersectingLightcones(sim.lightcone[i], s[j+1], s[j], domain, vertex);

				if (n > 0 && sim.out_lightcone[i] & MASK_PHI)
				{
					xsim.initialize(phi->lattice());
					xlc.initialize((*lclat[LIGHTCONE_PHI_OFFSET]));

					for (xsim.first(), xlc.first(); xsim.test(); xsim.next(), xlc.next())
					{
						for (p = 0; p < 3; p++)
							pos[p] = (double) xsim.coord(p) / (double) sim.numpts;

						for (p = 0; p < n; p++)
						{
							if (pointInShell(pos, sim.lightcone[i], s[j+1], s[j], vertex[p]))
							{
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_PHI_OFFSET+j])(xlc) = (*phi)(xsim);
								break;
							}
						}
					}
				}

				if (n > 0 && sim.out_lightcone[i] & MASK_CHI)
				{
					xsim.initialize(chi->lattice());
					xlc.initialize((*lclat[LIGHTCONE_CHI_OFFSET]));

					for (xsim.first(), xlc.first(); xsim.test(); xsim.next(), xlc.next())
					{
						for (p = 0; p < 3; p++)
							pos[p] = (double) xsim.coord(p) / (double) sim.numpts;

						for (p = 0; p < n; p++)
						{
							if (pointInShell(pos, sim.lightcone[i], s[j+1], s[j], vertex[p]))
							{
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_CHI_OFFSET+j])(xlc) = (*chi)(xsim);
								break;
							}
						}
					}
				}

				if (sim.gr_flag == 0 && sim.out_lightcone[i] & MASK_B && done_B == 0)
				{
					plan_Bi->execute(FFT_BACKWARD);
					Bi->updateHalo();
					done_B = 1;
				}

				if (n > 0 && sim.out_lightcone[i] & MASK_B)
				{
					xsim.initialize(Bi->lattice());
					xlc.initialize((*lclat[LIGHTCONE_B_OFFSET]));

					for (xsim.first(), xlc.first(); xsim.test(); xsim.next(), xlc.next())
					{
#ifdef LIGHTCONE_INTERPOLATE
						for (p = 0; p < 3; p++)
							pos[p] = (double) xsim.coord(p) / (double) sim.numpts;

						for (p = 0; p < n; p++)
						{
							if (pointInShell(pos, sim.lightcone[i], s[j+1], s[j], vertex[p]))
							{
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_B_OFFSET+3*j])(xlc) = ((*Bi)(xsim, 0) + (*Bi)(xsim-0, 0)) / (2. * a * a * sim.numpts);
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_B_OFFSET+3*j+1])(xlc) = ((*Bi)(xsim, 1) + (*Bi)(xsim-1, 1)) / (2. * a * a * sim.numpts);
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_B_OFFSET+3*j+2])(xlc) = ((*Bi)(xsim, 2) + (*Bi)(xsim-2, 2)) / (2. * a * a * sim.numpts);
								break;
							}
						}
#else
						pos[0] = (0.5 + (double) xsim.coord(0)) / (double) sim.numpts;
						pos[1] = (double) xsim.coord(1) / (double) sim.numpts;
						pos[2] = (double) xsim.coord(2) / (double) sim.numpts;

						for (p = 0; p < n; p++)
						{
							if (pointInShell(pos, sim.lightcone[i], s[j+1], s[j], vertex[p]))
							{
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_B_OFFSET+3*j])(xlc) = (*Bi)(xsim, 0) / (a * a * sim.numpts);
								break;
							}
						}

						pos[0] = (double) xsim.coord(0) / (double) sim.numpts;
						pos[1] = (0.5 + (double) xsim.coord(1)) / (double) sim.numpts;

						for (p = 0; p < n; p++)
						{
							if (pointInShell(pos, sim.lightcone[i], s[j+1], s[j], vertex[p]))
							{
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_B_OFFSET+3*j+1])(xlc) = (*Bi)(xsim, 1) / (a * a * sim.numpts);
								break;
							}
						}

						pos[1] = (double) xsim.coord(1) / (double) sim.numpts;
						pos[2] = (0.5 + (double) xsim.coord(2)) / (double) sim.numpts;

						for (p = 0; p < n; p++)
						{
							if (pointInShell(pos, sim.lightcone[i], s[j+1], s[j], vertex[p]))
							{
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_B_OFFSET+3*j+2])(xlc) = (*Bi)(xsim, 2) / (a * a * sim.numpts);
								break;
							}
						}
#endif
					}
				}

				if (sim.out_lightcone[i] & MASK_HIJ && done_hij == 0)
				{
					projectFTtensor(*SijFT, *SijFT);
					plan_Sij->execute(FFT_BACKWARD);
					Sij->updateHalo();
					done_hij = 1;
				}

				if (n > 0 && sim.out_lightcone[i] & MASK_HIJ)
				{
					xsim.initialize(Sij->lattice());
					xlc.initialize((*lclat[LIGHTCONE_HIJ_OFFSET]));

					for (xsim.first(), xlc.first(); xsim.test(); xsim.next())
					{
#ifdef LIGHTCONE_DOWNGRADE
						if ((xsim.coord(0) % LIGHTCONE_DOWNGRADE) > 0 || (xsim.coord(1) % LIGHTCONE_DOWNGRADE) > 0 || (xsim.coord(2) % LIGHTCONE_DOWNGRADE) > 0)
							continue;
#endif

						for (p = 0; p < 3; p++)
							pos[p] = (double) xsim.coord(p) / (double) sim.numpts;

#ifdef LIGHTCONE_INTERPOLATE
						for (p = 0; p < n; p++)
						{
							if (pointInShell(pos, sim.lightcone[i], s[j+1], s[j], vertex[p]))
							{
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j])(xlc) = (*Sij)(xsim, 0, 0);
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j+1])(xlc) = ((*Sij)(xsim, 0, 1) + (*Sij)(xsim-0, 0, 1) + (*Sij)(xsim-1, 0, 1) + (*Sij)(xsim-0-1, 0, 1)) / 4.;
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j+2])(xlc) = ((*Sij)(xsim, 0, 2) + (*Sij)(xsim-0, 0, 2) + (*Sij)(xsim-2, 0, 2) + (*Sij)(xsim-0-2, 0, 2)) / 4.;
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j+3])(xlc) = (*Sij)(xsim, 1, 1);
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j+4])(xlc) = ((*Sij)(xsim, 1, 2) + (*Sij)(xsim-1, 1, 2) + (*Sij)(xsim-2, 1, 2) + (*Sij)(xsim-1-2, 1, 2)) / 4.;
								break;
							}
						}
#else
#ifdef LIGHTCONE_DOWNGRADE
						for (p = 0; p < 3; p++)
							pos[p] += 0.5 / (double) sim.numpts;

						for (p = 0; p < n; p++)
						{
							if (pointInShell(pos, sim.lightcone[i], s[j+1], s[j], vertex[p]))
							{
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j])(xlc) = ((*Sij)(xsim, 0, 0) + (*Sij)(xsim+0, 0, 0) + (*Sij)(xsim+1, 0, 0) + (*Sij)(xsim+1+0, 0, 0) + (*Sij)(xsim+2, 0, 0) + (*Sij)(xsim+2+0, 0, 0) + (*Sij)(xsim+2+1, 0, 0) + (*Sij)(xsim+2+1+0, 0, 0)) / 8.;
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j+1])(xlc) = ((*Sij)(xsim, 0, 1) + (*Sij)(xsim+2, 0, 1)) / 2.;
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j+2])(xlc) = ((*Sij)(xsim, 0, 2) + (*Sij)(xsim+1, 0, 2)) / 2.;
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j+3])(xlc) = ((*Sij)(xsim, 1, 1) + (*Sij)(xsim+0, 1, 1) + (*Sij)(xsim+1, 1, 1) + (*Sij)(xsim+1+0, 1, 1) + (*Sij)(xsim+2, 1, 1) + (*Sij)(xsim+2+0, 1, 1) + (*Sij)(xsim+2+1, 1, 1) + (*Sij)(xsim+2+1+0, 1, 1)) / 8.;
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j+4])(xlc) = ((*Sij)(xsim, 1, 2) + (*Sij)(xsim+0, 1, 2)) / 2.;
								break;
							}
						}
#else
						for (p = 0; p < n; p++)
						{
							if (pointInShell(pos, sim.lightcone[i], s[j+1], s[j], vertex[p]))
							{
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j])(xlc) = (*Sij)(xsim, 0, 0);
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j+3])(xlc) = (*Sij)(xsim, 1, 1);
								break;
							}
						}

						pos[0] = (0.5 + (double) xsim.coord(0)) / (double) sim.numpts;
						pos[1] = (0.5 + (double) xsim.coord(1)) / (double) sim.numpts;

						for (p = 0; p < n; p++)
						{
							if (pointInShell(pos, sim.lightcone[i], s[j+1], s[j], vertex[p]))
							{
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j+1])(xlc) = (*Sij)(xsim, 0, 1);
								break;
							}
						}

						pos[1] = (double) xsim.coord(1) / (double) sim.numpts;
						pos[2] = (0.5 + (double) xsim.coord(2)) / (double) sim.numpts;

						for (p = 0; p < n; p++)
						{
							if (pointInShell(pos, sim.lightcone[i], s[j+1], s[j], vertex[p]))
							{
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j+2])(xlc) = (*Sij)(xsim, 0, 2);
								break;
							}
						}

						pos[0] = (double) xsim.coord(0) / (double) sim.numpts;
						pos[1] = (0.5 + (double) xsim.coord(1)) / (double) sim.numpts;

						for (p = 0; p < n; p++)
						{
							if (pointInShell(pos, sim.lightcone[i], s[j+1], s[j], vertex[p]))
							{
								(*lcbuffer[LIGHTCONE_THICKNESS*LIGHTCONE_HIJ_OFFSET+5*j+4])(xlc) = (*Sij)(xsim, 1, 2);
								break;
							}
						}
#endif
#endif
						xlc.next();
					}
				}
			}
#endif // HAVE_HEALPIX
		}
		else if (parallel.isRoot() && outfile != NULL)
			fclose(outfile);

		if (sim.out_lightcone[i] & MASK_GADGET && sim.lightcone[i].distance[0] > d - tau + 0.5 * dtau_old && sim.lightcone[i].distance[1] <= d - tau + 0.5 * dtau_old && d - tau + 0.5 * dtau_old > 0.)
		{
			n = findIntersectingLightcones(sim.lightcone[i], d - tau + (0.5 + maxvel) * dtau_old, d - tau - 0.5 * dtau, domain, vertex);

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

			pcls_cdm->saveGadget2(h5filename + filename + "_cdm", hdr, sim.lightcone[i], d - tau + (0.5 + maxvel) * dtau_old, d - tau - 0.5 * dtau, vertex, n, sim.tracer_factor[0]);

			if (sim.baryon_flag)
			{
				hdr.mass[1] = (double) sim.tracer_factor[1] * C_RHO_CRIT * cosmo.Omega_b * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[1] / GADGET_MASS_CONVERSION;
				pcls_b->saveGadget2(h5filename + filename + "_b", hdr, sim.lightcone[i], d - tau + (0.5 + maxvel) * dtau_old, d - tau - 0.5 * dtau, vertex, n, sim.tracer_factor[1]);
			}
			for (p = 0; p < cosmo.num_ncdm; p++)
			{
				if (sim.numpcl[1+sim.baryon_flag+p] == 0) continue;
				sprintf(buffer, "_ncdm%d", p);
				hdr.mass[1] = (double) sim.tracer_factor[i+1+sim.baryon_flag] * C_RHO_CRIT * cosmo.Omega_ncdm[i] * sim.boxsize * sim.boxsize * sim.boxsize / sim.numpcl[i+1+sim.baryon_flag] / GADGET_MASS_CONVERSION;
				pcls_ncdm[p].saveGadget2(h5filename + filename + buffer, hdr, sim.lightcone[i], d - tau + (0.5 + maxvel) * dtau_old, d - tau - 0.5 * dtau, vertex, n, sim.tracer_factor[i+1+sim.baryon_flag]);
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
//   pcls_cdm       pointer to (uninitialized) particle handler for CDM
//   pcls_b         pointer to (uninitialized) particle handler for baryons
//   pcls_ncdm      array of (uninitialized) particle handlers for
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
background & class_background, perturbs & class_perturbs, spectra & class_spectra, icsettings & ic,
#endif
Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_cdm, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_b, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_ncdm, Field<Real> * phi, Field<Real> * chi, Field<Real> * Bi, Field<Real> * source, Field<Real> * Sij, Field<Cplx> * scalarFT, Field<Cplx> * BiFT, Field<Cplx> * SijFT, PlanFFT<Cplx> * plan_phi, PlanFFT<Cplx> * plan_chi, PlanFFT<Cplx> * plan_Bi, PlanFFT<Cplx> * plan_source, PlanFFT<Cplx> * plan_Sij, Field<Real> * Bi_check = NULL, Field<Cplx> * BiFT_check = NULL, PlanFFT<Cplx> * plan_Bi_check = NULL)
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
			projection_T00_project(class_background, class_perturbs, class_spectra, *source, *scalarFT, plan_source, sim, ic, cosmo, fourpiG, a);
			if (sim.out_pk & MASK_DELTA)
			{
				Omega_ncdm = 0;
				for (i = 0; i < cosmo.num_ncdm; i++)
				{
					if (a < 1. / (sim.z_switch_deltancdm[i] + 1.) && cosmo.Omega_ncdm[i] > 0)
						Omega_ncdm += bg_ncdm(a, cosmo, i);
				}
				plan_source->execute(FFT_FORWARD);
				extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
				sprintf(filename, "%s%s%03d_deltaclass.dat", sim.output_path, sim.basename_pk, pkcount);
				writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * ((a < 1. / (sim.z_switch_deltarad + 1.) ? sim.radiation_flag : 0) * cosmo.Omega_rad / a + sim.fluid_flag * cosmo.Omega_fld / pow(a, 3. * cosmo.w0_fld) + Omega_ncdm) * ((a < 1. / (sim.z_switch_deltarad + 1.) ? sim.radiation_flag : 0) * cosmo.Omega_rad / a + sim.fluid_flag * cosmo.Omega_fld / pow(a, 3. * cosmo.w0_fld) + Omega_ncdm), filename, "power spectrum of delta for linear fields (CLASS)", a, sim.z_pk[pkcount]);
        for (kFT.first(); kFT.test(); kFT.next())
	           (*SijFT)(kFT, 0) = (*scalarFT)(kFT);
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

    #ifdef HAVE_CLASS
    extractCrossSpectrum(*scalarFT, *SijFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
    sprintf(filename, "%s%s%03d_deltaclass_deltam.dat", sim.output_path, sim.basename_pk, pkcount);
    writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * (cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) * ((a < 1. / (sim.z_switch_deltarad + 1.) ? sim.radiation_flag : 0) * cosmo.Omega_rad / a + sim.fluid_flag * cosmo.Omega_fld / pow(a, 3. * cosmo.w0_fld) + Omega_ncdm)  , filename, "Cross power spectrum of deltakess * delta_m", a, sim.z_pk[pkcount]);
    #endif

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

		if (cosmo.num_ncdm > 0 || sim.baryon_flag || sim.radiation_flag > 0 || sim.fluid_flag > 0)
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

	free(kbin);
	free(power);
	free(kscatter);
	free(pscatter);
	free(occupation);
}

#endif



void writeSpectra_phi_prime(metadata & sim, cosmology & cosmo, const double fourpiG, const double a, const int pkcount, Field<Real> * phi_prime ,Field<Cplx> * phi_prime_scalarFT , PlanFFT<Cplx> * phi_prime_plan)
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


// Phi_prime is dimensionful so we divide  to Hconf to make dimensionless!
if (sim.out_pk & MASK_PHI_PRIME)
  {
    phi_prime_plan->execute(FFT_FORWARD);
    extractPowerSpectrum(*phi_prime_scalarFT , kbin, power, kscatter, pscatter, occupation, sim.numbins, true, KTYPE_LINEAR);
    sprintf(filename, "%s%s%03d_phi_prime.dat", sim.output_path, sim.basename_pk, pkcount);
    writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI * (Hconf(a,fourpiG,cosmo) * Hconf(a,fourpiG,cosmo)), filename, "power spectrum of phi_prime/H (dimensionless)", a, sim.z_pk[pkcount]);
  }

free(kbin);
free(power);
free(kscatter);
free(pscatter);
free(occupation);
}
