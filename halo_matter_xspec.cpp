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

#define CHUNKSIZE 33554432

int project_halos(Field<Real> * halodens, float * posdata, const int count, const double boxsize);
int project_pcls(Field<Real> * matterdens, float * posdata, float * veldata, const int count, const double a, const double mass, const double boxsize);

int main(int argc, char **argv)
{
	char * filebase = NULL;
	int n = 0;
	int m = 0;
	int i, tot = 0, mcount = 0;
	int metaread = 0;
	int numpts = 0;
	int numbins = 0;
	icsettings ic;
	gadget2_header hdr;
	int box[3];
	long nhalos = 0, npcls, count = 0;
	long numpts3d;
	double h = P_HUBBLE;
	double boxsize;
	double a;
	double masslimit = 0.;
	double masslimit2 = 1.0e20;
	gsl_spline * phispline = NULL;
	gsl_spline * dummyspline = NULL;
	float * posdata;
	float * veldata;
	float * massdata;
	long * iddata;
	FILE * fileptr;
	FILE * outfile;
	char filename[1024];
	char line[32768];
	uint32_t blocksize;
	int64_t fastforward, backtrack;
	double renorm = 0.;
	double mass;
//	string fnstr;
//	double Omega_m;
	
	//fnstr.reserve(PARAM_MAX_LENGTH);
	//hdr.npart[1] = 0;

	Real * kbin;
	Real * power;
	Real * kscatter;
	Real * pscatter;
	int * occupation;

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
			case 'l':
				masslimit = atof(argv[++i]); //halo mass limit
				break;
			case 'x':
				masslimit2 = atof(argv[++i]);
				break;
			case 's':
				ic.seed = atoi(argv[++i]); //seed
				break;
			case 'b':
				numbins = atoi(argv[++i]); //number of Pk bins
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

	COUT << " Halo-matter cross spectrum tool" << endl;
	COUT << " Ngrid = " << numpts << "; Pk bins = " << numbins << "; halo mass limit = " << masslimit << ", " << masslimit2 << "; seed = " << ic.seed << endl << endl;

	box[0] = numpts;
	box[1] = numpts;
	box[2] = numpts;

	numpts3d = (long) numpts * (long) numpts * (long) numpts;

	ic.A_s = P_SPECTRAL_AMP;
	ic.n_s = P_SPECTRAL_INDEX;
	ic.k_pivot = P_PIVOT_SCALE;

	Lattice lat(3,box,1);
	Lattice latFT;
	latFT.initializeRealFFT(lat,0);

	Site x(lat);
#ifdef READ_PHIFILE
    int box2[3];
    box2[0] = 1024;
    box2[1] = 1024;
    box2[2] = 1024; 
    Lattice lat2(3,box2,1);
    Site x2(lat2);
    int xp, yp, zp;
    Real temp;
#else
	double * temp;
#endif


	Field<Real> halodens;
	Field<Real> matterdens;
	Field<Real> phi;
	Field<Cplx> scalarFT;
	Field<Cplx> scalarFT2;

	halodens.initialize(lat,1);
	matterdens.initialize(lat,1);
	scalarFT.initialize(latFT,1);
	scalarFT2.initialize(latFT,1);
	PlanFFT<Cplx> plan_halodens(&halodens, &scalarFT);
	PlanFFT<Cplx> plan_matterdens(&matterdens, &scalarFT2);
#ifdef READ_PHIFILE
	phi.initialize(lat2,1);
	phi.alloc();
#else
	phi.initialize(lat,1);
	PlanFFT<Cplx> plan_phi(&phi, &scalarFT);
#endif

/*	Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> pcls_cdm;
	Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> pcls_ncdm;
	
	part_simple_info pcls_cdm_info;
	part_simple_dataType pcls_cdm_dataType;
	part_simple_info pcls_ncdm_info;
	part_simple_dataType pcls_ncdm_dataType;
	
	Real boxSize[3] = {1.,1.,1.};
	
	strcpy(pcls_cdm_info.type_name, "part_simple");
	pcls_cdm_info.mass = 0.;
	pcls_cdm_info.relativistic = false;
	
	strcpy(pcls_ncdm_info.type_name, "part_simple");
	pcls_ncdm_info.mass = 0.;
	pcls_ncdm_info.relativistic = true;
	
	
	pcls_cdm.initialize(pcls_cdm_info, pcls_cdm_dataType, &(matterdens.lattice()), boxSize);
	pcls_ncdm.initialize(pcls_ncdm_info, pcls_ncdm_dataType, &(matterdens.lattice()), boxSize); */

	posdata = (float *) malloc(3 * sizeof(float) * CHUNKSIZE);
	veldata = (float *) malloc(3 * sizeof(float) * CHUNKSIZE);

	projection_init(&halodens);
	projection_init(&matterdens);

	if (parallel.isRoot())
	{
		sprintf(filename, "%s_halos.dat", filebase);

		fileptr = fopen(filename, "r");

		massdata = (float *) malloc(sizeof(float) * CHUNKSIZE);
		iddata = (long *) malloc(sizeof(long) * CHUNKSIZE);

		if (fileptr == NULL)
		{
			cout << " error: could not open halo file!" << endl;
			parallel.abortForce();
		}

		while (!feof(fileptr) && !ferror(fileptr))
		{
			if(fgets(line, 32768, fileptr) == NULL) break;

			if (line[0] == '\0') continue;

			if (line[0] == '#')
			{
				if (sscanf(line, "#a =  %lf", &a) == 1)
				{
					COUT << " a = " << a << endl;
					metaread++;
				}
				else if (sscanf(line, "#Om = %*f; Ol = %*f; h = %lf", &h) == 1)
				{
					COUT << " h = " << h << endl;
					metaread++;
				}
				else if (sscanf(line, "#Box size: %lf", &boxsize) == 1)
				{
					COUT << " boxsize = " << boxsize << "Mpc/h" << endl;
					metaread++;
				}
				
				if (metaread == 3)
				{
					parallel.broadcast<double>(a, 0);
					parallel.broadcast<double>(h, 0);
					parallel.broadcast<double>(boxsize, 0);
					metaread++;
				}
			}
			else
			{
				if (metaread != 4)
				{
					cout << " error: metadata (a, h, boxsize) could not be read from halo file!" << endl;
					parallel.abortForce();
				}

				if (sscanf(line, "%ld %*d %*f %*f %*f %*f %*f %*d %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %f", iddata+mcount, posdata+(3*count), posdata+(3*count+1), posdata+(3*count+2), massdata+mcount) != 5)
				{
					cout << " error: could not parse position data for halo number " << count << endl;
					parallel.abortForce();
				}

				if (massdata[mcount] > masslimit && massdata[mcount] < masslimit2)
					count++;

				mcount++;

				if (mcount == CHUNKSIZE)
				{
#ifdef WRITE_MASSFILE
					sprintf(filename, "%s_halomasses.dat", filebase);

					outfile = fopen(filename, "a");
					
					if (outfile != NULL)
					{
						for (i = 0; i < mcount; i++)
							fprintf(outfile, " %ld %f\n", iddata[i], massdata[i]);
						fclose(outfile);
					}
#endif
					mcount = 0;
				}

				if (count == CHUNKSIZE)
				{
					parallel.broadcast(count, 0);
					parallel.broadcast<float>(posdata, 3*count, 0);
					tot += project_halos(&halodens, posdata, count, boxsize);
					count = 0;
				}
			}
		}

#ifdef WRITE_MASSFILE
		if (mcount > 0)
		{
			sprintf(filename, "%s_halomasses.dat", filebase);

			outfile = fopen(filename, "a");
					
			if (outfile != NULL)
			{
				for (i = 0; i < mcount; i++)
					fprintf(outfile, " %ld %e\n", iddata[i], massdata[i]);
				fclose(outfile);
			}
		}
#endif

		free(massdata);
		free(iddata);

		if (count > 0)
		{
			parallel.broadcast(count, 0);
			parallel.broadcast<float>(posdata, 3*count, 0);
			tot += project_halos(&halodens, posdata, count, boxsize);
			count = 0;
		}
		parallel.broadcast(count, 0);

		COUT << " " << tot << " halos used." << endl;
		
		for (i = 0; i < 8; i++)
		{
			COUT << " reading cdm file " << i << endl;
			sprintf(filename, "%s_cdm%d", filebase, i);

			fileptr = fopen(filename, "r");

			if (fileptr == NULL)
			{
				cout << " error: could not open cdm file!" << endl;
				parallel.abortForce();
			}
			
			fread(&blocksize, sizeof(uint32_t), 1, fileptr);
			fread(&hdr, sizeof(hdr), 1, fileptr);
			fread(&blocksize, sizeof(uint32_t), 1, fileptr);
			fread(&blocksize, sizeof(uint32_t), 1, fileptr);
			npcls = 0;
			fastforward = ftell(fileptr) + 2l * sizeof(uint32_t) + 3l * sizeof(float) * hdr.npart[1];
			mass = hdr.mass[1] / 8;
			boxsize = hdr.BoxSize;
			
			COUT << " header read successfully. npart = " << hdr.npart[1] << ", mass = " << mass << ", boxsize = " << boxsize << endl;
			
			parallel.broadcast<double>(mass, 0);
			parallel.broadcast<double>(boxsize, 0);
			
			while (npcls < hdr.npart[1])
			{
				count = (hdr.npart[1] - npcls > CHUNKSIZE) ? CHUNKSIZE : (hdr.npart[1] - npcls);
				fread(posdata, sizeof(float), 3 * count, fileptr);
				backtrack = ftell(fileptr);
				fseek(fileptr, fastforward, SEEK_SET);
				fread(veldata, sizeof(float), 3 * count, fileptr);
				fastforward = ftell(fileptr);
				fseek(fileptr, backtrack, SEEK_SET);
				parallel.broadcast(count, 0);
				parallel.broadcast<float>(posdata, 3*count, 0);
				parallel.broadcast<float>(veldata, 3*count, 0);
				renorm += count * mass;
				npcls += count;
				project_pcls(&matterdens, posdata, veldata, count, a, mass, boxsize);
				//COUT << " " << npcls << " particles..." << endl;
			}
			
			fclose(fileptr);
			
			count = 0;
			parallel.broadcast(count, 0);		
			COUT << " " << npcls << " particles projected." << endl;
		}
		
		for (i = 0; i < 4; i++)
		{
			COUT << " reading ncdm file " << i << endl;
			sprintf(filename, "%s_ncdm%d", filebase, i);

			fileptr = fopen(filename, "r");

			if (fileptr == NULL)
			{
				cout << " error: could not open ncdm file!" << endl;
				parallel.abortForce();
			}
			
			fread(&blocksize, sizeof(uint32_t), 1, fileptr);
			fread(&hdr, sizeof(hdr), 1, fileptr);
			fread(&blocksize, sizeof(uint32_t), 1, fileptr);
			fread(&blocksize, sizeof(uint32_t), 1, fileptr);
			npcls = 0;
			fastforward = ftell(fileptr) + 2l * sizeof(uint32_t) + 3l * sizeof(float) * hdr.npart[1];
			mass = hdr.mass[1];
			boxsize = hdr.BoxSize;
			
			COUT << " header read successfully. npart = " << hdr.npart[1] << ", mass = " << mass << ", boxsize = " << boxsize << endl;
			
			parallel.broadcast<double>(mass, 0);
			parallel.broadcast<double>(boxsize, 0);
			
			while (npcls < hdr.npart[1])
			{
				count = (hdr.npart[1] - npcls > CHUNKSIZE) ? CHUNKSIZE : (hdr.npart[1] - npcls);
				fread(posdata, sizeof(float), 3 * count, fileptr);
				backtrack = ftell(fileptr);
				fseek(fileptr, fastforward, SEEK_SET);
				fread(veldata, sizeof(float), 3 * count, fileptr);
				fastforward = ftell(fileptr);
				fseek(fileptr, backtrack, SEEK_SET);
				parallel.broadcast(count, 0);
				parallel.broadcast<float>(posdata, 3*count, 0);
				parallel.broadcast<float>(veldata, 3*count, 0);
				renorm += count * mass;
				npcls += count;
				project_pcls(&matterdens, posdata, veldata, count, a, mass, boxsize);
			}
			
			fclose(fileptr);
			
			count = 0;
			parallel.broadcast(count, 0);	
			COUT << " " << npcls << " particles projected." << endl;
		}
	}
	else
	{
		parallel.broadcast<double>(a, 0);
		parallel.broadcast<double>(h, 0);
		parallel.broadcast<double>(boxsize, 0);

		do
		{
			parallel.broadcast(count, 0);
			if (count > 0)
			{
				parallel.broadcast<float>(posdata, 3*count, 0);
				tot += project_halos(&halodens, posdata, count, boxsize);
			}
		}
		while (count > 0);
		
		for (i = 0; i < 12; i++)
		{
			/*if (parallel.rank() == 5)
			{
				cout << " proc#5: expecting data for file " << i << "..." << endl;
			}*/
			parallel.broadcast<double>(mass, 0);
			parallel.broadcast<double>(boxsize, 0);
			/*if (parallel.rank() == 5)
			{
				cout << " proc#5: received mass = " << mass << ", boxsize = " << boxsize << endl;
			}*/
			do
			{
				parallel.broadcast(count, 0);
				if (count > 0)
				{
					parallel.broadcast<float>(posdata, 3*count, 0);
					parallel.broadcast<float>(veldata, 3*count, 0);
					renorm += count * mass;
					project_pcls(&matterdens, posdata, veldata, count, a, mass, boxsize);
					/*if (parallel.rank() == 5)
					{
						cout << " proc#5: received " << count << " particles." << endl;
					}*/
				}
			}
			while (count > 0);
		}
	}

	projection_T00_comm(&halodens);

	free(posdata);
	free(veldata);
	
/*	for (i = 0; i < 8; i++)
	{
		COUT << " reading cdm file " << i << endl;
		sprintf(filename, "%s_cdm%d", filebase, i);
		fnstr.assign(filename);
		pcls_cdm.loadGadget2(fnstr, hdr);
		if (i == 0) Omega_m = hdr.Omega0;
		if (hdr.npart[1] == 0) break;
		numcdm += hdr.npart[1];
	}
	
	pcls_cdm.parts_info()->mass = (Omega_m - 0.2 / P_NCDM_MASS_OMEGA / h / h) / (Real) numcdm;

	COUT << " " << numcdm << " cdm particles read." << endl;
	
	projection_T00_project(&pcls_cdm, &matterdens);
	
	for (i = 0; i < 4; i++)
	{
		COUT << " reading ncdm file " << i << endl;
		sprintf(filename, "%s_ncdm%d", filebase, i);
		fnstr.assign(filename);
		pcls_ncdm.loadGadget2(fnstr, hdr);
		if (hdr.npart[1] == 0) break;
		numncdm += hdr.npart[1];
	}
	
	pcls_ncdm.parts_info()->mass = (0.2 / P_NCDM_MASS_OMEGA / h / h) / (Real) numncdm;

	COUT << " " << numncdm << " ncdm particles read." << endl;
	
	projection_T00_project(&pcls_ncdm, &matterdens); */
	
	projection_T00_comm(&matterdens);

#ifdef READ_PHIFILE
	string h5filename(filebase);
	phi.loadHDF5(h5filename + "_phi.h5");
	phi.updateHalo();

	COUT << " phi file read successfully." << endl;

	for (x.first(); x.test(); x.next())
	{
		temp = 0.;
		
		x2.setCoord(x.coord(0)/4, x.coord(1)/4, x.coord(2)/4);
		
		xp = x.coord(0) % 4;
		yp = x.coord(1) % 4;
		zp = x.coord(2) % 4;
		
		if (xp > 1)
		{
			if (yp > 1)
			{
				if (zp > 1)
				{
					temp = (5.5 - xp) * (5.5 - yp) * (5.5 - zp) * phi(x2);
					temp += (xp - 1.5) * (5.5 - yp) * (5.5 - zp) * phi(x2+0);
					temp += (5.5 - xp) * (yp - 1.5) * (5.5 - zp) * phi(x2+1);
					temp += (xp - 1.5) * (yp - 1.5) * (5.5 - zp) * phi(x2+0+1);
					temp += (5.5 - xp) * (5.5 - yp) * (zp - 1.5) * phi(x2+2);
					temp += (xp - 1.5) * (5.5 - yp) * (zp - 1.5) * phi(x2+0+2);
					temp += (5.5 - xp) * (yp - 1.5) * (zp - 1.5) * phi(x2+1+2);
					temp += (xp - 1.5) * (yp - 1.5) * (zp - 1.5) * phi(x2+0+1+2);
				}
				else
				{
					temp = (5.5 - xp) * (5.5 - yp) * (1.5 - zp) * phi(x2-2);
					temp += (xp - 1.5) * (5.5 - yp) * (1.5 - zp) * phi(x2+0-2);
					temp += (5.5 - xp) * (yp - 1.5) * (1.5 - zp) * phi(x2+1-2);
					temp += (xp - 1.5) * (yp - 1.5) * (1.5 - zp) * phi(x2+0+1-2);
					temp += (5.5 - xp) * (5.5 - yp) * (zp + 2.5) * phi(x2);
					temp += (xp - 1.5) * (5.5 - yp) * (zp + 2.5) * phi(x2+0);
					temp += (5.5 - xp) * (yp - 1.5) * (zp + 2.5) * phi(x2+1);
					temp += (xp - 1.5) * (yp - 1.5) * (zp + 2.5) * phi(x2+0+1);
				}
			}
			else
			{
				if (zp > 1)
				{
					temp = (5.5 - xp) * (1.5 - yp) * (5.5 - zp) * phi(x2-1);
					temp += (xp - 1.5) * (1.5 - yp) * (5.5 - zp) * phi(x2+0-1);
					temp += (5.5 - xp) * (yp + 2.5) * (5.5 - zp) * phi(x2);
					temp += (xp - 1.5) * (yp + 2.5) * (5.5 - zp) * phi(x2+0);
					temp += (5.5 - xp) * (1.5 - yp) * (zp - 1.5) * phi(x2+2-1);
					temp += (xp - 1.5) * (1.5 - yp) * (zp - 1.5) * phi(x2+0+2-1);
					temp += (5.5 - xp) * (yp + 2.5) * (zp - 1.5) * phi(x2+2);
					temp += (xp - 1.5) * (yp + 2.5) * (zp - 1.5) * phi(x2+0+2);
				}
				else
				{
					temp = (5.5 - xp) * (1.5 - yp) * (1.5 - zp) * phi(x2-1-2);
					temp += (xp - 1.5) * (1.5 - yp) * (1.5 - zp) * phi(x2+0-1-2);
					temp += (5.5 - xp) * (yp + 2.5) * (1.5 - zp) * phi(x2-2);
					temp += (xp - 1.5) * (yp + 2.5) * (1.5 - zp) * phi(x2+0-2);
					temp += (5.5 - xp) * (1.5 - yp) * (zp + 2.5) * phi(x2-1);
					temp += (xp - 1.5) * (1.5 - yp) * (zp + 2.5) * phi(x2+0-1);
					temp += (5.5 - xp) * (yp + 2.5) * (zp + 2.5) * phi(x2);
					temp += (xp - 1.5) * (yp + 2.5) * (zp + 2.5) * phi(x2+0);
				}
			}
		}
		else
		{
			if (yp > 1)
			{
				if (zp > 1)
				{
					temp = (1.5 - xp) * (5.5 - yp) * (5.5 - zp) * phi(x2-0);
					temp += (xp + 2.5) * (5.5 - yp) * (5.5 - zp) * phi(x2);
					temp += (1.5 - xp) * (yp - 1.5) * (5.5 - zp) * phi(x2-0+1);
					temp += (xp + 2.5) * (yp - 1.5) * (5.5 - zp) * phi(x2+1);
					temp += (1.5 - xp) * (5.5 - yp) * (zp - 1.5) * phi(x2-0+2);
					temp += (xp + 2.5) * (5.5 - yp) * (zp - 1.5) * phi(x2+2);
					temp += (1.5 - xp) * (yp - 1.5) * (zp - 1.5) * phi(x2-0+1+2);
					temp += (xp + 2.5) * (yp - 1.5) * (zp - 1.5) * phi(x2+1+2);
				}
				else
				{
					temp = (1.5 - xp) * (5.5 - yp) * (1.5 - zp) * phi(x2-0-2);
					temp += (xp + 2.5) * (5.5 - yp) * (1.5 - zp) * phi(x2-2);
					temp += (1.5 - xp) * (yp - 1.5) * (1.5 - zp) * phi(x2-0+1-2);
					temp += (xp + 2.5) * (yp - 1.5) * (1.5 - zp) * phi(x2+1-2);
					temp += (1.5 - xp) * (5.5 - yp) * (zp + 2.5) * phi(x2-0);
					temp += (xp + 2.5) * (5.5 - yp) * (zp + 2.5) * phi(x2);
					temp += (1.5 - xp) * (yp - 1.5) * (zp + 2.5) * phi(x2-0+1);
					temp += (xp + 2.5) * (yp - 1.5) * (zp + 2.5) * phi(x2+1);
				}
			}
			else
			{
				if (zp > 1)
				{
					temp = (1.5 - xp) * (1.5 - yp) * (5.5 - zp) * phi(x2-0-1);
					temp += (xp + 2.5) * (1.5 - yp) * (5.5 - zp) * phi(x2-1);
					temp += (1.5 - xp) * (yp + 2.5) * (5.5 - zp) * phi(x2-0);
					temp += (xp + 2.5) * (yp + 2.5) * (5.5 - zp) * phi(x2);
					temp += (1.5 - xp) * (1.5 - yp) * (zp - 1.5) * phi(x2-0+2-1);
					temp += (xp + 2.5) * (1.5 - yp) * (zp - 1.5) * phi(x2+2-1);
					temp += (1.5 - xp) * (yp + 2.5) * (zp - 1.5) * phi(x2-0+2);
					temp += (xp + 2.5) * (yp + 2.5) * (zp - 1.5) * phi(x2+2);
				}
				else
				{
					temp = (1.5 - xp) * (1.5 - yp) * (1.5 - zp) * phi(x2-0-1-2);
					temp += (xp + 2.5) * (1.5 - yp) * (1.5 - zp) * phi(x2-1-2);
					temp += (1.5 - xp) * (yp + 2.5) * (1.5 - zp) * phi(x2-0-2);
					temp += (xp + 2.5) * (yp + 2.5) * (1.5 - zp) * phi(x2-2);
					temp += (1.5 - xp) * (1.5 - yp) * (zp + 2.5) * phi(x2-0-1);
					temp += (xp + 2.5) * (1.5 - yp) * (zp + 2.5) * phi(x2-1);
					temp += (1.5 - xp) * (yp + 2.5) * (zp + 2.5) * phi(x2-0);
					temp += (xp + 2.5) * (yp + 2.5) * (zp + 2.5) * phi(x2);
				}
			}
		}
	
		halodens(x) *= (1. + 3. * temp / 64.) / (Real) tot;
		matterdens(x) *= (1. + 3. * temp / 64.) / renorm;
	}
#else
	sprintf(filename, "%s_tk.dat", filebase);

	loadTransferFunctions(filename, phispline, dummyspline, "phi", boxsize, h);

	temp = (double *) malloc(phispline->size * sizeof(double));

	for (i = 0; i < phispline->size; i++)
		temp[i] = -phispline->y[i] * M_PI * sqrt(Pk_primordial(phispline->x[i] * h / boxsize, ic) / phispline->x[i]) / phispline->x[i];
	
	gsl_spline_free(phispline);
	phispline = gsl_spline_alloc(gsl_interp_cspline, dummyspline->size);
	gsl_spline_init(phispline, dummyspline->x, temp, dummyspline->size);
	gsl_spline_free(dummyspline);

	generateRealization(scalarFT, 0., phispline, (unsigned int) ic.seed, ICFLAG_KSPHERE, 0);
	
	gsl_spline_free(phispline);
	free(temp);

	plan_phi.execute(FFT_BACKWARD);
	
#ifdef WRITE_PHIFILE
	string h5filename(filebase);
	phi.saveHDF5_coarseGrain3D(h5filename + "_linearphi.h5", 4);
#endif

	for (x.first(); x.test(); x.next())
	{
		halodens(x) *= (1. + 3. * phi(x)) / (Real) tot;
		matterdens(x) *= (1. + 3. * phi(x)) / Omega_m;
	}
#endif

	plan_halodens.execute(FFT_FORWARD);
	plan_matterdens.execute(FFT_FORWARD);

	kbin = (Real *) malloc(numbins * sizeof(Real));
	power = (Real *) malloc(numbins * sizeof(Real));
	kscatter = (Real *) malloc(numbins * sizeof(Real));
	pscatter = (Real *) malloc(numbins * sizeof(Real));
	occupation = (int *) malloc(numbins * sizeof(int));

	extractCrossSpectrum(scalarFT, scalarFT2, kbin, power, kscatter, pscatter, occupation, numbins, true, KTYPE_LINEAR);

	sprintf(filename, "%s_haloxmatter_pk.dat", filebase);
	writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, numbins, boxsize, 2. * M_PI * M_PI, filename, "cross spectrum of halos x matter", a);

	free(kbin);
	free(power);
	free(kscatter);
	free(pscatter);
	free(occupation);
}

int project_halos(Field<Real> * halodens, float * posdata, const int count, const double boxsize)
{
	Site x(halodens->lattice());

	int i, xh, yh, zh, numpts;
	double xr, yr, zr, dummy;

	numpts = halodens->lattice().size(0);

	for (i = 0; i < count; i++)
	{
		xr = modf((double) numpts * (double) posdata[3*i]/boxsize, &dummy);
		xh = (int) dummy;
		yr = modf((double) numpts * (double) posdata[3*i+1]/boxsize, &dummy);
		yh = (int) dummy;
		zr = modf((double) numpts * (double) posdata[3*i+2]/boxsize, &dummy);
		zh = (int) dummy;
		
		if (x.setCoord(xh, yh, zh))
		{
			(*halodens)(x) += (1.-xr) * (1.-yr) * (1.-zr);
			(*halodens)(x+0) += xr * (1.-yr) * (1.-zr);
			(*halodens)(x+1) += (1.-xr) * yr * (1.-zr);
			(*halodens)(x+2) += (1.-xr) * (1.-yr) * zr;
			(*halodens)(x+0+1) += xr * yr * (1.-zr);
			(*halodens)(x+0+2) += xr * (1.-yr) * zr;
			(*halodens)(x+1+2) += (1.-xr) * yr * zr;
			(*halodens)(x+0+1+2) += xr * yr * zr;
		}
	}

	return count;
}

int project_pcls(Field<Real> * matterdens, float * posdata, float * veldata, const int count, const double a, const double mass, const double boxsize)
{
	Site x(matterdens->lattice());

	int i, xh, yh, zh, numpts;
	double xr, yr, zr, dummy, e;

	numpts = matterdens->lattice().size(0);

	for (i = 0; i < count; i++)
	{
		xr = modf((double) numpts * (double) posdata[3*i]/boxsize, &dummy);
		xh = (int) dummy;
		yr = modf((double) numpts * (double) posdata[3*i+1]/boxsize, &dummy);
		yh = (int) dummy;
		zr = modf((double) numpts * (double) posdata[3*i+2]/boxsize, &dummy);
		zh = (int) dummy;
		
		e = mass * a * sqrt(veldata[3*i]*veldata[3*i]*GADGET_VELOCITY_CONVERSION*GADGET_VELOCITY_CONVERSION*a + veldata[3*i+1]*veldata[3*i+1]*GADGET_VELOCITY_CONVERSION*GADGET_VELOCITY_CONVERSION*a + veldata[3*i+2]*veldata[3*i+2]*GADGET_VELOCITY_CONVERSION*GADGET_VELOCITY_CONVERSION*a + 1);
		
		if (x.setCoord(xh, yh, zh))
		{
			(*matterdens)(x) += (1.-xr) * (1.-yr) * (1.-zr) * e;
			(*matterdens)(x+0) += xr * (1.-yr) * (1.-zr) * e;
			(*matterdens)(x+1) += (1.-xr) * yr * (1.-zr) * e;
			(*matterdens)(x+2) += (1.-xr) * (1.-yr) * zr * e;
			(*matterdens)(x+0+1) += xr * yr * (1.-zr) * e;
			(*matterdens)(x+0+2) += xr * (1.-yr) * zr * e;
			(*matterdens)(x+1+2) += (1.-xr) * yr * zr * e;
			(*matterdens)(x+0+1+2) += xr * yr * zr * e;
		}
	}

	return count;
}



