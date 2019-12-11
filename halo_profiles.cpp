#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <cstdint>

using namespace std;

#define BUFFER 33554432

struct gadget2_header
{
  uint32_t npart[6];
  double mass[6];
  double time;
  double redshift;
  int32_t flag_sfr;
  int32_t flag_feedback;
  uint32_t npartTotal[6];
  int32_t flag_cooling;
  uint32_t num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  int32_t flag_age;
  int32_t flag_metals;
  uint32_t npartTotalHW[6];
  char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8 - 2 * 4 - 6 * 4];       /* fills to 256 Bytes */
};

int main(int argc, char **argv)
{
	char * filebase = NULL;
	char filename[100];
	char line[8192];
	FILE * snapfile;
	FILE * halofile;
    unsigned int i, j, k, blocksize1, blocksize2, num_read, batch;
	gadget2_header filehdr;
	float * posdata;
	float * halonudens;
	float * halocdmdens;
	long * haloID;
	float * halopos;
	int nhalo = 0;
	int nbins = 10;
	int species = 0;
	int snapnum = 0;
	double masslimit = 0.;
	double mass;
	float dr = 1.0;
	float dx, dy, dz, r;
	long npart, npartcdm = 0, npartnu = 0;
	
	for (i=1 ; i < argc ; i++ ){
		if ( argv[i][0] != '-' )
			continue;
		switch(argv[i][1]) {
			case 'f':
				filebase = argv[++i];
				break;
			case 'M':
				masslimit = atof(argv[++i]); //halo mass limit
				break;
			case 'b':
				nbins = atoi(argv[++i]);
				break;
			case 'd':
				dr = atof(argv[++i]);
		}
	}
	
	if (nbins < 2 || !isfinite(nbins))
	{
		cout << " error: number of bins not set properly!" << endl;
		return -1;
	}
	
	if (dr <= 0. || !isfinite(dr))
	{
		cout << " error: bin width not set properly!" << endl;
		return -1;
	}
		
	if (filebase == NULL)
	{
	    cout << " error: no file base specified!" << endl;
	    return -1;
	}
	
	line[8191] = 0;

	sprintf(filename, "%s_halos.dat", filebase);
	
	halofile = fopen(filename, "r");
	if (halofile == NULL)
	{
		cout << " error: could not open halo file!" << endl;
		return -1;
	}
	
	while (!feof(halofile) && !ferror(halofile))
	{
		fgets(line, 8192, halofile);
		if (line[8191] != 0)
		{
			cerr << " error reading halo file! Character limit (8191/line) exceeded." << endl;
		    fclose(halofile);
			return -1;
		}
		    
		if (line[0] != '#' && !feof(halofile) && !ferror(halofile))
		{
			if (sscanf(line, "%*d %*d %*f %*f %*f %*f %*f %*d %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %lf", &mass) == 1)
			{
				if (mass > masslimit)
					nhalo++;
			}
		}
	}
	
	rewind(halofile);
	
	halonudens = (float *) malloc(nhalo * nbins * sizeof(float));
	halocdmdens = (float *) malloc(nhalo * nbins * sizeof(float));
	haloID = (long *) malloc(nhalo * sizeof(long));
	halopos = (float *) malloc(nhalo * 3 * sizeof(float));
	
	for (i = 0; i < nhalo * nbins; i++)
	{
		halonudens[i] = 0.;
		halocdmdens[i] = 0.;
	}
	
	i = 0;
	
	while (!feof(halofile) && !ferror(halofile) && i < nhalo)
	{
		fgets(line, 8192, halofile);
		if (line[0] != '#' && !feof(halofile) && !ferror(halofile))
		{
			if (sscanf(line, "%ld %*d %*f %*f %*f %*f %*f %*d %f %f %f %*f %*f %*f %*f %*f %*f %*f %*f %lf", haloID+i, halopos+(3*i), halopos+(3*i)+1, halopos+(3*i)+2, &mass) == 5)
			{
				if (mass > masslimit)
					i++;
			}
		}
	}
	
	fclose(halofile);
	
	if (i < nhalo)
		nhalo = i;
		
	cout << " " << nhalo << " halos found in halo file." << endl;
	
	posdata = (float *) malloc(3 * sizeof(float) * BUFFER);

	while (true)
	{	
		sprintf(filename, "%s_cdm%d", filebase, snapnum);
		snapfile = fopen(filename, "r");
		if (snapfile == NULL && snapnum == 0)
		{
			cout << " error: could not open snapshot file!" << endl;
			return -1;
		}
		else if (snapfile == NULL)
			break;
        
		fread(&blocksize1, sizeof(int), 1, snapfile);
		if (blocksize1 != sizeof(filehdr))
		{
			cout << " error: unknown file format - header not recognized." << endl;
			fclose(snapfile);
			return -1;
		}		
    	fread(&filehdr, sizeof(filehdr), 1, snapfile);
    	fread(&blocksize2, sizeof(int), 1, snapfile);
    	if (blocksize1 != blocksize2)
    	{
			cout << " error: unknown file format - block size mismatch while reading header." << endl;
			fclose(snapfile);
			return -1;
		}
		
		cout << " reading snapshot file " << snapnum << " for cdm. HEADER information:" << endl;
		
		cout << " * number of files: " << filehdr.num_files << endl;
		cout << " * number of particles (total): ";
		for (i = 0; i < 5; i++)
			cout << filehdr.npart[i] << " (" << filehdr.npartTotal[i] << "), ";
		cout << filehdr.npart[i] << " (" << filehdr.npartTotal[i] << ")" << endl;
		cout << " * BoxSize: " << filehdr.BoxSize << endl;
		cout << " * time: " << filehdr.time << endl;
		cout << " * redshift: " << filehdr.redshift << endl;
		cout << " * Hubble parameter: " << filehdr.HubbleParam << endl;
		cout << " * Omega0: " << filehdr.Omega0 << endl;
		cout << " * OmegaLambda: " << filehdr.OmegaLambda << endl;
	
		batch = (filehdr.npart[1] < BUFFER) ? filehdr.npart[1] : BUFFER;
		fread(&blocksize1, sizeof(int), 1, snapfile);
		
		npart = 0;
		
		while (npart < filehdr.npart[1])
		{
			batch = (filehdr.npart[1] - npart < BUFFER) ? (filehdr.npart[1] - npart) : BUFFER;
			num_read = fread(posdata, sizeof(float), 3 * batch, snapfile);
			if (num_read != 3 * batch)
			{
				cout << " error: could not read position data." << endl;
				fclose(snapfile);
				return -1;
			}
	
			for (i = 0; i < batch; i++)
			{
				for (j = 0; j < nhalo; j++)
				{
					dx = posdata[3*i] - halopos[3*j]*1000.;
					if (dx > filehdr.BoxSize / 2.) dx -= filehdr.BoxSize;
					if (dx < -filehdr.BoxSize / 2.) dx += filehdr.BoxSize;
					dy = posdata[3*i+1] - halopos[3*j+1]*1000.;
					if (dy > filehdr.BoxSize / 2.) dy -= filehdr.BoxSize;
					if (dy < -filehdr.BoxSize / 2.) dy += filehdr.BoxSize;
					dz = posdata[3*i+2] - halopos[3*j+2]*1000.;
					if (dz > filehdr.BoxSize / 2.) dz -= filehdr.BoxSize;
					if (dz < -filehdr.BoxSize / 2.) dz += filehdr.BoxSize;
						
					r = sqrt(dx*dx + dy*dy + dz*dz);
						
					if (r < nbins * dr)
					{
						k = (int) floor(r / dr);
						halocdmdens[j*nbins+k] += 1.;
					}
				}
			}
		
			npart += batch;

			cout << " batch of " << batch << " particles processed (" << npart << " so far for current file)" << endl;
		}
		
		cout << " " << npart << " particles read from snapshot file." << endl;
		fclose(snapfile);

		npartcdm += npart;
		snapnum++;
	}

	mass = filehdr.mass[1];

	cout << " " << npartcdm << " cdm particles processed in total." << endl;
	
	while (true)
	{
		sprintf(filename, "%s_ncdm%d", filebase, species);
		snapfile = fopen(filename, "r");
		if (snapfile == NULL)
			break;
        
		fread(&blocksize1, sizeof(int), 1, snapfile);
		if (blocksize1 != sizeof(filehdr))
		{
			cout << " error: unknown file format - header not recognized." << endl;
			fclose(snapfile);
			return -1;
		}		
    	fread(&filehdr, sizeof(filehdr), 1, snapfile);
    	fread(&blocksize2, sizeof(int), 1, snapfile);
    	if (blocksize1 != blocksize2)
    	{
			cout << " error: unknown file format - block size mismatch while reading header." << endl;
			fclose(snapfile);
			return -1;
		}
		
		cout << " reading snapshot file for species " << species << ". HEADER information:" << endl;
		
		cout << " * number of files: " << filehdr.num_files << endl;
		cout << " * number of particles (total): ";
		for (i = 0; i < 5; i++)
			cout << filehdr.npart[i] << " (" << filehdr.npartTotal[i] << "), ";
		cout << filehdr.npart[i] << " (" << filehdr.npartTotal[i] << ")" << endl;
		cout << " * BoxSize: " << filehdr.BoxSize << endl;
		cout << " * time: " << filehdr.time << endl;
		cout << " * redshift: " << filehdr.redshift << endl;
		cout << " * Hubble parameter: " << filehdr.HubbleParam << endl;
		cout << " * Omega0: " << filehdr.Omega0 << endl;
		cout << " * OmegaLambda: " << filehdr.OmegaLambda << endl;
	
		batch = (filehdr.npart[1] < BUFFER) ? filehdr.npart[1] : BUFFER;
		fread(&blocksize1, sizeof(int), 1, snapfile);
		
		npart = 0;
		
		while (npart < filehdr.npart[1])
		{
			batch = (filehdr.npart[1] - npart < BUFFER) ? (filehdr.npart[1] - npart) : BUFFER;
			num_read = fread(posdata, sizeof(float), 3 * batch, snapfile);
			if (num_read != 3 * batch)
			{
				cout << " error: could not read position data." << endl;
				fclose(snapfile);
				return -1;
			}
			for (i = 0; i < batch; i++)
			{
				for (j = 0; j < nhalo; j++)
				{
					dx = posdata[3*i] - halopos[3*j]*1000.;
					if (dx > filehdr.BoxSize / 2.) dx -= filehdr.BoxSize;
					if (dx < -filehdr.BoxSize / 2.) dx += filehdr.BoxSize;
					dy = posdata[3*i+1] - halopos[3*j+1]*1000.;
					if (dy > filehdr.BoxSize / 2.) dy -= filehdr.BoxSize;
					if (dy < -filehdr.BoxSize / 2.) dy += filehdr.BoxSize;
					dz = posdata[3*i+2] - halopos[3*j+2]*1000.;
					if (dz > filehdr.BoxSize / 2.) dz -= filehdr.BoxSize;
					if (dz < -filehdr.BoxSize / 2.) dz += filehdr.BoxSize;
					
					r = sqrt(dx*dx + dy*dy + dz*dz);
					
					if (r < nbins * dr)
					{
						k = (int) floor(r / dr);
						halonudens[j*nbins+k] += 1.;
					}
				}
			}
		
			npart += batch;

			cout << " batch of " << batch << " particles processed (" << npart << " so far for current file)" << endl;
		}
		
		cout << " " << npart << " particles read from snapshot file." << endl;
		fclose(snapfile);
		species++;

		npartnu += npart;
	}
	
	for (i = 0; i < nhalo; i++)
	{
		for (j = 0; j < nbins; j++)
		{
			halocdmdens[i*nbins+j] /= (3 * j * (j+1) + 1) * dr * dr * dr / mass;
			halonudens[i*nbins+j] /= (3 * j * (j+1) + 1) * dr * dr * dr / filehdr.mass[1];
		}
	}
	
	sprintf(filename, "%s_halocdmdens.dat", filebase);
	halofile = fopen(filename, "w");
	
	if (halofile == NULL)
	{
		cout << " error: could not open file for output!" << endl;
	}
	else
	{
		fprintf(halofile, "ID");
		for (j = 0; j < nbins; j++)
			fprintf(halofile, " %f", (j + 0.5) * dr);
		fprintf(halofile, "\n");
		for (i = 0; i < nhalo; i++)
		{
			fprintf(halofile, " %ld", haloID[i]);
			for (j = 0; j < nbins; j++)
				fprintf(halofile, " %e", halocdmdens[i*nbins+j]);
			fprintf(halofile, "\n");
		}
		fclose(halofile);
	}

	if (species > 0)
	{	
		sprintf(filename, "%s_halonudens.dat", filebase);
		halofile = fopen(filename, "w");
	
		if (halofile == NULL)
		{
			cout << " error: could not open file for output!" << endl;
		}
		else
		{
			fprintf(halofile, "ID");
			for (j = 0; j < nbins; j++)
				fprintf(halofile, " %f", (j + 0.5) * dr);
			fprintf(halofile, "\n");
			for (i = 0; i < nhalo; i++)
			{
				fprintf(halofile, " %ld", haloID[i]);
				for (j = 0; j < nbins; j++)
					fprintf(halofile, " %e", halonudens[i*nbins+j]);
				fprintf(halofile, "\n");
			}
			fclose(halofile);
		}
	}	

	free(halonudens);
	free(halocdmdens);
	free(haloID);
	free(halopos);
	
    free(posdata);

	return 0;
}

