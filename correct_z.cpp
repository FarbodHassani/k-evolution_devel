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
  char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];       /* fills to 256 Bytes */
};

int main(int argc, char **argv)
{
	char * filename = NULL;
	char outfilename[80];
	FILE * snapfile;
	FILE * outfile;
    uint32_t i, blocksize1, blocksize2, num_read, batch;
	gadget2_header filehdr;
	float * posdata;
	float * veldata;
	long npart = 0;
	long backtrack, fastfwd;
	double dz = 0.;
	double dx_u;
	
	for (i=1 ; i < argc ; i++ ){
		if ( argv[i][0] != '-' )
			continue;
		switch(argv[i][1]) {
			case 's':
				filename = argv[++i];
				break;
		    case 'z':
				dz = atof(argv[++i]);
		}
	}
	
	if (filename == NULL)
	{
	    cout << " error: no snapshot file specified!" << endl;
	    return -1;
	}
	
	posdata = (float *) malloc(3 * sizeof(float) * BUFFER);
	veldata = (float *) malloc(3 * sizeof(float) * BUFFER);
	
	snapfile = fopen(filename, "r+");
	if (snapfile == NULL)
	{
		cout << " error: could not open snapshot file!" << endl;
		return -1;
	}
        
	cout << " reading header information..." << endl;
        
	fread(&blocksize1, sizeof(uint32_t), 1, snapfile);
	if (blocksize1 != sizeof(filehdr))
	{
		cout << " error: unknown file format - header not recognized." << endl;
		fclose(snapfile);
		return -1;
	}
	backtrack = ftell(snapfile);		
    fread(&filehdr, sizeof(filehdr), 1, snapfile);
    fread(&blocksize2, sizeof(uint32_t), 1, snapfile);
    if (blocksize1 != blocksize2)
    {
		cout << " error: unknown file format - block size mismatch while reading header." << endl;
		fclose(snapfile);
		return -1;
	}
		
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
	
	dx_u = -10.0 * dz * filehdr.time / sqrt(((filehdr.Omega0 / filehdr.time) / filehdr.time) + (filehdr.OmegaLambda * filehdr.time) + ((1. - (filehdr.Omega0 + filehdr.OmegaLambda)) / (filehdr.time * filehdr.time * filehdr.time)));
	
	cout << endl << " the dx/u for this update is: " << dx_u << endl;
	
	filehdr.redshift += dz;
	filehdr.time = 1. / (filehdr.redshift + 1.);
	fastfwd = ftell(snapfile);
	
	cout << " updating header..." << endl;
	if (fseek(snapfile, backtrack, SEEK_SET) != 0)
	{
		cout << " error: could not rewind to HEADER block." << endl;
		fclose(snapfile);
		return -1;
	}
	fwrite(&filehdr, sizeof(filehdr), 1, snapfile);
		
	cout << " reading particle data..." << endl;
	if (fseek(snapfile, fastfwd, SEEK_SET) != 0)
	{
		cout << " error: could not fast forward to POSITIONS block." << endl;
		fclose(snapfile);
		return -1;
	}
	batch = (filehdr.npart[1] < BUFFER) ? filehdr.npart[1] : BUFFER;
	fread(&blocksize1, sizeof(uint32_t), 1, snapfile);
	cout << " * size of POSITIONS block: " << blocksize1 << " bytes" << endl;
	backtrack = ftell(snapfile);
	num_read = fread(posdata, sizeof(float), 3 * batch, snapfile);
	if (num_read != 3 * batch)
	{
		cout << " error: could not read position data." << endl;
		fclose(snapfile);
		return -1;
	}
	if (fseek(snapfile, 3 * sizeof(float) * (filehdr.npart[1] - batch), SEEK_CUR) != 0)
	{
		cout << " error: could not fast forward to VELOCITIES block." << endl;
		fclose(snapfile);
		return -1;
	}
	fread(&blocksize2, sizeof(uint32_t), 1, snapfile);
	if (blocksize1 != blocksize2)
    {
		cout << " error: unknown file format - block size mismatch while reading POSITIONS." << endl;
		fclose(snapfile);
		return -1;
	}
	fread(&blocksize2, sizeof(uint32_t), 1, snapfile);
	cout << " * size of VELOCITIES block: " << blocksize2 << " bytes" << endl;
	num_read = fread(veldata, sizeof(float), 3 * batch, snapfile);
	if (num_read != 3 * batch)
	{
		cout << " error: could not read velocity data." << endl;
		fclose(snapfile);
		return -1;
	}
	fastfwd = ftell(snapfile);
	for (i = 0; i < 3 * batch; i++)
	{
		    posdata[i] += veldata[i] * dx_u;
		    if (posdata[i] >= filehdr.BoxSize) posdata[i] -= filehdr.BoxSize;
		    if (posdata[i] < 0.) posdata[i] += filehdr.BoxSize;
	}
	if (fseek(snapfile, backtrack, SEEK_SET) != 0)
	{
		cout << " error: could not rewind to POSITIONS block." << endl;
		fclose(snapfile);
		return -1;
	}
	num_read = fwrite(posdata, sizeof(float), 3 * batch, snapfile);
	if (num_read != 3 * batch)
	{
		cout << " error: could not write position data." << endl;
		fclose(snapfile);
		return -1;
	}
		
	npart += batch;
		
	while (npart < filehdr.npart[1])
	{
		batch = (filehdr.npart[1] - npart < BUFFER) ? (filehdr.npart[1] - npart) : BUFFER;
		backtrack = ftell(snapfile);
		num_read = fread(posdata, sizeof(float), 3 * batch, snapfile);
		if (num_read != 3 * batch)
		{
			cout << " error: could not read position data." << endl;
			fclose(snapfile);
			return -1;
		}
		if (fseek(snapfile, fastfwd, SEEK_SET) != 0)
		{
			cout << " error: could not fast forward to VELOCITIES block." << endl;
			fclose(snapfile);
			return -1;
		}
		num_read = fread(veldata, sizeof(float), 3 * batch, snapfile);
		if (num_read != 3 * batch)
		{
			cout << " error: could not read velocity data." << endl;
			fclose(snapfile);
			return -1;
		}
		fastfwd = ftell(snapfile);
		for (i = 0; i < 3 * batch; i++)
		{
			    posdata[i] += veldata[i] * dx_u;
			    if (posdata[i] >= filehdr.BoxSize) posdata[i] -= filehdr.BoxSize;
			    if (posdata[i] < 0.) posdata[i] += filehdr.BoxSize;
		}
		if (fseek(snapfile, backtrack, SEEK_SET) != 0)
		{
			cout << " error: could not rewind to POSITIONS block." << endl;
			fclose(snapfile);
			return -1;
		}
		num_read = fwrite(posdata, sizeof(float), 3 * batch, snapfile);
		if (num_read != 3 * batch)
		{
			cout << " error: could not write position data." << endl;
			fclose(snapfile);
			return -1;
		}
		
		npart += batch;
	}
		
	cout << " " << npart << " particles updated." << endl;
		
	fclose(snapfile);
    
    free(posdata);
    free(veldata);

	return 0;
}

