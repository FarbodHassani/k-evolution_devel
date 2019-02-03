#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

using namespace std;

#define BUFFER 131072

struct gadget2_header
{
  unsigned int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  unsigned int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];       /* fills to 256 Bytes */
};

int main(int argc, char **argv)
{
	char * filename = NULL;
	FILE * snapfile;
    unsigned int i, blocksize1, blocksize2, num_read, batch;
	gadget2_header filehdr;
	float * posdata;
	float * veldata;
	long npart = 0;
	long backtrack, fastfwd;
	double x = 0., y = 0., z = 0., d = 0.;
	
	for (i=1 ; i < argc ; i++ ){
		if ( argv[i][0] != '-' )
			continue;
		switch(argv[i][1]) {
			case 'f':
				filename = argv[++i];
				break;
			case 'x':
				x = atof(argv[++i]);
				break;
			case 'y':
				y = atof(argv[++i]);
				break;
		    case 'z':
				z = atof(argv[++i]);
				break;
			case 'd':
				d = atof(argv[++i]);
		}
	}
	
	if (filename == NULL)
	{
	    cout << " error: no snapshot file specified!" << endl;
	    return -1;
	}
	
	posdata = (float *) malloc(3 * sizeof(float) * BUFFER);
	veldata = (float *) malloc(3 * sizeof(float) * BUFFER);
	
	snapfile = fopen(filename, "r");
	if (snapfile == NULL)
	{
		cout << " error: could not open snapshot file!" << endl;
		return -1;
	}
        
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
	
	batch = (filehdr.npart[1] < BUFFER) ? filehdr.npart[1] : BUFFER;
	fread(&blocksize1, sizeof(int), 1, snapfile);
	num_read = fread(posdata, sizeof(float), 3 * batch, snapfile);
	if (num_read != 3 * batch)
	{
		cout << " error: could not read position data." << endl;
		fclose(snapfile);
		return -1;
	}
	backtrack = ftell(snapfile);
	if (fseek(snapfile, 3 * sizeof(float) * (filehdr.npart[1] - batch), SEEK_CUR) != 0)
	{
		cout << " error: could not fast forward to VELOCITIES block." << endl;
		fclose(snapfile);
		return -1;
	}
	fread(&blocksize2, sizeof(int), 1, snapfile);
	if (blocksize1 != blocksize2)
    {
		cout << " error: unknown file format - block size mismatch while reading POSITIONS." << endl;
		fclose(snapfile);
		return -1;
	}
	fread(&blocksize2, sizeof(int), 1, snapfile);
	num_read = fread(veldata, sizeof(float), 3 * batch, snapfile);
	if (num_read != 3 * batch)
	{
		cout << " error: could not read velocity data." << endl;
		fclose(snapfile);
		return -1;
	}
	fastfwd = ftell(snapfile);
	for (i = 0; i < batch; i++)
	{
		    if (((posdata[3*i] > x-d && posdata[3*i] < x+d) || (posdata[3*i] > filehdr.BoxSize+x-d) || (posdata[3*i] < x+d-filehdr.BoxSize))
		    	&& ((posdata[3*i+1] > y-d && posdata[3*i+1] < y+d) || (posdata[3*i+1] > filehdr.BoxSize+y-d) || (posdata[3*i+1] < y+d-filehdr.BoxSize))
		    	&& ((posdata[3*i+2] > z-d && posdata[3*i+2] < z+d) || (posdata[3*i+2] > filehdr.BoxSize+z-d) || (posdata[3*i+2] < z+d-filehdr.BoxSize)))
		    	cout << posdata[3*i] << "  " << posdata[3*i+1] << "  " << posdata[3*i+2] << "  " << veldata[3*i] << "  " << veldata[3*i+1] << "  " << veldata[3*i+2] << endl;
	}
	if (fseek(snapfile, backtrack, SEEK_SET) != 0)
	{
		cout << " error: could not rewind to POSITIONS block." << endl;
		fclose(snapfile);
		return -1;
	}
		
	npart += batch;
		
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
		backtrack = ftell(snapfile);
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
		for (i = 0; i < batch; i++)
		{
			    if (((posdata[3*i] > x-d && posdata[3*i] < x+d) || (posdata[3*i] > filehdr.BoxSize+x-d) || (posdata[3*i] < x+d-filehdr.BoxSize))
		    		&& ((posdata[3*i+1] > y-d && posdata[3*i+1] < y+d) || (posdata[3*i+1] > filehdr.BoxSize+y-d) || (posdata[3*i+1] < y+d-filehdr.BoxSize))
		    		&& ((posdata[3*i+2] > z-d && posdata[3*i+2] < z+d) || (posdata[3*i+2] > filehdr.BoxSize+z-d) || (posdata[3*i+2] < z+d-filehdr.BoxSize)))
		    		cout << posdata[3*i] << "  " << posdata[3*i+1] << "  " << posdata[3*i+2] << "  " << veldata[3*i] << "  " << veldata[3*i+1] << "  " << veldata[3*i+2] << endl;
		}
		if (fseek(snapfile, backtrack, SEEK_SET) != 0)
		{
			cout << " error: could not rewind to POSITIONS block." << endl;
			fclose(snapfile);
			return -1;
		}
		
		npart += batch;
	}
		
	fclose(snapfile);
    
    free(posdata);
    free(veldata);

	return 0;
}

