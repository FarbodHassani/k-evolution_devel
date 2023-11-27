#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "metadata.hpp"

using namespace std;

int main(int argc, char **argv)
{
	char * iname;
	char * oname;
	char filename[1024];
	char ofilename[1024];
	FILE * infile;
	FILE * outfile;

	uint64_t numpart_tot = 0;
	uint64_t numpart_write = 0;
	uint32_t blocksize = 0;
	int numfiles = 1;
	int numfiles2 = 1;
	long numread = 0;
	long numwrite = 0;
	gadget2_header hdr;
	gadget2_header outhdr;

	long backtrack;
	long fastforward;
	float * posbatch = NULL;
	float * velbatch = NULL;
	uint32_t batch;
#if GADGET_ID_BYTES == 8
	uint64_t * IDbatch = NULL;
#else
	uint32_t * IDbatch = NULL;
#endif

	double offset = 0.;
	
	if (argc < 2)
	{
		cout << " redistributes Gadget-2 binaries" << endl << endl;
		
		cout << " List of command-line options:" << endl;
		cout << " -i <filebase>       : input file base (mandatory)" << endl;
		cout << " -o <filebase>       : output file base (mandatory)" << endl;
		cout << " -n <numfiles>       : number of Gadget-2 binaries to distribute the" << endl;
		cout << "                       output over (optional, default 1)" << endl;
		cout << " The output will be written to <numfiles> approximately" << endl;
		cout << " equal-sized Gadget-2 binaries." << endl;
		return 0;
	}

	for (int i = 1 ; i < argc ; i++ ){
		if ( argv[i][0] != '-' )
			continue;
		switch(argv[i][1]) {
			case 'n':
				numfiles = atoi(argv[++i]); // number of output files
				break;
			case 'o':
				oname = argv[++i];
				break;
			case 'i':
				iname = argv[++i];
				break;
			default:
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unknown command-line parameter " << argv[i] << endl << " call redistribute without arguments to display help" << endl;
				return -1;
		}
	}

	cout << " reading particle header..." << endl << endl;
	
	sprintf(filename, "%s.0", iname);

	infile = fopen(filename, "rb");

	if (infile != NULL)
	{
		fread(&blocksize, sizeof(uint32_t), 1, infile);

		if (blocksize != sizeof(hdr))
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unknown file format " << filename << "!" << endl;
			fclose(infile);
		}

		if(fread(&hdr, sizeof(hdr), 1, infile) != 1)
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read header block from " << filename << "!" << endl;
			fclose(infile);
		}

		numpart_tot += (uint64_t) hdr.npartTotal[1] + ((uint64_t) hdr.npartTotalHW[1] << 32);
		fclose(infile);
	}
	else
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to open " << filename << "!" << endl;
		return -1;
	}

	cout << " particle header read successfully. Total number of particles = " << numpart_tot << endl << endl;

	for (int i = 0; i < 6; i++)
	{
		outhdr.npart[i] = 0;
		outhdr.mass[i] = 0.;
		outhdr.npartTotal[i] = hdr.npartTotal[i];
		outhdr.npartTotalHW[i] = hdr.npartTotalHW[i];
	}
	for (int i = 0; i < 256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8 - 2 * 4 - 6 * 4; i++)
		outhdr.fill[i] = 0;

	outhdr.Omega0 = hdr.Omega0;
	outhdr.OmegaLambda = hdr.OmegaLambda;
	outhdr.HubbleParam = hdr.HubbleParam;
	outhdr.BoxSize = hdr.BoxSize;
	outhdr.flag_sfr = hdr.flag_sfr;
	outhdr.flag_cooling = hdr.flag_cooling;
	outhdr.flag_feedback = hdr.flag_feedback;
	outhdr.flag_age = hdr.flag_age;
	outhdr.flag_metals = hdr.flag_metals;
	outhdr.time = hdr.time;
	outhdr.redshift = hdr.redshift;

	outhdr.npart[1] = (uint32_t) ((((int64_t) numpart_tot) / numfiles) % (1ll << 32));
	outhdr.mass[1] = hdr.mass[1];
	outhdr.num_files = numfiles;
	
	numfiles2 = hdr.num_files;
	
	posbatch = (float *) malloc(3 * (outhdr.npart[1] + numfiles) * sizeof(float));
	velbatch = (float *) malloc(3 * (outhdr.npart[1] + numfiles) * sizeof(float));
#if GADGET_ID_BYTES == 8
	IDbatch = (uint64_t *) malloc((outhdr.npart[1] + numfiles) * sizeof(uint64_t));
#else
	IDbatch = (uint32_t *) malloc((outhdr.npart[1] + numfiles) * sizeof(uint32_t));
#endif

	if (posbatch == NULL || velbatch == NULL || IDbatch == NULL)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to allocate memory for data read!" << endl;
		return -1;
	}

	cout << " redistributing particles..." << endl << endl;
	
	numread = 0;

	for (int nf = 0; nf < numfiles2; nf++)
	{
		cout << " file " << nf << " ..." << endl;

		sprintf(filename, "%s.%d", iname, nf);
			
		infile = fopen(filename, "rb");

		if (infile != NULL)
		{
			fread(&blocksize, sizeof(uint32_t), 1, infile);

			if (blocksize != sizeof(hdr))
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unknown file format " << filename << "!" << endl;
				fclose(infile);
				continue;
			}

			if(fread(&hdr, sizeof(hdr), 1, infile) != 1)
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read header block from " << filename << "!" << endl;
				fclose(infile);
				continue;
			}
				
			fread(&blocksize, sizeof(uint32_t), 1, infile);
			fread(&blocksize, sizeof(uint32_t), 1, infile);

			long blockoffset = 3l * sizeof(float) * (long) hdr.npart[1] + 2l * sizeof(uint32_t);

			backtrack = ftell(infile);

			if (fseek(infile, 2 * blockoffset, SEEK_CUR))
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to fast forward to ID block in " << filename << "!" << endl;
				fclose(infile);
				continue;
			}

			for (int64_t p = 0; p < (int64_t) hdr.npart[1]; p += batch)
			{
				batch = ((int64_t) hdr.npart[1] - p >= (int64_t) outhdr.npart[1] - numread) ? (outhdr.npart[1] - (uint32_t) numread) : (uint32_t) ((int64_t) hdr.npart[1] - p);

				if (
#if GADGET_ID_BYTES == 8
				fread(IDbatch+numread, sizeof(uint64_t), batch, infile)
#else
				fread(IDbatch+numread, sizeof(uint32_t), batch, infile)
#endif
				!= batch)
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read ID batch from " << filename << "!" << endl;
					fclose(infile);
					return -1;
				}

				fastforward = ftell(infile);

				if (fseek(infile, backtrack, SEEK_SET))
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to rewind to positions block in " << filename << "!" << endl;
					fclose(infile);
					return -1;
				}

				if (fread(posbatch+3l*numread, sizeof(float), 3l*batch, infile) != 3l*batch)
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read position data from " << filename << "!" << endl;
					fclose(infile);
					return -1;
				}

				backtrack = ftell(infile);

				if (fseek(infile, blockoffset - 3l * batch * sizeof(float), SEEK_CUR))
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to fast forward to velocities block in " << filename << "!" << endl;
					fclose(infile);
					return -1;
				}

				if (fread(velbatch+3l*numread, sizeof(float), 3l*batch, infile) != 3l*batch)
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read velocity data from " << filename << "!" << endl;
					fclose(infile);
					return -1;
				}

				if (fseek(infile, fastforward, SEEK_SET))
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to fast forward to ID block in " << filename << "!" << endl;
					fclose(infile);
					return -1;
				}

				numread += batch;

				if ((numread == outhdr.npart[1] && numwrite < numfiles-1) || (uint64_t) numread + numpart_write == numpart_tot)
				{
					if (numfiles > 1)
						sprintf(ofilename, "%s.%ld", oname, numwrite);
					else
						sprintf(ofilename, "%s", oname);
						
					outfile = fopen(ofilename, "wb");
	
					if (outfile == NULL)
					{
						fclose(infile);
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to open file " << ofilename << " for output!" << endl;
						return -1;
					}
	
					blocksize = sizeof(outhdr);

					cout << COLORTEXT_CYAN << " writing" << COLORTEXT_RESET << " output file " << ofilename << " ..." << endl;

					fwrite(&blocksize, sizeof(uint32_t), 1, outfile);
					fwrite(&outhdr, sizeof(outhdr), 1, outfile);
					fwrite(&blocksize, sizeof(uint32_t), 1, outfile);

					blocksize = 3l * outhdr.npart[1] * sizeof(float);

					fwrite(&blocksize, sizeof(uint32_t), 1, outfile);

					if (fwrite(posbatch, sizeof(float), 3l * outhdr.npart[1], outfile) != 3l * outhdr.npart[1])
					{
						fclose(infile);
						fclose(outfile);
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to write position block!" << endl;
						return -1;
					}

					fwrite(&blocksize, sizeof(uint32_t), 1, outfile);
					fwrite(&blocksize, sizeof(uint32_t), 1, outfile);

					if (fwrite(velbatch, sizeof(float), 3l * outhdr.npart[1], outfile) != 3l * outhdr.npart[1])
					{
						fclose(infile);
						fclose(outfile);
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to write velocity block!" << endl;
						return -1;
					}

					fwrite(&blocksize, sizeof(uint32_t), 1, outfile);

					blocksize = outhdr.npart[1] * GADGET_ID_BYTES;

					fwrite(&blocksize, sizeof(uint32_t), 1, outfile);

					if (fwrite(IDbatch, GADGET_ID_BYTES, outhdr.npart[1], outfile) != outhdr.npart[1])
					{
						fclose(infile);
						fclose(outfile);
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to write ID block!" << endl;
						return -1;
					}
					
					fwrite(&blocksize, sizeof(uint32_t), 1, outfile);

					fclose(outfile);

					cout << " output file written, contains " << outhdr.npart[1] << " particles." << endl;
	
					numpart_write += outhdr.npart[1];
					numwrite++;

					if (numwrite == numfiles-1)
						outhdr.npart[1] += (uint32_t) (((long) numpart_tot) % numfiles);

					numread = 0;
				}
			}
				
			fclose(infile);
		}
	}

	cout << endl << COLORTEXT_GREEN << " particle redistribution complete." << endl << COLORTEXT_RESET << endl;

	free(IDbatch);
	free(posbatch);
	free(velbatch);

	cout << endl << COLORTEXT_GREEN << " normal completion." << COLORTEXT_RESET << endl;

	return 0;
}

