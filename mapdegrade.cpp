#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "chealpix.h"

// g++ mapdegrade.cpp -o mapdegrade -std=c++11 -O3 -I/astro/adamek/local/include -L/astro/adamek/local/lib -lchealpix

using namespace std;

struct healpix_header
{
	uint32_t Nside;
	uint32_t Npix;
	uint32_t precision;
	uint32_t Ngrid;
	double direction[3];
	double distance;
	double boxsize;
	char fill[256 - 4 * 4 - 5 * 8]; /* fills to 256 Bytes */
	uint32_t header_blocksize;
	uint32_t data_blocksize; /* to get rid of two fread commands per header/data pair */
};

int main(int argc, char **argv)
{
	FILE * infile = NULL;
	FILE * outfile = NULL;
	healpix_header inhdr;
	healpix_header outhdr;
	uint32_t blocksize[2];
	float * inpix = NULL;
	float * outpix = NULL;
	int * count = NULL;
	int par;
	int Nside_min = 4;
	long Npix_prev = 0;
	uint32_t Nside_prev = 0;
	int64_t index;
	double vec[3];
	
	outhdr.Npix = 0;
	outhdr.precision = 4;
	outhdr.header_blocksize = 256;
	
	if (argc < 3) return -1;

	if (argc > 3)
		Nside_min = atoi(argv[3]);

	infile = fopen(argv[1], "rb");
	
	if (infile == NULL) return -1;

	if (fread(blocksize+1, sizeof(uint32_t), 1, infile) != 1)
	{
		fclose(infile);
		return -1;
	}
	
	outfile = fopen(argv[2], "wb");

	if (outfile == NULL) return -1;

	while (!feof(infile) && !ferror(infile))
	{
		if (blocksize[1] != 256)
		{
			cout << " header block size mismatch! " << blocksize[1] << endl;
			fclose(infile);
			fclose(outfile);
			return -1;
		}
		
		fread(&inhdr, sizeof(healpix_header), 1, infile);
		
		if (inhdr.data_blocksize != 4 * inhdr.Npix)
		{
			cout << " data block size mismatch!" << endl;
			fclose(infile);
			fclose(outfile);
			return -1;
		}
		
		if (fmod(inhdr.distance * inhdr.Ngrid, 2.0) < 0.5 || fmod(inhdr.distance * inhdr.Ngrid, 2.0) > 1.5) // even shell
		{
			outhdr.Nside = (inhdr.Nside / 2 > Nside_min) ? inhdr.Nside / 2 : inhdr.Nside;
			outhdr.Ngrid = inhdr.Ngrid / 2;
			outhdr.boxsize = inhdr.boxsize;
			outhdr.distance = inhdr.distance;
			outhdr.direction[0] = inhdr.direction[0];
			outhdr.direction[1] = inhdr.direction[1];
			outhdr.direction[2] = inhdr.direction[2];
			if ((int64_t) inhdr.Npix >= 12l * (int64_t) inhdr.Nside * (int64_t) inhdr.Nside) outhdr.Npix = 12 * outhdr.Nside * outhdr.Nside;
			else if (outhdr.Nside == inhdr.Nside) outhdr.Npix = inhdr.Npix;
			else
			{
				outhdr.Npix = 0;
				for (int i = 0; (outhdr.Npix - 4 * i) * 4 < inhdr.Npix; i++) outhdr.Npix += 4 * (i + 1);
			}
			outhdr.data_blocksize = 4 * outhdr.Npix;
			
			outpix = (float *) malloc(outhdr.Npix * sizeof(float));
			count = (int *) malloc(outhdr.Npix * sizeof(int));
			
			for (int i = 0; i < outhdr.Npix; i++)
			{
				outpix[i] = 0.;
				count[i] = 0;
			}
			
			for (int64_t i = 0; i < Npix_prev; i++)
			{
				if (inpix[i] <= -1e30) continue;
				pix2vec_ring64(Nside_prev, i, vec);
				vec2pix_ring64(outhdr.Nside, vec, &index);
				if (index < outhdr.Npix)
				{
					if (inhdr.Nside > Nside_prev)
					{
						outpix[index] += 4. * inpix[i];
						count[index] += 4;
					}
					else
					{
						outpix[index] += inpix[i];
						count[index]++;
					}
				}
			}
			
			free(inpix);
			
			inpix = (float *) malloc(inhdr.data_blocksize);
			
			fread(inpix, sizeof(float), inhdr.Npix, infile);
			
			for (int64_t i = 0; i < inhdr.Npix; i++)
			{
				if (inpix[i] <= -1e30) continue;
				pix2vec_ring64(inhdr.Nside, i, vec);
				vec2pix_ring64(outhdr.Nside, vec, &index);
				if (index < outhdr.Npix)
				{
					outpix[index] += 2. * inpix[i];
					count[index] += 2;
				}
			}
			
			free(inpix);
			Npix_prev = 0;
		}
		else // odd shell
		{
			inpix = (float *) malloc(inhdr.data_blocksize);
			
			fread(inpix, sizeof(float), inhdr.Npix, infile);
			
			Npix_prev = inhdr.Npix;
			
			if (outhdr.Npix > 0)
			{
				if (inhdr.Nside > Nside_prev)
					for (int i = 0; i < outhdr.Npix; i++)
					{
						outpix[i] *= 4.;
						count[i] *= 4;
					}
					
				for (int64_t i = 0; i < inhdr.Npix; i++)
				{
					if (inpix[i] <= -1e30) continue;
					pix2vec_ring64(inhdr.Nside, i, vec);
					vec2pix_ring64(outhdr.Nside, vec, &index);
					if (index < outhdr.Npix)
					{
						outpix[index] += inpix[i];
						count[index]++;
					}
				}
				
				for (int i = 0; i < outhdr.Npix; i++)
				{
					if (count[i] > 0) outpix[i] /= count[i];
					else outpix[i] = -1.6375e30;
				}
				
				fwrite(&(outhdr.header_blocksize), sizeof(uint32_t), 1, outfile);
				fwrite(&outhdr, sizeof(healpix_header), 1, outfile);
				fwrite(outpix, sizeof(float), outhdr.Npix, outfile);
				fwrite(&(outhdr.data_blocksize), sizeof(uint32_t), 1, outfile);
				
				free(outpix);
				free(count);
			}
		}
		
		Nside_prev = inhdr.Nside;
		
		fread(blocksize, sizeof(uint32_t), 2, infile);
	}

	if (Npix_prev > 0) free(inpix);
	else if (outhdr.Npix > 0)
	{
		free(outpix);
		free(count);
	}

	fclose(infile);
	fclose(outfile);	
}

