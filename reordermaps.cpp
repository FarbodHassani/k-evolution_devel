#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <cstring>
#include "chealpix.h"

// g++ reordermaps.cpp -o reordermaps -std=c++11 -O3 -I/astro/adamek/local/include -L/astro/adamek/local/lib -lchealpix -lcfitsio

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
	uint32_t Nside_ring;
	char fill[256 - 5 * 4 - 5 * 8]; /* fills to 256 Bytes */
};

int main(int argc, char **argv)
{
	FILE * infile = NULL;
	healpix_header hdr;
	uint32_t blocksize;
	char * pix;
	char * temp;
	long backtrack, fastforward;
	int64_t j, q;
	int ring;
	int pixbatch_delim[3];
	int pixbatch_size[3] = {0, 0, 0};

	if (argc < 2) return -1;

	infile = fopen(argv[1], "r+b");

	if (infile == NULL) return -1;

	for (int i = 0; !feof(infile) && !ferror(infile); i++)
	{
		if (fread(&blocksize, sizeof(uint32_t), 1, infile) != 1) break;

		if (blocksize != 256)
		{
			cout << " header block size mismatch!" << endl;
			fclose(infile);
			return -1;
		}

		backtrack = ftell(infile);
		fread(&hdr, sizeof(hdr), 1, infile);
		fread(&blocksize, sizeof(uint32_t), 1, infile);
		fread(&blocksize, sizeof(uint32_t), 1, infile);

		if (blocksize != hdr.precision * hdr.Npix)
		{
			cout << " data block size mismatch!" << endl;
			fclose(infile);
			return -1;
		}

		if (hdr.Nside_ring > 0 && hdr.Nside_ring < hdr.Nside)
		{
			cout << " pixels will be reordered for shell " << i << endl;
		
			if ((long) hdr.Npix <= 2 * (long) hdr.Nside * (hdr.Nside + 1))
				ring = (int) floor((sqrt(2. * hdr.Npix + 1.01) - 1.) / 2.);
			else if ((long) hdr.Npix <= 2 * (long) hdr.Nside * (hdr.Nside + 1) + 4 * (2 * hdr.Nside - 1) * (long) hdr.Nside)
				ring = ((hdr.Npix - 2 * hdr.Nside * (hdr.Nside + 1)) / 4 / hdr.Nside) + hdr.Nside;
			else if ((long) hdr.Npix < 12 * (long) hdr.Nside * hdr.Nside)
			{
				ring = 12 * (long) hdr.Nside * hdr.Nside - (long) hdr.Npix;
				ring = (int) floor((sqrt(2. * ring + 1.01) - 1.) / 2.);
				ring = 4 * hdr.Nside - 1 - ring;
			}
			else
				ring = 4 * hdr.Nside - 1;
			
			pixbatch_size[0] = (hdr.Nside / hdr.Nside_ring);
				
			pixbatch_delim[1] = ring / pixbatch_size[0];
			pixbatch_delim[0] = (pixbatch_delim[1] > 0) ? pixbatch_delim[1]-1 : 0;
			pixbatch_delim[2] = pixbatch_delim[1]+1;
			pixbatch_size[1] = (pixbatch_size[0] * (pixbatch_size[0]+1) + (2*pixbatch_size[0] - 1 - ring%pixbatch_size[0]) * (ring%pixbatch_size[0])) / 2;
			//pixbatch_size[2] = ((pixbatch_size[0] - ring%pixbatch_size[0] - 1) * (pixbatch_size[0] - ring%pixbatch_size[0])) / 2;
			pixbatch_size[2] = ((ring%pixbatch_size[0] + 1) * (ring%pixbatch_size[0])) / 2;
			pixbatch_size[0] *= pixbatch_size[0];
			for (int p = 0; p < 3; p++)
			{
				if (pixbatch_delim[p] <= hdr.Nside_ring)
					pixbatch_delim[p] = 2 * pixbatch_delim[p] * (pixbatch_delim[p]+1);
				else if (pixbatch_delim[p] <= 3 * hdr.Nside_ring)
					pixbatch_delim[p] = 2 * hdr.Nside_ring * (hdr.Nside_ring+1) + (pixbatch_delim[p]-hdr.Nside_ring) * 4 * hdr.Nside_ring;
				else if (pixbatch_delim[p] < 4 * hdr.Nside_ring)
					pixbatch_delim[p] = 12 * hdr.Nside_ring * hdr.Nside_ring - 2 * (4 * hdr.Nside_ring - 1 - pixbatch_delim[p]) * (4 * hdr.Nside_ring - pixbatch_delim[p]);
				else
					pixbatch_delim[p] = 12 * hdr.Nside_ring * hdr.Nside_ring;
			}
		
			cout << " full patches (" << pixbatch_size[0] << " pixels): " << pixbatch_delim[0] << ", cut patches (" << pixbatch_size[1] << " pixels): " << pixbatch_delim[1]-pixbatch_delim[0] << ", cut patches (" << pixbatch_size[2] << " pixels): " << pixbatch_delim[2]-pixbatch_delim[1] << endl;
		}
		
		pix = (char *) malloc(hdr.precision * hdr.Npix);
		
		if (pixbatch_size[0] > 1)
		{
			temp = (char *) malloc(hdr.precision * hdr.Npix);
			fread(temp, hdr.precision, hdr.Npix, infile);
			for (int p = 0; p < pixbatch_delim[0]; p++)
			{
				for (int k = 0; k < pixbatch_size[0]; k++)
				{
					ring2nest64(hdr.Nside_ring, p, &j);
					j = j*pixbatch_size[0] + k;
					nest2ring64(hdr.Nside, j, &j);
					memcpy(pix+(j*hdr.precision), temp+((pixbatch_size[0]*p + k)*hdr.precision), hdr.precision);
				}
			}
			for (int p = pixbatch_delim[0]; p < pixbatch_delim[1]; p++)
			{
				q = 0;
				for (int k = 0; k < pixbatch_size[0]; k++)
				{
					ring2nest64(hdr.Nside_ring, p, &j);
					j = j*pixbatch_size[0] + k;
					nest2ring64(hdr.Nside, j, &j);
					if (j < hdr.Npix)
					{
						memcpy(pix+(j*hdr.precision), temp+((pixbatch_size[0]*pixbatch_delim[0] + pixbatch_size[1]*(p-pixbatch_delim[0]) + q)*hdr.precision), hdr.precision);
						q++;
					}
				}
			}
			for (int p = pixbatch_delim[1]; p < pixbatch_delim[2]; p++)
			{
				q = 0;
				for (int k = 0; k < pixbatch_size[0]; k++)
				{
					ring2nest64(hdr.Nside_ring, p, &j);
					j = j*pixbatch_size[0] + k;
					nest2ring64(hdr.Nside, j, &j);
					if (j < hdr.Npix)
					{
						memcpy(pix+(j*hdr.precision), temp+((pixbatch_size[0]*pixbatch_delim[0] + pixbatch_size[1]*(pixbatch_delim[1]-pixbatch_delim[0]) + pixbatch_size[2]*(p-pixbatch_delim[1]) + q)*hdr.precision), hdr.precision);
						q++;
					}
				}
			}
			free(temp);
		}
		else
			fread(pix, hdr.precision, hdr.Npix, infile);

		fread(&blocksize, sizeof(uint32_t), 1, infile);
		
		fastforward = ftell(infile);
		
		hdr.Nside_ring = 0;
		
		if(fseek(infile, backtrack, SEEK_SET) != 0)
		{
			cout << " unable to backtrack!" << endl;
			free(pix);
			break;
		}
		
		fwrite(&hdr, sizeof(hdr), 1, infile);
		
		fread(&blocksize, sizeof(uint32_t), 1, infile);
		fread(&blocksize, sizeof(uint32_t), 1, infile);
		
		fwrite(pix, hdr.precision, hdr.Npix, infile);
		
		free(pix);
		
		if(fseek(infile, fastforward, SEEK_SET) != 0)
		{
			cout << " unable to fastforward!" << endl;
			break;
		}
	}
	
	fclose(infile);

	if (ferror(infile))
	{
		cout << " an error occurred!" << endl;
		return -1;
	}

	return 0;
}

