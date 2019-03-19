#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "chealpix.h"

// g++ map2fits.cpp -o map2fits -std=c++11 -O3 -I/astro/adamek/local/include -L/astro/adamek/local/lib -lchealpix -lcfitsio

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
	float * pix;
	uint32_t blocksize;
	double * dpix;
	float * fpix;
	float * rpix;
	char coordsys = 'G';
	int skip = 0;
	int Nside = -1;
	double vec[3];
	int64_t j, q;
	int ring;
	int pixbatch_delim[3];
	int pixbatch_size[3] = {0, 0, 0};

	if (argc < 3) return -1;
	else if (argc > 3) skip = atoi(argv[3]);

	if (argc > 4) Nside = atoi(argv[4]);

	infile = fopen(argv[1], "rb");

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

		fread(&hdr, sizeof(hdr), 1, infile);
		fread(&blocksize, sizeof(uint32_t), 1, infile);
		fread(&blocksize, sizeof(uint32_t), 1, infile);

		if (blocksize != hdr.precision * hdr.Npix)
		{
			cout << " data block size mismatch!" << endl;
			fclose(infile);
			return -1;
		}

		if (i == skip) break;

		if (fseek(infile, blocksize, SEEK_CUR))
		{
			cout << " unable to skip shell number " << i << "!" << endl;
			fclose(infile);
			return -1;
		}

		fread(&blocksize, sizeof(uint32_t), 1, infile);
	}

	if (feof(infile))
	{
		if (skip) cout << " unexpected EOF occurred when skipping " << skip << " shells!" << endl;
		else cout << " unexpected EOF occurred!" << endl;
		fclose(infile);
		return -1;
	}
	else if (ferror(infile))
	{
		cout << " an error occurred!" << endl;
		fclose(infile);
		return -1;
	}

	if (Nside > hdr.Nside || Nside < 0) Nside = hdr.Nside;
	
	if (hdr.Nside_ring > 0 && hdr.Nside_ring < hdr.Nside)
	{
		cout << " pixels will be reordered" << endl;
		
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

	pix = (float *) malloc(Nside * Nside * 12 * sizeof(float));

	if (hdr.precision == 4)
	{
		if (Nside == hdr.Nside)
		{
			if (pixbatch_size[0] > 1)
			{
				rpix = (float *) malloc (hdr.Npix * sizeof(float));
				fread(rpix, sizeof(float), hdr.Npix, infile);
				for (int p = 0; p < pixbatch_delim[0]; p++)
				{
					for (int i = 0; i < pixbatch_size[0]; i++)
					{
						ring2nest64(hdr.Nside_ring, p, &j);
						j = j*pixbatch_size[0] + i;
						nest2ring64(hdr.Nside, j, &j);
						pix[j] = rpix[pixbatch_size[0]*p + i];
					}
				}
				for (int p = pixbatch_delim[0]; p < pixbatch_delim[1]; p++)
				{
					q = 0;
					for (int i = 0; i < pixbatch_size[0]; i++)
					{
						ring2nest64(hdr.Nside_ring, p, &j);
						j = j*pixbatch_size[0] + i;
						nest2ring64(hdr.Nside, j, &j);
						if (j < hdr.Npix)
						{
							pix[j] = rpix[pixbatch_size[0]*pixbatch_delim[0] + pixbatch_size[1]*(p-pixbatch_delim[0]) + q];
							q++;
						}
					}
				}
				for (int p = pixbatch_delim[1]; p < pixbatch_delim[2]; p++)
				{
					q = 0;
					for (int i = 0; i < pixbatch_size[0]; i++)
					{
						ring2nest64(hdr.Nside_ring, p, &j);
						j = j*pixbatch_size[0] + i;
						nest2ring64(hdr.Nside, j, &j);
						if (j < hdr.Npix)
						{
							pix[j] = rpix[pixbatch_size[0]*pixbatch_delim[0] + pixbatch_size[1]*(pixbatch_delim[1]-pixbatch_delim[0]) + pixbatch_size[2]*(p-pixbatch_delim[1]) + q];
							q++;
						}
					}
				}
				free(rpix);
			}
			else
				fread(pix, sizeof(float), hdr.Npix, infile);
		}
		else
		{
			cout << " downgrading from Nside = " << hdr.Nside << " to Nside = " << Nside << endl;
			fpix = (float *) malloc (hdr.Npix * sizeof(float));
			if (pixbatch_size[0] > 1)
			{
				rpix = (float *) malloc (hdr.Npix * sizeof(float));
				fread(rpix, sizeof(float), hdr.Npix, infile);
				for (int p = 0; p < pixbatch_delim[0]; p++)
				{
					for (int i = 0; i < pixbatch_size[0]; i++)
					{
						ring2nest64(hdr.Nside_ring, p, &j);
						j = j*pixbatch_size[0] + i;
						nest2ring64(hdr.Nside, j, &j);
						fpix[j] = rpix[pixbatch_size[0]*p + i];
					}
				}
				for (int p = pixbatch_delim[0]; p < pixbatch_delim[1]; p++)
				{
					q = 0;
					for (int i = 0; i < pixbatch_size[0]; i++)
					{
						ring2nest64(hdr.Nside_ring, p, &j);
						j = j*pixbatch_size[0] + i;
						nest2ring64(hdr.Nside, j, &j);
						if (j < hdr.Npix)
						{
							fpix[j] = rpix[pixbatch_size[0]*pixbatch_delim[0] + pixbatch_size[1]*(p-pixbatch_delim[0]) + q];
							q++;
						}
					}
				}
				for (int p = pixbatch_delim[1]; p < pixbatch_delim[2]; p++)
				{
					q = 0;
					for (int i = 0; i < pixbatch_size[0]; i++)
					{
						ring2nest64(hdr.Nside_ring, p, &j);
						j = j*pixbatch_size[0] + i;
						nest2ring64(hdr.Nside, j, &j);
						if (j < hdr.Npix)
						{
							fpix[j] = rpix[pixbatch_size[0]*pixbatch_delim[0] + pixbatch_size[1]*(pixbatch_delim[1]-pixbatch_delim[0]) + pixbatch_size[2]*(p-pixbatch_delim[1]) + q];
							q++;
						}
					}
				}
				free(rpix);
			}
			else
				fread(fpix, sizeof(float), hdr.Npix, infile);
			cout << " a few pixel values: ";
			for (int i = 0; i < 12; i++) cout << fpix[i] << " ";
			cout << endl;
			for (int i = 0; i < 12 * Nside * Nside; i++) pix[i] = 0.;
			for (int i = 0; i < hdr.Npix; i++)
			{
				pix2vec_ring64(hdr.Nside, (long) i, vec);
				vec2pix_ring64(Nside, vec, &j);
				pix[j] += fpix[i] / (hdr.Nside / Nside) / (hdr.Nside / Nside);
			}
			free(fpix);
			for (int i = 0; i < hdr.Npix / (hdr.Nside / Nside) / (hdr.Nside / Nside); i++)
			{
				if (pix[i] <= -1e30 / (hdr.Nside / Nside) / (hdr.Nside / Nside)) pix[i] = -1.6375e30;
			}
		}
		fclose(infile);
	}
	else if (hdr.precision == 8)
	{
		dpix = (double *) malloc (hdr.Npix * sizeof(double));
		fread(dpix, sizeof(double), hdr.Npix, infile);
		fclose(infile);

		if (pixbatch_size[0] > 1)
		{
			for (int p = 0; p < pixbatch_delim[0]; p++)
			{
				for (int i = 0; i < pixbatch_size[0]; i++)
				{
					ring2nest64(hdr.Nside_ring, p, &j);
					j = j*pixbatch_size[0] + i;
					nest2ring64(hdr.Nside, j, &j);
					pix[j] = dpix[pixbatch_size[0]*p + i];
				}
			}
			for (int p = pixbatch_delim[0]; p < pixbatch_delim[1]; p++)
			{
				q = 0;
				for (int i = 0; i < pixbatch_size[0]; i++)
				{
					ring2nest64(hdr.Nside_ring, p, &j);
					j = j*pixbatch_size[0] + i;
					nest2ring64(hdr.Nside, j, &j);
					if (j < hdr.Npix)
					{
						pix[j] = dpix[pixbatch_size[0]*pixbatch_delim[0] + pixbatch_size[1]*(p-pixbatch_delim[0]) + q];
						q++;
					}
				}
			}
			for (int p = pixbatch_delim[1]; p < pixbatch_delim[2]; p++)
			{
				q = 0;
				for (int i = 0; i < pixbatch_size[0]; i++)
				{
					ring2nest64(hdr.Nside_ring, p, &j);
					j = j*pixbatch_size[0] + i;
					nest2ring64(hdr.Nside, j, &j);
					if (j < hdr.Npix)
					{
						pix[j] = dpix[pixbatch_size[0]*pixbatch_delim[0] + pixbatch_size[1]*(pixbatch_delim[1]-pixbatch_delim[0]) + pixbatch_size[2]*(p-pixbatch_delim[1]) + q];
						q++;
					}
				}
			}
		}
		else
		{
			for (int i = 0; i < hdr.Npix; i++) pix[i] = dpix[i];
		}

		free(dpix);
	}
	else
	{
		cout << " precision " << hdr.precision << " bytes not supported!" << endl;
		fclose(infile);
		return -1;
	}

	for (int i = hdr.Npix / (hdr.Nside / Nside) / (hdr.Nside / Nside); i < 12 * Nside * Nside; i++)
		pix[i] = -1.6375e30;

	write_healpix_map(pix, (long) Nside, argv[2], 0, &coordsys);

	free(pix);

	return 0;
}

