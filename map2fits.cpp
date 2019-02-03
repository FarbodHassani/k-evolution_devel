#include <stdlib.h>
#include <iostream>
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
	char fill[256 - 4 * 4 - 5 * 8]; /* fills to 256 Bytes */
};

int main(int argc, char **argv)
{
	FILE * infile = NULL;
	healpix_header hdr;
	float * pix;
	uint32_t blocksize;
	double * dpix;
	float * fpix;
	char coordsys = 'G';
	int skip = 0;
	int Nside = -1;
	double vec[3];
	long j;

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

	pix = (float *) malloc(Nside * Nside * 12 * sizeof(float));

	if (hdr.precision == 4)
	{
		if (Nside == hdr.Nside)
			fread(pix, sizeof(float), hdr.Npix, infile);
		else
		{
			cout << " downgrading from Nside = " << hdr.Nside << " to Nside = " << Nside << endl;
			fpix = (float *) malloc (hdr.Npix * sizeof(float));
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

		for (int i = 0; i < hdr.Npix; i++) pix[i] = dpix[i];

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

