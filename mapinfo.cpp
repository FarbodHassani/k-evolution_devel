#include <stdlib.h>
#include <iostream>

// g++ mapinfo.cpp -o mapinfo -std=c++11 -O3

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
	uint32_t header_blocksize;
	uint32_t data_blocksize; /* to get rid of two fread commands per header/data pair */
};

int main(int argc, char **argv)
{
	FILE * infile = NULL;
	healpix_header * hdr;
	uint32_t blocksize[2];
	int count = 0;
#ifdef BUFFEREDIO
	char * buffer;
	size_t bytes_read = 0;
	long offset = 0;

	buffer = (char *) malloc(BUFFEREDIO+4);
#else
	hdr = (healpix_header *) malloc(sizeof(healpix_header));
#endif

	if (argc < 2) return -1;

	infile = fopen(argv[1], "rb");

	if (infile == NULL) return -1;

#ifndef BUFFEREDIO
	setbuf(infile, NULL);

	if (fread(blocksize+1, sizeof(uint32_t), 1, infile) != 1)
	{
		fclose(infile);
		return -1;
	}

	for (count = 0; !feof(infile) && !ferror(infile); count++)
	{
#else
	for (count = 0; true; count++)
	{
		if (bytes_read - offset < sizeof(healpix_header) + sizeof(uint32_t))
		{
			if (feof(infile) || ferror(infile)) break;

			for (long i = 0; i < bytes_read - offset; i++)
				buffer[i] = buffer[offset+i];

			bytes_read += fread(buffer+(bytes_read-offset), sizeof(char), BUFFEREDIO-(bytes_read-offset), infile) - offset;

			if (bytes_read < sizeof(healpix_header) + sizeof(uint32_t)) break;

			cout << " buffering: " << bytes_read << " bytes buffered." << endl;

			offset = 0;
		}

		blocksize[1] = *((uint32_t *) (buffer + offset));
		hdr = (healpix_header *) (buffer + offset + sizeof(uint32_t));

		offset += sizeof(healpix_header) + sizeof(uint32_t);
#endif
		if (blocksize[1] != 256)
		{
			cout << " header block size mismatch! " << blocksize[1] << endl;
			fclose(infile);
			return -1;
		}

#ifndef BUFFEREDIO
		fread(hdr, sizeof(healpix_header), 1, infile);
#endif

		cout << " shell number " << count << ":" << endl;
		cout << "  Nside     = " << hdr->Nside << endl;
		cout << "  Nside_ring= " << hdr->Nside_ring << endl;
		cout << "  Npix      = " << hdr->Npix << endl;
		cout << "  distance  = " << hdr->distance * hdr->boxsize << " Mpc/h" << endl;
		cout << "  direction = (" << hdr->direction[0] << ", " << hdr->direction[1] << ", " << hdr->direction[2] << ")" << endl << endl;

		//fread(&blocksize, sizeof(uint32_t), 1, infile);
		//fread(&blocksize, sizeof(uint32_t), 1, infile);

		if (hdr->data_blocksize != hdr->precision * hdr->Npix)
		{
			cout << " data block size mismatch!" << endl;
			fclose(infile);
			return -1;
		}

#ifndef BUFFEREDIO
		if (fseek(infile, hdr->data_blocksize, SEEK_CUR))
		{
			cout << " unable to skip data block number " << count << "!" << endl;
			fclose(infile);
			return -1;
		}

		fread(blocksize, sizeof(uint32_t), 2, infile);
#else
		blocksize[0] = hdr->data_blocksize;

		if (bytes_read - offset < blocksize[0] + sizeof(uint32_t))
		{
			for (long i = 0; i < bytes_read - offset; i++)
				buffer[i] = buffer[offset+i];

			bytes_read += fread(buffer+(bytes_read-offset), sizeof(char), BUFFEREDIO-(bytes_read-offset), infile) - offset;

			cout << " buffering: " << bytes_read << " bytes buffered." << endl;

			if (bytes_read < blocksize[0] + sizeof(uint32_t))
			{
				cout << " unable to skip data block number " << count << "!" << endl;
				fclose(infile);
				return -1;
			}

			offset = 0;
		}
		
		offset += blocksize[0] + sizeof(uint32_t);
#endif
	}
	fclose(infile);

#ifdef BUFFEREDIO
	free(buffer);
#else
	free(hdr);
#endif

	cout << " file contains " << count << " shells in total." << endl;
}

