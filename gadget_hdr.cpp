#include <stdlib.h>
#include <stdio.h>
#include <cstdint>

using namespace std;

struct io_header_1
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
} header1;


int main(int argc, char **argv)
{
  FILE * pclfile;
  uint32_t blocksize1, blocksize2, i;
  uint64_t num64bit;
  
  if (argc < 2)
  {
  	printf("too few parameters!\n");
  	return 0;
  }
  
  pclfile = fopen(argv[1], "r+");
  
  if (pclfile == NULL)
  {
    printf("file not found!\n");
    return 0;
  }
  else
  {
    printf("reading %s ...\n", argv[1]);
    fread(&blocksize1, sizeof(uint32_t), 1, pclfile);
    fread(&header1, sizeof(header1), 1, pclfile);
    fread(&blocksize2, sizeof(uint32_t), 1, pclfile);
    
    if (blocksize1 == blocksize2)
    {
      printf("header block (%d bytes) read successfully.\n", blocksize1);
      
      printf(" -- HEADER --\n");
      printf(" npart = %u\t%u\t%u\t%u\t%u\t%u\n", header1.npart[0], header1.npart[1], header1.npart[2], header1.npart[3], header1.npart[4], header1.npart[5]);
      printf(" mass  = %f\t%f\t%f\t%f\t%f\t%f\n", header1.mass[0], header1.mass[1], header1.mass[2], header1.mass[3], header1.mass[4], header1.mass[5]);
      printf(" time = %f\n", header1.time);
      printf(" redshift = %f\n", header1.redshift);
      printf(" flag_sfr = %d\n", header1.flag_sfr);
      printf(" flag_feedback = %d\n", header1.flag_feedback);
      printf(" npartTotal = %u\t%u\t%u\t%u\t%u\t%u\n", header1.npartTotal[0], header1.npartTotal[1], header1.npartTotal[2], header1.npartTotal[3], header1.npartTotal[4], header1.npartTotal[5]);
      printf(" flag_cooling = %d\n", header1.flag_cooling);
      printf(" num_files = %u\n", header1.num_files);
      printf(" BoxSize = %f\n", header1.BoxSize);
      printf(" Omega0 = %f\n", header1.Omega0);
      printf(" OmegaLambda = %f\n", header1.OmegaLambda);
      printf(" HubbleParam = %f\n", header1.HubbleParam);
      printf(" flag_age = %d\n", header1.flag_age);
      printf(" flag_metals = %d\n", header1.flag_metals);
      printf(" npartTotalHW = %u\t%u\t%u\t%u\t%u\t%u\n", header1.npartTotalHW[0], header1.npartTotalHW[1], header1.npartTotalHW[2], header1.npartTotalHW[3], header1.npartTotalHW[4], header1.npartTotalHW[5]);
      printf(" fill = ");
      for (i = 0; i < 256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8 - 2 * 4 - 6 * 4; i++) printf("%02X", header1.fill[i]);
      printf("\n");
    }
    else
    {
      printf("error reading header block; block sizes where %d and %d, respectively.\n", blocksize1, blocksize2);
      fclose(pclfile);
    }
    
    if (argc > 2 && argv[2][0] == 'r')
    {
    	header1.OmegaLambda = 1. - header1.Omega0;
    	printf("setting OmegaLambda = 1 - Omega0\n");
    }
	else if (argc > 2 && argv[2][0] >= '1' && argv[2][0] < '9')
	{
		for (i = 0; i < 6; i++)
			header1.mass[i] *= (double) header1.num_files / (double) (argv[2][0]-'0');

		for (i = 0; i < 6; i++)
		{
			num64bit = (uint64_t) (argv[2][0]-'0') * (uint64_t) header1.npart[i];
			header1.npartTotal[i] = (uint32_t) (num64bit % (uint64_t) (1l << 32));
			header1.npartTotalHW[i] = (uint32_t) (num64bit >> 32);
		}

		header1.num_files = (uint32_t) (argv[2][0]-'0');
		printf("setting number of files = %d\n", header1.num_files);
	}
    else
    {
    	printf("don't know what to do, exiting.\n");
    	fclose(pclfile);
    	return 0;
    }
    rewind(pclfile);
    fread(&blocksize1, sizeof(uint32_t), 1, pclfile);
    fwrite(&header1, sizeof(header1), 1, pclfile);
    fclose(pclfile);
  }

  return 0;
}

