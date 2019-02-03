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
  int32_t num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];       /* fills to 256 Bytes */
} header1;


int main(int argc, char **argv)
{
  FILE * infile;
  FILE * outfile;
  uint32_t blocksize1, blocksize2;
  float * fchunk;
  int64_t * lchunk;
  char filename[256];
  int64_t pos, vel, id, n, count, batch, rem[6];
  int i, j, numfile;
  
  if (argc < 3)
  {
  	printf("too few parameters!\n");
  	return 0;
  }
  
  numfile = atoi(argv[2]);
  if (numfile < 1)
  {
  	printf("number of files not recognized!\n");
  	return 0;
  }

  infile = fopen(argv[1], "r");
  
  if (infile == NULL)
  {
    printf("input file not found!\n");
    return 0;
  }
  
  printf("reading %s ...\n", argv[1]);
  fread(&blocksize1, sizeof(uint32_t), 1, infile);
  fread(&header1, sizeof(header1), 1, infile);
  fread(&blocksize2, sizeof(uint32_t), 1, infile);
    
  if (blocksize1 == blocksize2)
  {
    printf("header block (%u bytes) read successfully.\n", blocksize1);
      
    printf(" -- HEADER --\n");
    printf(" npart = %u\t%u\t%u\t%u\t%u\t%u\n", header1.npart[0], header1.npart[1], header1.npart[2], header1.npart[3], header1.npart[4], header1.npart[5]);
    printf(" mass  = %f\t%f\t%f\t%f\t%f\t%f\n", header1.mass[0], header1.mass[1], header1.mass[2], header1.mass[3], header1.mass[4], header1.mass[5]);
    printf(" time = %f\n", header1.time);
    printf(" redshift = %f\n", header1.redshift);
    printf(" flag_sfr = %d\n", header1.flag_sfr);
    printf(" flag_feedback = %d\n", header1.flag_feedback);
    printf(" npartTotal = %u\t%u\t%u\t%u\t%u\t%u\n", header1.npartTotal[0], header1.npartTotal[1], header1.npartTotal[2], header1.npartTotal[3], header1.npartTotal[4], header1.npartTotal[5]);
    printf(" flag_cooling = %d\n", header1.flag_cooling);
    printf(" num_files = %d\n", header1.num_files);
    printf(" BoxSize = %f\n", header1.BoxSize);
    printf(" Omega0 = %f\n", header1.Omega0);
    printf(" OmegaLambda = %f\n", header1.OmegaLambda);
    printf(" HubbleParam = %f\n", header1.HubbleParam);
    printf(" fill = ");
    for (i = 0; i < 256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8; i++) printf("%2X", header1.fill[i]);
    printf("\n");
  }
  else
  {
    printf("error reading header block; block sizes where %u and %u, respectively.\n", blocksize1, blocksize2);
    fclose(infile);
    return 0;
  }

  for (i = 0; i < 6; i++)
  {
	rem[i] = header1.npart[i] % numfile;
  	header1.npart[i] /= numfile;
  }
  header1.num_files *= numfile;
  
  fread(&blocksize2, sizeof(uint32_t), 1, infile);
  pos = ftell(infile);
  if (fseek(infile, (long) header1.npartTotal[1] * (long) sizeof(float) * 3l + (long) sizeof(uint32_t) * 2l, SEEK_CUR) != 0)
  {
  	printf("could not fast-forward to velocities block!\n");
  	fclose(infile);
  	return 0;
  }
  vel = ftell(infile);
  if (fseek(infile, (long) header1.npartTotal[1] * (long) sizeof(float) * 3l + (long) sizeof(uint32_t) * 2l, SEEK_CUR) != 0)
  {
  	printf("could not fast-forward to ID block!\n");
  	fclose(infile);
  	return 0;
  }
  id = ftell(infile);

  fchunk = (float *) malloc(402653184 * sizeof(float));
  lchunk = (int64_t *) malloc(134217728 * sizeof(int64_t));
  
  for (j = 0; j < numfile; j++)
  {  
    sprintf(filename, "%s.%d", argv[1], j);
    
    outfile = fopen(filename, "w");
    if (outfile == NULL)
  	{
      printf("could not create output file!\n");
      fclose(infile);
      free(fchunk);
      free(lchunk);
      return 0;
  	}

	if (j == numfile-1)
	  header1.npart[1] += rem[1];

 	blocksize1 = (uint32_t) sizeof(header1);
 	fwrite(&blocksize1, sizeof(uint32_t), 1, outfile);
 	fwrite(&header1, sizeof(header1), 1, outfile);
 	fwrite(&blocksize1, sizeof(uint32_t), 1, outfile);
 	blocksize1 = (uint32_t) header1.npart[1] * sizeof(float) * 3;
 	fwrite(&blocksize1, sizeof(uint32_t), 1, outfile);
 	if (fseek(infile, pos, SEEK_SET) != 0)
  	{
  	  printf("could not backtrack to positions block!\n");
  	  fclose(infile);
  	  fclose(outfile);
      free(fchunk);
      free(lchunk);
  	  return 0;
    }
 	
    printf("copying positions (%u bytes)...\n", blocksize1);
    for (i = 0; i < header1.npart[1];)
    {
      	if (header1.npart[1] - i > 134217728) batch = 402653184;
   		else batch = (header1.npart[1] - i) * 3;
      	fread(fchunk, sizeof(float), batch, infile);
        fwrite(fchunk, sizeof(float), batch, outfile);
        i += batch / 3;
    }
    pos = ftell(infile);
    fwrite(&blocksize1, sizeof(uint32_t), 1, outfile);
    fwrite(&blocksize1, sizeof(uint32_t), 1, outfile);
    if (fseek(infile, vel, SEEK_SET) != 0)
  	{
  	  printf("could not fast-forward to velocities block!\n");
  	  fclose(infile);
  	  fclose(outfile);
      free(fchunk);
      free(lchunk);
  	  return 0;
    }
 	
    printf("copying velocities (%u bytes)...\n", blocksize1);
    for (i = 0; i < header1.npart[1];)
    {
      	if (header1.npart[1] - i > 134217728) batch = 402653184;
   		else batch = (header1.npart[1] - i) * 3;
      	fread(fchunk, sizeof(float), batch, infile);
        fwrite(fchunk, sizeof(float), batch, outfile);
        i += batch / 3;
    }
    vel = ftell(infile);
    fwrite(&blocksize1, sizeof(uint32_t), 1, outfile);
    blocksize1 = (uint32_t) header1.npart[1] * sizeof(int64_t);
 	fwrite(&blocksize1, sizeof(uint32_t), 1, outfile);
 	if (fseek(infile, id, SEEK_SET) != 0)
  	{
  	  printf("could not fast-forward to ID block!\n");
  	  fclose(infile);
  	  fclose(outfile);
      free(fchunk);
      free(lchunk);
  	  return 0;
    }
 	
    printf("copying IDs (%u bytes)...\n", blocksize1);
    for (i = 0; i < header1.npart[1];)
    {
      	if (header1.npart[1] - i > 134217728) batch = 134217728;
   		else batch = (header1.npart[1] - i);
      	fread(lchunk, sizeof(int64_t), batch, infile);
        fwrite(lchunk, sizeof(int64_t), batch, outfile);
        i += batch;
    }
    id = ftell(infile);
    fwrite(&blocksize1, sizeof(uint32_t), 1, outfile);
    fclose(outfile);
    printf("%s written successfully.\n", filename);
  }
  
  fclose(infile);

  free(fchunk);
  free(lchunk);

  return 0;
}

