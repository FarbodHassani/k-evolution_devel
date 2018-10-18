#include <stdlib.h>
#include "LATfield2.hpp"

using namespace LATfield2;

int main(int argc, char **argv)
{
int i,n,m;
string filename;
string outname;
string ids_file;

	for (i=1 ; i < argc ; i++ ){
		if ( argv[i][0] != '-' )
			continue;
		switch(argv[i][1]) {
			case 'n':
				n = atoi(argv[++i]); //size of the dim 1 of the processor grid
				break;
			case 'm':
				m =  atoi(argv[++i]); //size of the dim 2 of the processor grid
				break;
			case 's':
				filename = argv[++i]; //snap file name
				break;
			case 'o':
				outname = argv[++i]; //output file name
				break;
			case 'i':
				ids_file = argv[++i]; //output file name
				break;

		}
	}

FILE *f_in;

f_in = fopen(ids_file.c_str(), "r");
if (f_in == NULL) {
  printf("The file doesn't exist!\n");
  exit(0);}

int rnk;
int npts = n*m;
long int itpp[n*m];

for(i=0;i<n*m;i++)	itpp[i]=0;

parallel.initialize(n,m);

Lattice lat_part(3,npts,0);

part_simple_info particles_global_info;
part_simple_dataType particles_dataType;

part_simple part;

struct fileDsc fd[1];
get_fileDsc_global(filename + ".h5",fd[0]);
get_partInfo(filename+".h5",particles_global_info,particles_dataType);
   
// Here is the iteration over all particles and you can adjust every given condition

// position coordinates
double x,y,z;
// Box length
double lx = fd[0].boxSize[0];
double ly = fd[0].boxSize[1];
double lz = fd[0].boxSize[2];

//Determining the number of data
long int npart_lim=0;
char check_char;

while ((check_char = fgetc(f_in)) != EOF){
  if (check_char == '\n'){
    ++npart_lim;  }}
rewind(f_in);

long int idlst[npart_lim];

i=0;
long int c1;

for(i=0;i<npart_lim;i++){
fscanf(f_in, "%ld", &c1);
idlst[i]=c1;
}

fclose(f_in);

parallel.barrier();

Particles<part_simple,part_simple_info,part_simple_dataType> parts3;
parts3.initialize(particles_global_info,particles_dataType,&lat_part,fd[0].boxSize);

Site xpart(parts3.lattice());   
std::list<part_simple>::iterator it;

Particles<part_simple,part_simple_info,part_simple_dataType> parts4;
parts4.initialize(particles_global_info,particles_dataType,&lat_part,fd[0].boxSize);

parts3.loadHDF5(filename,1);
Site xpart3(parts3.lattice());   

int id;

for(xpart3.first();xpart3.test();xpart3.next())
{
	if(parts3.field()(xpart3).size!=0)
	{
		for(it=parts3.field()(xpart3).parts.begin();it != parts3.field()(xpart3).parts.end();++it)
		{
			id = (*it).ID;
			rnk = parallel.rank();
			for(i=0;i<npart_lim;i++)
			{
				if (id == idlst[i])
				{
						part.ID = id;
						part.pos[0]= (*it).pos[0];
						part.pos[1]= (*it).pos[1];
						part.pos[2]= (*it).pos[2];
						part.vel[0] = (*it).vel[0];
						part.vel[1] = (*it).vel[1];
						part.vel[2] = (*it).vel[2];
						parts4.addParticle_global(part);
//cout<<(*it)<<endl;
				}
			}
		}
	}
}


parts4.saveHDF5(outname,1);
//cout<<particles_global_info.mass<<endl;
//cout<<particles_global_info.relativistic<<endl;
//cout<<fd[0].boxSize[0]<<endl;
}
