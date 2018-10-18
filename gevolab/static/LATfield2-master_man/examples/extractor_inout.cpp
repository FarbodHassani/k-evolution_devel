#include <stdlib.h>
#include "LATfield2.hpp"

using namespace LATfield2;

int main(int argc, char **argv)
{
int n,m;

string filename1 = "./files/lcdm_snap000_part";
string filename2 = "./files/lcdm_snap003_part";
string outname = "./files/out";

n = 2;
m = 2;
int i,rnk;
int npts = n*m;
long int itpp[n*m];

for(i=0;i<n*m;i++)	itpp[i]=0;

parallel.initialize(n,m);

Lattice lat_part(3,npts,0);

part_simple_info particles_global_info;
part_simple_dataType particles_dataType;

part_simple part;

struct fileDsc fd[1];
get_fileDsc_global(filename1 + ".h5",fd[0]);
get_partInfo(filename1+".h5",particles_global_info,particles_dataType);

Particles<part_simple,part_simple_info,part_simple_dataType> parts;
parts.initialize(particles_global_info,particles_dataType,&lat_part,fd[0].boxSize);

Particles<part_simple,part_simple_info,part_simple_dataType> parts2;
parts2.initialize(particles_global_info,particles_dataType,&lat_part,fd[0].boxSize);

parts.loadHDF5(filename1,1);

Site xpart(parts.lattice());   
std::list<part_simple>::iterator it;

long int npart_tot=0;
for(xpart.first();xpart.test();xpart.next())
{
	if(parts.field()(xpart).size!=0)
	{
		for(it=parts.field()(xpart).parts.begin();it != parts.field()(xpart).parts.end();++it)
		{
			npart_tot++;
		}
	}
}
parallel.barrier();
parallel.sum(npart_tot);

COUT<<"Total number of particles: "<<npart_tot<<endl;

long int idlist[npart_tot][n*m];
   
// Here is the iteration over all particles and you can adjust every given condition

// position coordinates
double x,y,z;
long int npart_lim=0 ;
// Box length
double lx = fd[0].boxSize[0];
double ly = fd[0].boxSize[1];
double lz = fd[0].boxSize[2];

for(xpart.first();xpart.test();xpart.next())
{
	if(parts.field()(xpart).size!=0)
	{
		for(it=parts.field()(xpart).parts.begin();it != parts.field()(xpart).parts.end();++it)
		{
			x = (*it).pos[0]-0.5;
			y = (*it).pos[1]-0.5;
			z = (*it).pos[2]-0.5;

			if ( sqrt(x*x + y*y + z*z) < 0.2 ){
				part.ID = (*it).ID;
//				part.pos[0]= (*it).pos[0];
//				part.pos[1]= (*it).pos[1];
//				part.pos[2]= (*it).pos[2];
//				part.vel[0] = (*it).vel[0];
//				part.vel[1] = (*it).vel[1];
//				part.vel[2] = (*it).vel[2];
//				parts2.addParticle_global(part);
				npart_lim++;

				rnk = parallel.rank();
				idlist[itpp[rnk]][rnk] = (*it).ID;
				itpp[rnk]++;

			}
		}
	}
}

parallel.barrier();
parallel.sum(npart_lim);
COUT<<"Number of limmited particles: "<<npart_lim<<endl;
parallel.barrier();


long int idlst[npart_lim];
parallel.barrier();

int j = 0;

for (rnk=0;rnk<n*m;rnk++){
for(i=0;i<itpp[rnk];i++){
idlst[j] = idlist[i][rnk];
j++;
}
}

parallel.barrier();

Particles<part_simple,part_simple_info,part_simple_dataType> parts3;
parts3.initialize(particles_global_info,particles_dataType,&lat_part,fd[0].boxSize);
particles_global_info.mass *= ((double)(npart_tot)/npart_lim);

Particles<part_simple,part_simple_info,part_simple_dataType> parts4;
parts4.initialize(particles_global_info,particles_dataType,&lat_part,fd[0].boxSize);

parts3.loadHDF5(filename2,1);
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
