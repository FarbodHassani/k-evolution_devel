#include <stdlib.h>
#include "LATfield2.hpp"

using namespace LATfield2;

int main(int argc, char **argv)
{
int i,n,m;
string filename;
string outname;

m = 2;
n = 2;

	for (i=1 ; i < argc ; i++ ){
		if ( argv[i][0] != '-' )
			continue;
		switch(argv[i][1]) {
			case 's':
				filename = argv[++i]; //snap file name
				break;
			case 'o':
				outname = argv[++i]; //snap file name
				break;
		}
	}

parallel.initialize(n,m);

Lattice lat_part(3,16,0);

part_simple_info particles_global_info;
part_simple_dataType particles_dataType;

struct fileDsc fd[1];
get_fileDsc_global(filename + ".h5",fd[0]);
get_partInfo(filename+".h5",particles_global_info,particles_dataType);

Particles<part_simple,part_simple_info,part_simple_dataType> parts;
parts.initialize(particles_global_info,particles_dataType,&lat_part,fd[0].boxSize);

parts.loadHDF5(filename,1);

int num=0;

Site xpart(parts.lattice());   
std::list<part_simple>::iterator it;

for(xpart.first();xpart.test();xpart.next())
{
	if(parts.field()(xpart).size!=0)
	{
		for(it=parts.field()(xpart).parts.begin();it != parts.field()(xpart).parts.end();++it)
		{
 num++;

		}
	}
}

//cout<<num<<endl;

parts.saveHDF5(outname,1);

}
