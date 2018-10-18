#include <stdlib.h>
#include "LATfield2.hpp"

using namespace LATfield2;

int main(int argc, char **argv)
{


    int n,m;
    int io_groupe_size,io_size;
    string str_filename;
    int npts = 4;
    int numparts = 4;
    Real  latresolution =0.1;


    n = 2;
    m = 2;


    parallel.initialize(n,m);

        int dim=3;
        int halo=1;
        int khalo=1;


        Lattice lat_part(dim,npts,0);
        Lattice lat(dim,npts,halo);

        Field<Real> phi(lat,1);
        Field<Real> B(lat,3);
        Field<Real> Tij(lat,3,3,LATfield2::symmetric);

        Real boxSize[3];
//        for(int i=0;i<3;i++)boxSize[i] = latresolution * lat_part.size(i);

        for(int i=0;i<3;i++)boxSize[i] = 1;

        double timerRef;
        double timerWrite,timerLoad,timerWriteServer;

        double timerProjScalar,timerProjVector,timerProjTensor;
        double timerCommScalar,timerCommVector,timerCommTensor;

        double timerMove,timerVel;

        Site x(lat);

        //cout<<"start part init done"<<endl;

        part_simple_info particles_global_info;
        part_simple_dataType particles_dataType;


        particles_global_info.mass=0.1;
        particles_global_info.relativistic=false;
        set_parts_typename(&particles_global_info,"part_simple");

Particles<part_simple,part_simple_info,part_simple_dataType> parts;
parts.initialize(particles_global_info,particles_dataType,&lat_part,boxSize);


        //cout<<"init done"<<endl;

        part_simple part;

        long index =0;
        /*
        for(int i=0;i<numparts;i++)
            for(int j=0;j<numparts;j++)
                for(int k=0;k<numparts;k++){
                    part.ID=index;
                    part.pos[0]= (Real)i * (Real)boxSize[0] / (Real)numparts;
                    part.pos[1]= (Real)j * (Real)boxSize[1] / (Real)numparts;
                    part.pos[2]= (Real)k * (Real)boxSize[2] / (Real)numparts;
                    part.vel[0]=1.0;
                    part.vel[1]=1.0;
                    part.vel[2]=1.0;
                    //part.mass=0.22;
                    parts.addParticle_global(part);
                    index++;
        }
        */
        int ratio = numparts/npts;
        //ratio;

//        Site xp(lat_part);
//        for(xp.first();xp.test();xp.next())
//        {
//            for(int i=0;i<ratio;i++)
//                for(int j=0;j<ratio;j++)
//                    for(int k=0;k<ratio;k++){

//            part.ID=index;
//            part.pos[0]= (Real)xp.coord(0) * (Real)boxSize[0] / (Real)npts;
//            part.pos[1]= (Real)xp.coord(1) * (Real)boxSize[1] / (Real)npts;
//            part.pos[2]= (Real)xp.coord(2) * (Real)boxSize[2] / (Real)npts;
//            part.vel[0]=0.01;
//            part.vel[1]=0.01;
//            part.vel[2]=0.01;
//            //part.mass=0.22;
//            parts.addParticle_global(part);
//            index++;
//            }
//        }

double r = 0.25;
int n_p = 7;
//int i;
int num;
double un;
int i;

//for(i=0;i<n_p;i++){
//parts.count(&num);
//parallel.sum(num);
//part.ID=num;
//un = (i-n_p/2)/(double)n_p;
//part.pos[0]= r*(1-un*un);
//part.pos[1]= 0.5;
//part.pos[2]= r*un;
//part.vel[0]=0.0;
//part.vel[1]=0.2;
//part.vel[2]=0.0;
//parts.addParticle_global(part);}


part.ID=0;
part.pos[0]= 0.5;
part.pos[1]= 0.52;
part.pos[2]= 0.4;
part.vel[0]=0.2;
part.vel[1]=0.0;
part.vel[2]=0.0;
parts.addParticle_global(part);
part.ID=1;
part.pos[0]= 0.5;
part.pos[1]= 0.57;
part.pos[2]= 0.45;
part.vel[0]=0.2;
part.vel[1]=0.0;
part.vel[2]=0.0;
parts.addParticle_global(part);
part.ID=2;
part.pos[0]= 0.5;
part.pos[1]= 0.6;
part.pos[2]= 0.5;
part.vel[0]=0.2;
part.vel[1]=0.0;
part.vel[2]=0.0;
parts.addParticle_global(part);
part.ID=3;
part.pos[0]= 0.5;
part.pos[1]= 0.57;
part.pos[2]= 0.55;
part.vel[0]=0.2;
part.vel[1]=0.0;
part.vel[2]=0.0;
parts.addParticle_global(part);
part.ID=4;
part.pos[0]= 0.5;
part.pos[1]= 0.52;
part.pos[2]= 0.6;
part.vel[0]=0.2;
part.vel[1]=0.0;
part.vel[2]=0.0;
parts.addParticle_global(part);

//        cout<<"implementation done"<<endl;

//parts.saveHDF5("bench_part",1);

//        parts.loadHDF5("bench_part",1);


double tau=0;
double tau_end =20;
double dtau=0.1;
int n_time=(int)(tau_end/dtau)+1;
int cycle = 0;
double pos[n_time][3];

int world_rank;
MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
int world_size;
MPI_Comm_size(MPI_COMM_WORLD, &world_size);
MPI_Barrier(MPI_COMM_WORLD);

//  FILE* Result;
//  Result=fopen("path","w");
//	char text[20];


while (true)    
{		if (tau > tau_end ) break; // simulation complete

         Site xpart(parts.lattice());   
         std::list<part_simple>::iterator it;
             
         for(xpart.first();xpart.test();xpart.next())
         {
             if(parts.field()(xpart).size!=0)
             {
                 for(it=parts.field()(xpart).parts.begin();it != parts.field()(xpart).parts.end();++it)
                 {
//                     cout<< "argarg: "<< (*it).ID <<endl;
//											index2++;
cout<<(*it).ID<<" "<<tau<<" "<<(*it).pos[0]<<" "<<(*it).pos[1]<<" "<<(*it).pos[2]<<endl;

//fprintf(Result,"%.7lf %.7lf %.7lf %.7lf\n",tau,(*it).pos[0],(*it).pos[1],(*it).pos[2] );

                 }

             }
         }

//        parts.updateVel(&updateVel_simple,1.0);
        parts.updateVel(&updateVel_round,dtau);

        parts.moveParticles(&move_particles_round,dtau);

tau += dtau;
cycle++;
}

////cout<<cycle<<endl;

//MPI_Barrier(MPI_COMM_WORLD);

//string H5FILE_NAME="path.h5",DATASETNAME="path";

//int RANK=2;
//hid_t       file, dataset; /* file and dataset handles */
//hid_t       datatype, dataspace,plist_id; /* handles */
//hsize_t     dimsf[2]; /* dataset dimensions */
//herr_t      status;

//plist_id = H5Pcreate(H5P_FILE_ACCESS);

//    file = H5Fcreate(H5FILE_NAME.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
//H5Pclose(plist_id);
//    /*
//     * Describe the size of the array and create the data space for fixed
//     * size dataset.
//     */
//    dimsf[0] = n_time;
//    dimsf[1] = 3;

//    dataspace = H5Screate_simple(RANK, dimsf, NULL);

//    /*
//     * Define datatype for the data in the file.
//     * We will store little endian INT numbers.
//     */
//    datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
//    status = H5Tset_order(datatype, H5T_ORDER_LE);

//    /*
//     * Create a new dataset within the file using defined dataspace and
//     * datatype and default dataset creation properties.
//     */
//    dataset = H5Dcreate2(file, DATASETNAME.c_str(), datatype, dataspace,
//			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

//    /*
//     * Write the data to the dataset using default transfer properties.
//     */
//    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos);

//    /*
//     * Close/release resources.
//     */
//    H5Sclose(dataspace);
//    H5Tclose(datatype);
//    H5Dclose(dataset);
//    H5Fclose(file);
}
