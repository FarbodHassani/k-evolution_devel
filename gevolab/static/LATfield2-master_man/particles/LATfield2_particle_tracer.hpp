#ifndef LATFIELD2_PARTICLE_DEF_HPP
#define LATFIELD2_PARTICLE_DEF_HPP

/*! \file LATfield2_particle_tracer.hpp
 \brief Description of the particle type: "part_tracer"

 */


/*type table

 datatype      memType         fileType
 int            H5T_NATIVE_INT   INT_TYPE_FILE
 long           H5T_NATIVE_LONG  LONG_TYPE_FILE
 float          H5T_NATIVE_FLOAT  FLOAT_TYPE_FILE
 double         H5T_NATIVE_DOUBLE  DOUBLE_TYPE_FILE
 bool           H5T_NATIVE_HBOOL  BOOL_TYPE_FILE
 */




/**
 * \addtogroup partdesc
 * @{
 */

/**
 * \addtogroup parttracer "part_tracer"
 * @{
 */

// part_tracer definition

/*! \struct part_tracer
 \brief individual properties of the particle type "part_tracer".

 "part_tracer" contains the minimal list of individual properties is: ID, position,velocity,
 */
struct part_tracer{

  long ID;
  LATfield2::Real time;
  LATfield2::Real pos[3];
  LATfield2::Real vel[3];
};
/*!
 \brief overloading of the << operator for individual property strucutre.
 \return ostream containing the ID, position and velocity of the particle.
 */
ostream& operator<<(ostream& os, const part_tracer& p)
{
    os << "ID: "<<p.ID<<" , Pos: ("<< p.pos[0]<<","<< p.pos[1]<<","<< p.pos[2]<<") , Vel: (" << p.vel[0]<<","<< p.vel[1]<<","<< p.vel[2]<<")";
    return os;
}

/*! \struct part_tracer_info
 \brief individual properties of the particle type "part_tracer".

 "part_tracer" contains the minimal list of global properties is: type_name and size of type_name. But also contain two additionnal global properties: "mass" and "relativistic".
 */
struct part_tracer_info{
    double  mass;
    int relativistic;
    int type_name_size;
    char  type_name[64];

};


#ifdef HDF5
/*! \struct part_tracer_dataType
 \brief properties datatype structure.

 Structure which contain the HDF5 datatype of all properties.
 */



struct part_tracer_dataType{
  hid_t part_memType;
  hid_t part_fileType;
  hid_t part_info_memType;
  hid_t part_info_fileType;

  part_tracer_dataType(){

    hid_t strtype = H5Tcopy (H5T_C_S1);
    H5Tset_size (strtype, 64);


    part_memType = H5Tcreate(H5T_COMPOUND, sizeof (part_tracer));
    H5Tinsert(part_memType, "ID", HOFFSET (part_tracer, ID), H5T_NATIVE_LONG);
    H5Tinsert(part_memType, "time", HOFFSET (part_tracer, time), REAL_TYPE);
    H5Tinsert(part_memType, "positionX", HOFFSET (part_tracer, pos[0]), REAL_TYPE);
    H5Tinsert(part_memType, "positionY", HOFFSET (part_tracer, pos[1]), REAL_TYPE);
    H5Tinsert(part_memType, "positionZ", HOFFSET (part_tracer, pos[2]), REAL_TYPE);
    H5Tinsert(part_memType, "velocityX", HOFFSET (part_tracer, vel[0]), REAL_TYPE);
    H5Tinsert(part_memType, "velocityY", HOFFSET (part_tracer, vel[1]), REAL_TYPE);
    H5Tinsert(part_memType, "velocityZ", HOFFSET (part_tracer, vel[2]), REAL_TYPE);



    part_info_memType = H5Tcreate(H5T_COMPOUND, sizeof (part_tracer_info));
    H5Tinsert(part_info_memType, "mass", HOFFSET (part_tracer_info, mass), H5T_NATIVE_DOUBLE);
    H5Tinsert(part_info_memType, "relativistic", HOFFSET (part_tracer_info, relativistic), INT_TYPE_FILE);
    H5Tinsert(part_info_memType, "type_name_size", HOFFSET (part_tracer_info, type_name_size),INT_TYPE_FILE );
    H5Tinsert(part_info_memType, "type_name", HOFFSET (part_tracer_info, type_name), strtype);


    part_fileType = H5Tcreate (H5T_COMPOUND, sizeof(long) + 7 * sizeof(Real) );
    H5Tinsert(part_fileType, "ID"       ,0  ,LONG_TYPE_FILE);
    H5Tinsert(part_fileType, "time",sizeof(long)                    ,REAL_TYPE_FILE);
    H5Tinsert(part_fileType, "positionX",sizeof(long) + 1 * sizeof(Real) ,REAL_TYPE_FILE);
    H5Tinsert(part_fileType, "positionY",sizeof(long) + 2 * sizeof(Real) ,REAL_TYPE_FILE);
    H5Tinsert(part_fileType, "positionZ",sizeof(long) + 3 * sizeof(Real) ,REAL_TYPE_FILE);
    H5Tinsert(part_fileType, "velocityX",sizeof(long) + 4 * sizeof(Real) ,REAL_TYPE_FILE);
    H5Tinsert(part_fileType, "velocityY",sizeof(long) + 5 * sizeof(Real) ,REAL_TYPE_FILE);
    H5Tinsert(part_fileType, "velocityZ",sizeof(long) + 6 * sizeof(Real) ,REAL_TYPE_FILE);

    part_info_fileType = H5Tcreate(H5T_COMPOUND, sizeof(double) +(2*sizeof(int))+ 64);
    H5Tinsert(part_info_fileType, "mass", 0 ,DOUBLE_TYPE_FILE );
    H5Tinsert(part_info_fileType, "relativistic",sizeof(double), INT_TYPE_FILE);
    H5Tinsert(part_info_fileType, "type_name_size", sizeof(double)+sizeof(int),INT_TYPE_FILE );
    H5Tinsert(part_info_fileType, "type_name",sizeof(double)+(2*sizeof(int)), strtype);

    H5Tclose (strtype);
  }

};
#endif




/**@}*/

/**@}*/


#endif
