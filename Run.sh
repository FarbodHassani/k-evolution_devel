#! /bin/bash

rm ./output/*
module load foss/2018b HDF5 GSL
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH://home/hassani/LightCone-kessence/Healpix_3.50/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/hassani/LightCone-kessence/cfitsio

make
time srun -n 4 ./gevolution -n 2 -m 2 -s setting.ini 
#-p pk_ref.pre

