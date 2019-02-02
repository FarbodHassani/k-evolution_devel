#! /bin/bash

make
rm ./output/lcdm_background.dat
#mpic++ -DSINGLE main.cpp -w -o gevolution -O3 -std=c++11 -DFFT3D -DHDF5  -DPHINONLINEAR -I/home/farbod/packages/LATfield2 -I/home/farbod/packages/hdf5-1.8.17/hdf5/include -L/home/farbod/packages/hdf5-1.8.17/hdf5/lib -L/home/farbod/packages/class_public-2.5.0 -I/home/farbod/packages/class_public-2.5.0/include  -lfftw3f -lm -lhdf5 -lgsl -lgslcblas $1 $2


LD_LIBRARY_PATH=/home/farbod/Packages/packages/hdf5-1.8.17/hdf5/lib
export LD_LIBRARY_PATH

mpirun -np 8 ./gevolution -n 4 -m 2 -s settings.ini


rm gevolution

