#! /bin/bash

make
LD_LIBRARY_PATH=/home/farbod/packages/hdf5-1.8.17/hdf5/lib
export LD_LIBRARY_PATH

mpirun -np 4 ./gevolution -n 2 -m 2 -s settings.ini

rm gevolution
