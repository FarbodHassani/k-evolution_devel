#! /bin/bash

#mpic++ -o getStart gettingStarted.cpp -I../
# -std=c++11

LD_LIBRARY_PATH=/home/farbod/packages/hdf5-1.8.17/hdf5/lib
export LD_LIBRARY_PATH
mpic++ -o mfour ./Wave_LatfielFourier.cpp -I../ -DFFT3D -DHDF5   -I/home/farbod/packages/hdf5-1.8.17/hdf5/include -I/home/farbod/packages/LATfield2 -L/home/farbod/packages/hdf5-1.8.17/hdf5/lib -DHDF5  -lfftw3 -lm -lhdf5 -lgsl -lgslcblas -O3 -std=c++11 
mpirun -np 4 ./mfour -n 2  -m 2
rm mfour
#mpirun -np $4 ./poissonSolver.cpp -n $2 -m $2
python2.7 Four_python_figs.py
