#! /bin/bash

code="rw.cpp"

mpic++ -o rw $code -DHDF5 -I../ -I/home/gf/work/forsat/geneva_works/gevolution/code/hdf5-1.8.18/hdf5/include -L/home/gf/work/forsat/geneva_works/gevolution/code/hdf5-1.8.18/hdf5/lib -std=c++11 -lm -lhdf5

LD_LIBRARY_PATH=/home/gf/work/forsat/geneva_works/gevolution/code/hdf5-1.8.18/hdf5/lib
export LD_LIBRARY_PATH

mpirun -np 4 ./rw -s ./files/ic -o ./files/out3
