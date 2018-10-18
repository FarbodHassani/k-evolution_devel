#! /bin/bash

code="extractor.cpp"

mpic++ -o extractor $code -DHDF5 -I../ -I/home/gf/work/forsat/geneva_works/gevolution/code/hdf5-1.8.18/hdf5/include -L/home/gf/work/forsat/geneva_works/gevolution/code/hdf5-1.8.18/hdf5/lib -std=c++11 -lm -lhdf5

LD_LIBRARY_PATH=/home/gf/work/forsat/geneva_works/gevolution/code/hdf5-1.8.18/hdf5/lib
export LD_LIBRARY_PATH

mpirun -np 4 ./extractor -m 2 -n 2 -s "./files/lcdm_snap000_part" -o "./files/out" -i "./files/ids"
