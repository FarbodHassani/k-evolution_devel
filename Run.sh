#! /bin/bash

rm ./output/*
make
LD_LIBRARY_PATH=/home/farbod/packages/hdf5-1.8.17/hdf5/lib
export LD_LIBRARY_PATH
export OMPI_MCA_rmaps_base_oversubscribe=1
export PATH=/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin:/Applications/Xcode.app/Contents/Developer/usr/bin:$PATH

mpirun -np 8 ./gevolution -n 4 -m 2 -s RiessSciama_test_cs7.ini

rm gevolution
