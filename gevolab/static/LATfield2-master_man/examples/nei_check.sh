#! /bin/bash

code="nei_check.cpp"

mpic++ -o nei_check $code -I../ -std=c++11 -lm -lfftw3 -DFFT3D

mpirun -np 4 ./nei_check -m 2 -n 2
