# programming environment
COMPILER     := mpic++
INCLUDE      := -I/usr/local/include/gsl  -I./../../LATfield2-master  -I /Users/farbod/Packages/hdf5-1.10.1/hdf5/include -I./../../class_public-2.7.1/include/
LIB          := -L./../../class_public-2.7.1/ -L/Users/farbod/Packages/hdf5-1.10.1/hdf5/lib  -lclass -lfftw3 -lm -lhdf5 -lgsl -lgslcblas

# target and source
EXEC         := gevolution
SOURCE       := main.cpp
HEADERS      := $(wildcard *.hpp)

# mandatory compiler settings (LATfield2)
DLATFIELD2   := -DFFT3D -DHDF5

# optional compiler settings (LATfield2)
#DLATFIELD2   += -DH5_HAVE_PARALLEL
#DLATFIELD2   += -DEXTERNAL_IO

# optional compiler settings (gevolution)
DGEVOLUTION  := -DPHINONLINEAR
DGEVOLUTION  += -DBENCHMARK
DGEVOLUTION  += -DEXACT_OUTPUT_REDSHIFTS
#DGEVOLUTION  += -DCHECK_B
DGEVOLUTION  += -DHAVE_CLASS # requires OPT -fopenmp and LIB -lclass

# further compiler options
OPT          := -O3 -std=c++11 -fopenmp

$(EXEC): $(SOURCE) $(HEADERS) makefile
	$(COMPILER) $< -o $@ $(OPT) $(DLATFIELD2) $(DGEVOLUTION) $(INCLUDE) $(LIB)

