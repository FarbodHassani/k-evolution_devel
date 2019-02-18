# programming environment
COMPILER     := mpic++
INCLUDE      := -I/usr/local/include/gsl  -I/home/farbod/Dropbox/Projects/k-evolution/LATfield2-master  -I/home/farbod/Packages/packages/hdf5-1.8.17/hdf5/include -L/home/farbod/Packages/packages/hdf5-1.8.17/hdf5/lib  -I/home/farbod/Dropbox/Projects/k-evolution/Analysis_21Jan2019/Tests/class_public-2.6.3_IC_Gev/include
LIB          := -L/home/farbod/Dropbox/Projects/k-evolution/Analysis_21Jan2019/Tests/class_public-2.6.3_IC_Gev -L/home/farbod/Packages/packages/hdf5-1.8.17/hdf5/lib  -lfftw3 -lm -lhdf5 -lgsl -lgslcblas

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
OPT          := -O3 -std=c++11

$(EXEC): $(SOURCE) $(HEADERS) makefile
	$(COMPILER) $< -o $@ $(OPT) $(DLATFIELD2) $(DGEVOLUTION) $(INCLUDE) $(LIB)

