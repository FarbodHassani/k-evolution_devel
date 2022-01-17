# programming environment
COMPILER     := /usr/local/bin/mpic++
INCLUDE      := -I/usr/local/Cellar/fftw/3.3.10/include -I/usr/local/Cellar/gsl/2.6/include -I/usr/local/Cellar/hdf5/1.12.1/include -I./../LATfield2   # add the path to LATfield2 and other libraries (if necessary)
LIB          :=-L/usr/local/Cellar/hdf5/1.12.1/lib -L/usr/local/Cellar/gsl/2.6/lib  -L/usr/local/Cellar/fftw/3.3.10/lib -lfftw3 -lm -lhdf5 -lgsl -lgslcblas  #-lclass

# target and sources
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
DGEVOLUTION  += -DBACKREACTION_TEST
DGEVOLUTION  += -DEXACT_OUTPUT_REDSHIFTS
DGEVOLUTION  += -DMAX_OUTPUTS=64
#DGEVOLUTION += -DCHECK_B
#DGEVOLUTION += -DHAVE_CLASS # requires OPT -fopenmp and LIB -lclass

# further compiler options
OPT          := -O3 -std=c++11

$(EXEC): $(SOURCE) $(HEADERS) makefile
	$(COMPILER) $< -o $@ $(OPT) $(DLATFIELD2) $(DGEVOLUTION) $(INCLUDE) $(LIB)

