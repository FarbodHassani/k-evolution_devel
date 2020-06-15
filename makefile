COMPILER     := mpic++
<<<<<<< HEAD
INCLUDE      := -I../LATfield2/ -I./../class_public-2.7.1/include
LIB          := -L./../class_public-2.7.1 -lclass -lfftw3 -lm -lhdf5 -lgsl -lgslcblas

# target and source
EXEC         := kevolution
=======
INCLUDE      := -I/usr/local/Cellar/gsl/2.6/include  -I/Users/farbod/Dropbox/Projects/Blowup-EFT-kessence-BH/LATfield2-master/  -I /Users/farbod/Packages/hdf5-1.10.1/hdf5/include -I/Users/farbod/include    # add the path to LATfield2 and other libraries (if necessary)
LIB          :=-L/Users/farbod/Packages/hdf5-1.10.1/hdf5/lib -L/Users/farbod/usr/lib -L/usr/local/Cellar/gsl/2.6/lib -lfftw3 -lm -lhdf5 -lgsl -lgslcblas

# target and sources
EXEC         := gevolution
>>>>>>> Blowup_k-evolution
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
<<<<<<< HEAD
=======
DGEVOLUTION  += -DBACKREACTION_TEST
>>>>>>> Blowup_k-evolution
DGEVOLUTION  += -DEXACT_OUTPUT_REDSHIFTS
#DGEVOLUTION  += -DHAVE_HEALPIX
#DGEVOLUTION += -DCHECK_B
DGEVOLUTION += -DHAVE_CLASS # requires OPT -fopenmp and LIB -lclass
# requires OPT -fopenmp and LIB -lclass

# further compiler options
OPT          := -O3 -std=c++11 -fopenmp

$(EXEC): $(SOURCE) $(HEADERS) makefile
	$(COMPILER) $< -o $@ $(OPT) $(DLATFIELD2) $(DGEVOLUTION) $(INCLUDE) $(LIB)
