# programming environment
COMPILER     := /usr/local/bin/mpic++
INCLUDE      := -I/Users/farbod/Documents/GitHub/hi_class_pub_devel/include -I/usr/local/Cellar/fftw/3.3.10/include -I/usr/local/Cellar/gsl/2.7.1/include -I/usr/local/Cellar/hdf5/1.12.2/include -I./../LATfield2 # add the path to LATfield2 and other libraries (if necessary)
LIB          := -L/Users/farbod/Documents/GitHub/hi_class_pub_devel/ -L/usr/local/Cellar/hdf5/1.12.2/lib -L/usr/local/Cellar/gsl/2.7.1/lib  -L/usr/local/Cellar/fftw/3.3.10/lib -lfftw3 -lm -lhdf5 -lgsl -lgslcblas -lclass
# /Users/farbod/Documents/GitHub/hi_class_pub_devel
# /Users/farbod/Documents/GitHub/class_public-2.7.1/
# target and source
EXEC         := gevolution
SOURCE       := main.cpp
HEADERS      := $(wildcard *.hpp)

# mandatory compiler settings (LATfield2)
DLATFIELD2   := -DFFT3D -DHDF5

# optional compiler settings (LATfield2)
#DLATFIELD2   += -DH5_HAVE_PARALLEL
#DLATFIELD2   += -DEXTERNAL_IO # enables I/O server (use with care)
#DLATFIELD2   += -DSINGLE      # switches to single precision, use LIB -lfftw3f

# optional compiler settings (gevolution)
DGEVOLUTION  := -DPHINONLINEAR
DGEVOLUTION  += -DBENCHMARK
#DGEVOLUTION  += -DBACKREACTION_TEST
DGEVOLUTION  += -DEXACT_OUTPUT_REDSHIFTS
#DGEVOLUTION  += -DVELOCITY      # enables velocity field utilities
DGEVOLUTION  += -DCOLORTERMINAL
#DGEVOLUTION  += -DCHECK_B
DGEVOLUTION  += -DHAVE_HICLASS    # -DHAVE_HICLASS  or -DHAVE_CLASS requires LIB -lclass. The initial conditions are provided by hiclass!
DGEVOLUTION  += -DHAVE_HICLASS_BG    # -DHAVE_HICLASS requires LIB -lclass. The BG quantities are provided by hiclass and also parameters like c_s^2,w ...
#DGEVOLUTION  += -DHAVE_HEALPIX  # requires LIB -lchealpix

CDBG +=
CFLAGS += $(CDBG)

# further compiler options
OPT          := -O3 -std=c++11

$(EXEC): $(SOURCE) $(HEADERS) makefile
	$(COMPILER) $< -o $@ $(OPT) $(DLATFIELD2) $(DGEVOLUTION) $(INCLUDE) $(LIB)
	
lccat: lccat.cpp
	$(COMPILER) $< -o $@ $(OPT) $(DGEVOLUTION) $(INCLUDE)
	
lcmap: lcmap.cpp
	$(COMPILER) $< -o $@ $(OPT) -fopenmp $(DGEVOLUTION) $(INCLUDE) $(LIB) $(HPXCXXLIB)

clean:
	-rm -f $(EXEC) lccat lcmap

