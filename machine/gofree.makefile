FC=mpiifort
FFLAGS= $(COMPILO) $(INCLUDE) $(HDF5) $(LAPACK) $(FFTW)
LIBS=$(LIBFFTW) $(LIB) $(LIBHDF5) $(LIBLAPACK)

#INTEL
COMPILO=-O2
#COMPILO=-O2 -g -warn all -check all -traceback -trace
#COMPILO=-g -trace
#COMPILO=-ftz -ftrapuv -g -debug all -debug-parameters all -traceback -O0
INCLUDE=
LIB=-lmpi
#HDF5=-I/applis/ciment/stow/x86_64/hdf5-1.8.4-patch1/include -DHDF5
#LIBHDF5=-L/applis/ciment/stow/x86_64/hdf5-1.8.4-patch1/lib
LAPACK=
LIBLAPACK=
FFTW=-I/applis/ciment/stow/x86_64/fftw-3.2.2-intel-12/include
LIBFFTW= -L/applis/ciment/stow/x86_64/fftw-3.2.2-intel-12/lib -lfftw3

##GNU
#COMPILO=-O2
##INCLUDE=-I/applis/ciment/stow/x86_64/gcc-4.6.1/include -I/applis/ciment/stow/x86_64/openmpi-1.4.1/include
##LIB=-L/applis/ciment/stow/x86_64/gcc-4.6.1/lib -L/applis/ciment/stow/x86_64/openmpi-1.4.1/lib
##HDF5=-I/applis/ciment/stow/x86_64/hdf5-1.8.4-patch1/include -DHDF5
##LIBHDF5=-L/applis/ciment/stow/x86_64/hdf5-1.8.4-patch1/lib
##LAPACK=
##LIBLAPACK=
##FFTW=-I/applis/ciment/stow/x86_64/fftw-3.2.2/include
##LIBFFTW=-L/applis/ciment/stow/x86_64/fftw-3.2.2/lib -lfftw3
#
##GNU
#COMPILO=-O2
#INCLUDE=-I/applis/ciment/v2/stow/x86_64/gcc_4.6.2/openmpi_1.6.4/include
#LIB=-L/applis/ciment/v2/stow/x86_64/gcc_4.6.2/openmpi_1.6.4/lib
#HDF5=-I -DHDF5
#LIBHDF5= -lhdf5
#LAPACK=
#LIBLAPACK=
#FFTW=-I/applis/ciment/v2/stow/x86_64/gcc_4.6.2/fftw_3.3/include
#LIBFFTW=-I/applis/ciment/v2/stow/x86_64/gcc_4.6.2/fftw_3.3/lib
