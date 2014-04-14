FC=mpif90
FFLAGS= $(COMPILO) $(INCLUDE) $(HDF5) $(LAPACK) $(FFTW)
LIBS=$(LIBFFTW) $(LIB) $(LIBHDF5) $(LIBLAPACK)
CC=mpicc
CCFLAGS=-O2 -g

COMPILO=-O2 -g -Wall -ffree-line-length-none -fbacktrace -fbounds-check
INCLUDE=-I/usr/include
LIB=-L/usr/lib64/gfortran/modules
HDF5=-I/opt/phdf5-1.8.10-for_gfortran/include -DHDF5
LIBHDF5=-L/opt/phdf5-1.8.10-for_gfortran/lib -lhdf5_fortran -lhdf5 -lz
LAPACK=-DLAPACK
LIBLAPACK=-l lapack
FFTW=-I/usr/include
LIBFFTW=-lfftw3
