FC=mpif90
#FFLAGS= $(COMPILO) $(INCLUDE) $(HDF5) $(LAPACK) $(FFTW)
#FFLAGS= $(COMPILO) $(INCLUDE) $(HDF5) $(FFTW)
FFLAGS= $(COMPILO) $(INCLUDE) $(HDF5) $(FFTW) -DBLOCKING_SEND -DBLOCKING_SEND_PLUS
LIBS=$(LIBFFTW) $(LIB) $(LIBHDF5) $(LIBLAPACK)
CC=mpicc
CCFLAGS=-O2 -g

COMPILO=-O2 -g -Wall -ffree-line-length-none -fbacktrace -fbounds-check
INCLUDE=-I/usr/include
LIB=-L/usr/lib64/gfortran/modules
HDF5=-I/opt/local/include -DHDF5
LIBHDF5=-lhdf5_fortran -lhdf5 -lz
LAPACK=#-DLAPACK
LIBLAPACK=#-llapack
FFTW=-I/opt/local/include
LIBFFTW=-lfftw3
