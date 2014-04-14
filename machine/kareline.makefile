FC=mpif90
FFLAGS= $(COMPILO) $(OMP) $(INCLUDE) $(HDF5) $(LAPACK) $(FFTW) $(GPROF) 
LIBS=$(LIBFFTW) $(LIB) $(LIBHDF5) $(LIBLAPACK) $(LIBOMP)
CC=mpicc
CCFLAGS=-O2
GPROF=-pg
LIBOMP=#-lfftw3_omp
OMP=-fopenmp
FFTWHOME=/home/marradi/CODE/fftw3

#Hard Prod
#COMPILO=-O2 -ipo -vec-report -xSSE4.2 
#Prod
COMPILO=-O3 -ffree-line-length-none
#COMPILO=-O2
#Debug
#COMPILO=-O2 -g -traceback -warn all -ftrapuv -fp-stack-check 
#COMPILO=-O2 -g -traceback -warn all -fp-stack-check 
#COMPILO=-O2 -g -traceback -check all -warn all
#hard debug
#COMPILO=-g -debug -traceback -implicitnone -warn all -fpe0 -fp-stack-check -ftrapuv -heap-arrays -gen-interface -warn interface
#COMPILO=-g -debug -traceback -check all -implicitnone -warn all -fpe0 -fp-stack-check -ftrapuv -heap-arrays -gen-interface -warn interface
INCLUDE=
LIB= 
HDF5=-I${PHDF5HOME}/include -DHDF5
LIBHDF5=-L${PHDF5HOME}/lib -lhdf5_fortran -lhdf5 -lz
LAPACK=-DLAPACK
LIBLAPACK=-llapack
FFTW=-I${FFTWHOME}/include
LIBFFTW=-L${FFTWHOME}/lib -lfftw3_omp -lfftw3 -lm

