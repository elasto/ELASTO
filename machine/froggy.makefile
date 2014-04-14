FC=mpif90
FFLAGS= $(COMPILO) $(INCLUDE) $(HDF5) $(LAPACK) $(FFTW)
LIBS=$(LIBFFTW) $(LIB) $(LIBHDF5) $(LIBLAPACK)
CC=mpicc
CCFLAGS=-O2

#INTEL
#Hard Prod
#COMPILO=-O2 -ipo -vec-report -xSSE4.2
#Prod
COMPILO=-O2 -ffree-line-length-none
#Debug
#COMPILO=-O2 -g -traceback -check all -warn all -ffree-line-length-none 
#hard debug
#COMPILO=-g -debug -traceback -check all -implicitnone -warn all -fpe0 -fp-stack-check -ftrapuv -heap-arrays -gen-interface -warn interface -ffree-line-length-none 
INCLUDE=
LIB=
HDF5=-I/home/vollant3a/opt_x86/phdf5/include -DHDF5 -I/home/vollant3a/opt_x86/zlib/include
LIBHDF5=-L/home/vollant3a/opt_x86/phdf5/lib -lhdf5_fortran -lhdf5 -L/home/vollant3a/opt_x86/zlib/lib -lz
#FFTW=-I/applis/ciment/v2/stow/x86_64/intel_13.0.1/fftw_3.3/include
#LIBFFTW=-L/applis/ciment/v2/stow/x86_64/intel_13.0.1/fftw_3.3/lib -lfftw3
FFTW=-I/applis/site/stow/gcc_4.4.6/fftw_3.3.3/include
LIBFFTW=-L/applis/site/stow/gcc_4.4.6/fftw_3.3/lib -lfftw3
#FFTW=${CFLAGS}
#LIBFFTW=${LDFLAGS}

#GFORTRAN
#Hard Prod
#COMPILO=-O3 -g -ffree-line-length-none
#Prod
#COMPILO=-O2 -ffree-line-length-none 
#Debug
#COMPILO=-O2 -g -traceback -check all -warn all -ffree-line-length-none 
#hard debug
#COMPILO=-g -debug -traceback -check all -implicitnone -warn all -fpe0 -fp-stack-check -ftrapuv -heap-arrays -gen-interface -warn interface -ffree-line-length-none 
INCLUDE=
LIB=
#FFTW=-I/applis/ciment/v2/stow/x86_64/gcc_4.6.2/fftw_3.3/include
#LIBFFTW=-L/applis/ciment/v2/stow/x86_64/gcc_4.6.2/fftw_3.3/lib -lfftw3
