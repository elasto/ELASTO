FC=mpif90
FFLAGS= $(COMPILO) $(INCLUDE) $(HDF5) $(LAPACK) $(FFTW) $(GPROF)
LIBS=$(LIBFFTW) $(LIB) $(LIBHDF5) $(LIBLAPACK)
CC=mpicc
GPROF=-pg


#INTEL
#Hard Prod
#COMPILO=-O2 -ipo -vec-report -xSSE4.2
#Prod
COMPILO=-O3 -ffree-line-length-none -fbounds-check
INCLUDE=
#HDF5=-I/home/vollant3a/opt_x86/phdf5/include -DHDF5 -I/home/vollant3a/opt_x86/zlib/include
#HDF5=-I/home/marradi/CODES/HDF5/include -DHDF5 
HDF5=-I${hdf5_DIR}/include -DHDF5 
LIBHDF5=-L${hdf5_DIR}/lib -lhdf5_fortran -lhdf5
FFTW=-I${fftw_DIR}/include
LIBFFTW=-L${fftw_DIR}/lib -lfftw3 -lm
BLAS=-I${blas_DIR}/include
LAPACK=-I${lapack_DIR}/include -DLAPACK
LIBLAPACK=-L${lapack_DIR}/lib -llapack ${blas_DIR}/lib/blas_LINUX.a
CCFLAGS=-O2
LIB=-L${zlib_DIR}/lib -lz
