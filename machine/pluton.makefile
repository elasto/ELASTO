# Debug symbols
GPROF=-pg

#INTEL
FC=mpiifort
CC=mpiicc
COMPILO=-O3 -xAVX2 

# HDF5
HDF5=-I${DIR}/include -DHDF5 
LIBHDF5=-L${DIR}/lib -Wl,--start-group ${DIR}/lib/libhdf5.a ${DIR}/lib/libhdf5_hl.a ${DIR}/lib/libhdf5hl_fortran.a \
                                       ${DIR}/lib/libhdf5_fortran.a ${DIR}/lib/libz.a -Wl,--end-group 

# FFTW
FFTW=-I${DIR}/include
LIBFFTW=-L${DIR}/lib -lfftw3 -lm

# BLAS
BLAS=-I${MKLINCLUDE}
LAPACK=-I${MKLINCLUDE} -DLAPACK
LIBLAPACK=-Wl,--start-group ${MKLPATH}/libmkl_lapack95_ilp64.a     \
			    ${MKLPATH}/libmkl_blas95_ilp64.a -Wl,--end-group
# ZLIB
LIB=-L${MKLPATH} -Wl,--start-group ${MKLPATH}/libmkl_intel_lp64.a ${MKLPATH}/libmkl_intel_thread.a                          \
		 ${MKLPATH}/libmkl_core.a -Wl,--end-group -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_thread -lmkl_core -liomp5   \
                 -lpthread -lm -ldl -L${DIR}/lib -lz

# FLAGS and LIBS
INCLUDE= $(BLAS) $(LAPACK) $(FFTW) $(HDF5) $(GPROF)
CCFLAGS=$(COMPILO)
FFLAGS= $(COMPILO) $(INCLUDE) 
LIBS=$(LIBFFTW) $(LIB) $(LIBHDF5) $(LIBLAPACK)
