FC=mpixlf2003
FFLAGS= $(OPTIM),$(OPTION)
#FFLAGS= $(DEBUG),$(OPTION)

#mpixlf2003 -c -O3 -g -WF,-DHDF5 advec_common_line.f90
OPTIM=-O3 -g -WF
#DEBUG = -O0 -g -C -qsigtrap -WF,-DHDF5
DEBUG=-O2 -g -qsigtrap -qoptdebug -qfullpath -C -qsigtrap -qkeepparm -WF,-DHDF5

HDF5=-DHDF5
#OPTION=$(HDF5),-DBLOCKING_SEND
OPTION=$(HDF5),-DBLOCKING_SEND,-DBLOCKING_SEND_PLUS

CC=mpixlc
CCFLAGS=-03 -g
