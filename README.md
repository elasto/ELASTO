ELASTO
======
Large scale parallelized 3d mesoscopic simulations of the mechanical response to shear in disordered media

## Todo
- Luca:	   put a tree representation of the folders
- Kirsten: make physical introduction of the model (Ref to papers, etc etc) 

## Requirement
```
	- FFTW: fftw-3.3.8.tar.gz
	- HDF5: hdf5-1.8.16.tar
	- ZLIB: zlib-1.2.11.tar.gz
```
In the ** .bashrc ** file add the line ** ulimit -s unlimit **. Such a command define, for each mpi process, 
an infinite stack. By default a process stack is fixed to 8MB.   

## Building
We define the following variables:
```
	- ROOT=<where you have ELASTO folder>
```

The parallel version of ELASTO can be compiled both on INTEL or GNU compiler. The env_script folder contains
couple of script:

	- elasto_intel.sh
	- elasto_gnu.sh

both prepare the proper compiling ELASTO enviromnent if we wish to use respectively the intel or gnu compiler. 

INTEL:
```
source $ROOT/env_script/elasto_intel.sh
cd $ROOT/builddir
make 
```
The executable filename elasto.exe is finally created in the $ROOT/builddir folder ready to be run.

## How to create Tags
It is possible to create the tags to browse the ELASTO code more efficiently.

Required:
```
	ctags
```
Vi/Vim:
```
	cd $ROOT
	ctags -R .
```
then a file called tags should appear in $ROOT folder. Once you open Vi/Vim set the tag as follow:
```
	:set tags=tags
```
press enter and hereafter when you move your cursor on a function you can jump directly to its 
definition typing:
```
	Ctrl+]
```
and jumping back:
```
	Ctrl+t
```

## Profiling
In the corresponding machinefile rebuild the binary with:
```
 GPTOF=-g
 COMPILO=-O0 
```
The previous option force the compiler to generate debug symbol (-g) and don't optimize code.
The last compiler action is necessary to most of the profiler to create one-to-one map between 
the source code and the hotfunction code part.
 

### Maqao
Download the compiled version for your machine arch from: (http://www.maqao.org/)
The *Documentation* section explain the different kind of profiling that maqao 
offers.

Below an example taken from the mo cluster:
1. collect the data running the maqao profiler
```
$MAQAO lprof --mpi-command="mpirun -np 4 -hosts mo78" -- $BIN -i ../input
```

2. generate the html view of the hotspots and loops:
```
$MAQAO lprof xp=$PWD/<maqao_profiler_folder> of=html
```

