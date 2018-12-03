ELASTO
======
Large scale parallelized 3d mesoscopic simulations of the mechanical response to shear in disordered media

## Todo
- Luca:	   put a tree representation of the folders
- Kirsten: make physical introduction of the model (Ref to papers, etc etc) 

## Installation

We define the following variables:
	- ROOT=<where you have ELASTO folder>

The parallel version of ELASTO can be compiled both on INTEL or GNU compiler. The env\_script folder contains
couple of script:

	- elasto\_intel.sh
	- elasto\_gnu.sh

both prepare the proper compiling ELASTO enviromnent if we wish to use respectively the intel or gnu compiler. 

INTEL:
```
source $ROOT/env\_script/elasto\_intel.sh
cd $ROOT/builddir
make 
```
The executable filename elasto.exe is finally created in the $ROOT/builddir folder ready to be run.


