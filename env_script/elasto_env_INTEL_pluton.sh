#!/bin/bash
# Remember: sourcing order between compiliers is NOT abelian 
#
# Ex:
# module load A
# source B
# has a different effect with respect to:
# source B
# module A

# TODO: scripting the module load and source research for latest version of compiler
module load openmpi/icc/mt/2.0.4.1
source /opt/intel/parallel_studio_xe_2018/compilers_and_libraries_2018.0.128/linux/bin/compilervars.sh intel64

export DIR="/home_nfs/marradil/Program/LIPHY/build_lib/BUILD"
export MKLPATH="${MKLROOT}/lib/intel64_lin"
export MKLINCLUDE="${MKLROOT}/include"
