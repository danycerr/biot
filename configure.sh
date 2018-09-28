#!/bin/bash
source /u/sw/etc/profile
module load  gcc-glibc
module load  getfem
module load  qhull
module load superlu


#export LD_LIBRARY_PATH=/u/archive/agip/cerroni/software/lapack/lib64/:$LD_LIBRARY_PATH
#export PATH=/u/archive/agip/cerroni/software/qhull/qhull/bin:$PATH
#export LD_LIBRARY_PATH=/u/archive/agip/cerroni/software/qhull/qhull/lib:$LD_LIBRARY_PATH
#export GET_FEM_DIR=/u/archive/agip/cerroni/software/getfem5/
#export QHULL_DIR=/u/archive/agip/cerroni/software/qhull/qhull/
#export LAPACK_DIR=/u/archive/agip/cerroni/software/lapack/
export SAMG=/opt/lib/samg
export OMP_NUM_THREADS=2

export LD_LIBRARY_PATH=$SAMG:$LD_LIBRARY_PATH
export SVD_LICENSE_FILE=@nisserver.mate.polimi.iexport SVD_LICENSE_FILE=@nisserver.mate.polimi.it
