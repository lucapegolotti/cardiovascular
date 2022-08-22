#/bin/bash
module purge
module load cpu
module load shared
module load slurm
module load sdsc
module load gcc/9.2.0
module load cmake
module load openmpi
export SLURM_MPI_TYPE=pmi2

mkdir build
cd build
cmake ..
make -j 2
