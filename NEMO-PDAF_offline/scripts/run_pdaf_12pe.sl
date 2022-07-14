#!/bin/bash
#SBATCH --job-name=PDAF_NORDIC
###SBATCH -p smp 
#SBATCH -p mpp
##SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=0:30:00
#SBATCH --mail-user=lars.nerger@awi.de
#SBATCH --mail-type=END
#
#
source ${MODULESHOME}/init/bash

module use /global/AWImodules/
module purge
module load intel.compiler/2021.3.0
module load intel.mpi/2021.3.0
module load netcdf-parallel/4.8.1_intel2021

ulimit -s unlimited
export LD_LIBRARY_PATH=/global/AWIsoft/hdf5-parallel/1.12.1_intel2021/lib:${LD_LIBRARY_PATH}":"
echo $LD_LIBRARY_PATH

wdir=`pwd`
echo ' '
echo 'Run directory: ' $wdir
export wdir

export OMP_NUM_THREADS=3
echo $OMP_NUM_THREADS

srun ./PDAF_offline 


