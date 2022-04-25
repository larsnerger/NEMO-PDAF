#!/bin/bash
#SBATCH --job-name=NEMO_NDC
####SBATCH -p mpp
##SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks=96
#SBATCH --ntasks-per-node=96
#SBATCH --time=0:30:00
#SBATCH --partition=standard96:test
#SBATCH --mail-user=lars.nerger@awi.de
#SBATCH --mail-type=END
#SBATCH -A zzz0002

module load intel/18.0.6
module load openmpi/intel/3.1.6
module load netcdf-parallel/ompi/intel.18/4.7.4

ulimit -s unlimited

export LD_LIBRARY_PATH=/sw/dataformats/netcdf-parallel/ompi/intel.18/4.7.4/skl/lib:/sw/dataformats/hdf5-parallel/ompi/intel.18/1.10.6/skl/lib:/sw/compiler/intel/compilers_and_libraries_2018.6.288/linux/compiler/lib/intel64_lin:/sw/compiler/intel/compilers_and_libraries_2018.6.288/linux/mkl/lib/intel64_lin

ls -l /sw/compiler/intel/compilers_and_libraries_2018.6.288/linux/mkl/lib/intel64_lin/libmkl_intel_lp64.so

np_nemo=92 #102    #number of PE's for nemo
np_xios=4    #number of PE's for xios

(( NCPU = np_nemo + np_xios ))
(( NEMO_LAST_CPU = np_nemo - 1 ))
(( XIOS_LAST_CPU = NCPU  - 1 ))
echo 0-${NEMO_LAST_CPU} ./nemo.exe > mpmd.conf 
echo ${np_nemo}-${XIOS_LAST_CPU} ./xios_server.exe >> mpmd.conf

srun -l --cpu_bind=cores --multi-prog mpmd.conf
