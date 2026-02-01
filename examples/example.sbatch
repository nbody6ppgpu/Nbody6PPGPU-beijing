#!/bin/bash

#SBATCH --account=madnuc
#SBATCH --job-name=1m20sb
#SBATCH --mail-user=kaiwu.astro@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err
#SBATCH --nodes=1
#SBATCH --gres=gpu:4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH --cuda-mps # multiple MPI on a single GPU; option name differs among clusters
#SBATCH --partition=booster

set -x

module purge
module load Stages/2024 GCC HDF5/1.14.2-serial ParaStationMPI git
module list

EXEC=$(readlink -f $(ls -rt ../Nbody6PPGPU-beijing/build/nbody6++.* | tail -1))

if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  omp_threads=$SLURM_CPUS_PER_TASK
else
  omp_threads=1
fi
export OMP_NUM_THREADS=$omp_threads
ulimit -s unlimited
ulimit -c unlimited
export OMP_STACKSIZE=20G
mpiprocs=$(($SLURM_NTASKS_PER_NODE * $SLURM_NNODES))
echo " OMP_NUM_THREADS=" $OMP_NUM_THREADS " mpiprocs=" $mpiprocs " OMP_STACKSIZE=" $OMP_STACKSIZE

export CUDA_VISIBLE_DEVICES=0,1,2,3
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=true
export OMP_PLACES=threads
input_file_fullname=`ls -rt *.inp | tail -1`
input_file_basename=`basename $input_file_fullname .inp`
time srun -n $mpiprocs $EXEC < $input_file_fullname 1> $input_file_basename.$SLURM_JOB_ID.out 2> $input_file_basename.$SLURM_JOB_ID.err
