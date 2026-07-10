#!/bin/bash

# Set account, partition, mail, and site modules for your cluster if needed.
#SBATCH --job-name=nbody-benchmark
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err
#SBATCH --nodes=1
#SBATCH --gres=gpu:4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00

set -euo pipefail

EXEC=${NBODY_EXEC:?Set NBODY_EXEC to the NBODY6++GPU executable}
omp_threads=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=$omp_threads
export OMP_STACKSIZE=${OMP_STACKSIZE:-20G}
ulimit -s unlimited

mpiprocs=$((SLURM_NTASKS_PER_NODE * SLURM_NNODES))
input_file=$(find . -maxdepth 1 -name '*.inp' -print -quit)
test -n "$input_file"
output_file="${input_file%.inp}.${SLURM_JOB_ID}.out"
error_file="${input_file%.inp}.${SLURM_JOB_ID}.err"

echo "OMP_NUM_THREADS=$OMP_NUM_THREADS mpiprocs=$mpiprocs OMP_STACKSIZE=$OMP_STACKSIZE"
time srun -n "$mpiprocs" "$EXEC" < "$input_file" > "$output_file" 2> "$error_file"
grep -Fq 'END RUN' "$output_file"
