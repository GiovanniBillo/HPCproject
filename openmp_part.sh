#!/bin/bash
# SBATCH --nodes=1
# SBATCH --cpus-per-task=112
# SBATCH --ntasks-per-node=1
# SBATCH --mem=0
# SBATCH --partition=dcgp_usr_prod
# SBATCH -A uTS25_Tornator_0
# SBATCH -t 00:30:00
# SBATCH --job-name=scalability_study

## in LEONARDO
# module load gcc/12.2.0
# module load openmpi/4.1.6--gcc--12.2.0
## in ORFEO
# salloc -A dssc -p EPYC --nodes=1 --cpus-per-task=112 --ntasks-per-node=1
module load openMPI

## clean compiled programs and redo everything within the machine
make clean && make

EXEC=./stencil_parallel

# List of threads to test
for threads in 1 2 4 8 16 32 56 84 112; do
    export OMP_NUM_THREADS=$threads
    export OMP_PLACES=cores
    export OMP_PROC_BIND=close

    echo "Running with OMP_NUM_THREADS=$threads"
    mpirun -np OMP_NUM_THREADS ${EXEC} -x 100 -y 100 > out_${threads}.txt
done

