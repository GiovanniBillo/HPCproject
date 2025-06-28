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

## PARAMS
module load openMPI  # Make sure this loads the right OpenMP-aware MPI if needed

# Clean and recompile
make clean && make

EXEC=./stencil_parallel

# Output CSV file
echo "threads,time" > results.csv

# List of threads to test
for threads in 1 2 4 8 16 32 56 84 112; do
    export OMP_NUM_THREADS=$threads
    export OMP_PLACES=cores
    export OMP_PROC_BIND=close

    echo "Running with OMP_NUM_THREADS=$threads..."

    # Run the program and capture runtime in temp file
    /usr/bin/time -f "%e" -o time_tmp.txt $EXEC -x 1000 -y 1000 -o 1 > out_${threads}.txt 2> warn_${threads}.txt

    # Append to CSV: threads,time
    runtime=$(cat time_tmp.txt)
    echo "${threads},${runtime}" >> results.csv
done

# Compute speedup and efficiency
awk -F, 'NR==1 {print "threads,speedup,efficiency"} NR==2 {t1=$2} NR>2 {s=t1/$2; e=s/$1; printf("%s,%.2f,%.2f\n", $1, s, e)}' results.csv > speedup.csv


# module load openMPI

# ## clean compiled programs and redo everything within the machine
# make clean && make

# EXEC=./stencil_parallel

# echo "threads,time" > results.csv
# # List of threads to test
# for threads in 1 2 4 8 16 32 56 84 112; do
#     export OMP_NUM_THREADS=$threads
#     printf "%d," "$threads" >> results.csv
#     export OMP_PLACES=cores
#     export OMP_PROC_BIND=close

#     echo "Running with OMP_NUM_THREADS=$threads"
#     { /usr/bin/time -f "%e" $EXEC -x 100 -y 100 -o 1 > out_${threads}.txt 2> warn_${threads}.txt ; } 2>> results.csv

# done

# awk -F, 'NR==2 {t1=$2} NR>1 {s=t1/$2; e=s/$1; printf("%s,%.2f,%.2f\n", $1, s, e)}' results.csv > speedup.csv
# gnuplot -persist plot_speedup.gp

 
