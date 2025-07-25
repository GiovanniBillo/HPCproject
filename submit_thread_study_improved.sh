#!/bin/bash
# === Load common SLURM options ===
source common_slurm_options.sh  # defines COMMON_OPTS

make clean
make -j$(nproc) OPENMP_SCHEDULE=1 PURGE_SERIAL=1   # Corresponds to -DOPTION1 -DOPTION2=value

# Submit the thread scaling study
echo "Submitting thread scalability study..."
sbatch \
  $COMMON_OPTS \
  $EXTRA_OPTS \
  --job-name=thread_study \
  --nodes=1 \
  --cpus-per-task=112 \
  --ntasks-per-node=1 \
  --mem=0 \
  --export=ENV=$ENV \
  --output=thread_scaling_%j.out \
  thread_study.sh

echo "Thread study submitted. Results will be in thread_results.csv and thread_speedup.csv"
