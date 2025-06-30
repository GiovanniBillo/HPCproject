#!/bin/bash
# === Load common SLURM options ===
source common_slurm_options.sh  # defines COMMON_OPTS
ENV=""

# === Detect cluster and set account/partition ===
if [[ "$PWD" == "/leonardo/home/userexternal/gbillo00"* ]]; then
  EXTRA_OPTS="--partition=dcgp_usr_prod -A uTS25_Tornator_0"
  ENV="LEONARDO"
  module load gcc/12.2.0
  module load openmpi/4.1.6--gcc--12.2.0
  echo "Detected Leonardo environment"
elif [[ "$PWD" == "/u/dssc/gbillo/HPCproject"* ]]; then
  EXTRA_OPTS="--partition=EPYC -A dssc"
  ENV="ORFEO"
  echo "Detected Orfeo environment"
  module load openMPI
else
  echo "❌ Unknown system environment (PWD=$PWD)" >&2
  exit 1
fi
make clean && make

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
