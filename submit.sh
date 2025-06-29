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
# === Read best_n_threads and best_n_tasks ===
if [[ ! -f best_thread_task ]]; then
  echo "❌ Missing 'best_thread_task' file" >&2
  exit 1
fi

read BEST_THREADS BEST_TASKS < <(awk -F, 'NR==1{next} {print $1, $2}' best_thread_task)
echo "✅ Read best config: $BEST_TASKS MPI tasks × $BEST_THREADS threads"

# === Submit job array for multinode scaling ===
NODE_COUNTS=(1 2 4 8 16)

for i in "${!NODE_COUNTS[@]}"; do
  NODES="${NODE_COUNTS[$i]}"
  # NTASKS=$((NODES*BEST_TASKS))
  echo "Starting job with $NODES nodes and $BEST_TASKS tasks per node, with $BEST_THREADS cpus per task..."
  sbatch \
    $COMMON_OPTS \
    $EXTRA_OPTS \
    --nodes=$NODES \
    --ntasks-per-node=$BEST_TASKS \
    --cpus-per-task=$BEST_THREADS \
    --export=ENV=$ENV,BEST_THREADS=$BEST_THREADS,BEST_TASKS=$BEST_TASKS,NODES=$NODES \
    --output=log_multinode_${NODES}nodes.out \
    multinode_study.sh
done

