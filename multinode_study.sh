#!/bin/bash
set -x 
if [[ $ENV -eq "LEONARDO" ]]; then
	echo "Set environment for LEONARDO"
        module load openmpi/4.1.6--gcc--12.2.0
elif [[ $ENV -eq "ORFEO" ]]; then
	echo "Set environment for ORFEO"
        module load openMPI
EXEC=./stencil_parallel

X_BASE=10000
Y_BASE=10000

x=$(($X_BASE * $NODES))
y=$(($Y_BASE * $NODES))

total_tasks=$(($NODES * $BEST_TASKS))

export OMP_NUM_THREADS=$BEST_THREADS
export OMP_PLACES=cores
export OMP_PROC_BIND=close

echo "Running multinode scaling with $NODES nodes ($total_tasks MPI tasks, $BEST_THREADS threads/task), input size ${x}x${y}" >> out_${NODES}nodes.txt 

mpirun -np $total_tasks $EXEC -x $x -y $y -o 1 \
  >> out_${NODES}nodes.txt 2>> err_${NODES}nodes.txt
# mpirun -np $total_tasks $EXEC -x $x -y $y -o 1 \
#   >> out_${NODES}nodes.txt 2 >> err_${NODES}nodes.txt
# srun --ntasks=$total_tasks --ntasks-per-node=$BEST_TASKS \
#      --cpus-per-task=$BEST_THREADS \
#      $EXEC -x $x -y $y -o 1 \
#      >> out_${NODES}nodes.txt 2>> err_${NODES}nodes.txt

if [[ $? -eq 0 ]]; then
  echo "$NODES,$x,$y,$(cat time_tmp.txt)" >> multinode_results.csv
else
  echo "$NODES,$x,$y,ERROR" >> multinode_results.csv
fi 

