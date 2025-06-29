#!/bin/bash
set -x

# Initialize environment and log file
LOG_FILE="out_${NODES}nodes.txt"
ERROR_FILE="err_${NODES}nodes.txt"
STRONG_TIME_FILE="time_${NODES}_STRONG.txt"
WEAK_TIME_FILE="time_${NODES}_WEAK.txt"
STRONG_RESULTS="multinode_results_STRONG.csv"
WEAK_RESULTS="multinode_results_WEAK.csv"

# Clear previous log files
> "$LOG_FILE"
> "$ERROR_FILE"

# Set environment
if [[ $ENV == "LEONARDO" ]]; then
    echo "[$(date)] Setting environment for LEONARDO" >> "$LOG_FILE"
    module load openmpi/4.1.6--gcc--12.2.0
elif [[ $ENV == "ORFEO" ]]; then
    echo "[$(date)] Setting environment for ORFEO" >> "$LOG_FILE"
    module load openMPI
fi

EXEC=./stencil_parallel

# Common parameters
X_BASE=10000
Y_BASE=10000
total_tasks=$((NODES * BEST_TASKS))

# Set OpenMP parameters
export OMP_NUM_THREADS=$BEST_THREADS
export OMP_PLACES=cores
export OMP_PROC_BIND=close

echo "=============================================================" >> "$LOG_FILE"
echo "[$(date)] Starting experiment with $NODES nodes" >> "$LOG_FILE"
echo "Configuration:" >> "$LOG_FILE"
echo "  - MPI tasks per node: $BEST_TASKS" >> "$LOG_FILE"
echo "  - Threads per task: $BEST_THREADS" >> "$LOG_FILE"
echo "  - Total MPI tasks: $total_tasks" >> "$LOG_FILE"
echo "=============================================================" >> "$LOG_FILE"

# Strong scaling test
echo "[$(date)] Starting STRONG SCALING test" >> "$LOG_FILE"
x=$((X_BASE * NODES))
y=$((Y_BASE * NODES))

echo "  Problem size: ${x}x${y} (scales with nodes)" >> "$LOG_FILE"
echo "  Running command: mpirun -np $total_tasks $EXEC -x $x -y $y -o 1" >> "$LOG_FILE"

start_time=$(date +%s.%N)
/usr/bin/time -f "%e" -o "$STRONG_TIME_FILE" mpirun -np $total_tasks $EXEC -x $x -y $y -o 1 >> "$LOG_FILE" 2>> "$ERROR_FILE"
strong_exit_code=$?
end_time=$(date +%s.%N)
elapsed=$(echo "$end_time - $start_time" | bc)

if [[ $strong_exit_code -eq 0 ]]; then
    strong_time=$(cat "$STRONG_TIME_FILE")
    echo "[$(date)] STRONG SCALING completed successfully in ${strong_time}s (${elapsed}s wall time)" >> "$LOG_FILE"
    echo "$NODES,$x,$y,$strong_time,$elapsed" >> "$STRONG_RESULTS"
else
    echo "[$(date)] ERROR in STRONG SCALING test (exit code: $strong_exit_code)" >> "$LOG_FILE"
    echo "$NODES,$x,$y,ERROR,$strong_exit_code" >> "$STRONG_RESULTS"
fi

# Weak scaling test
echo "[$(date)] Starting WEAK SCALING test" >> "$LOG_FILE"
x=$X_BASE
y=$Y_BASE

echo "  Problem size: ${x}x${y} (fixed)" >> "$LOG_FILE"
echo "  Running command: mpirun -np $total_tasks $EXEC -x $x -y $y -o 1" >> "$LOG_FILE"

start_time=$(date +%s.%N)
/usr/bin/time -f "%e" -o "$WEAK_TIME_FILE" mpirun -np $total_tasks $EXEC -x $x -y $y -o 1 >> "$LOG_FILE" 2>> "$ERROR_FILE"
weak_exit_code=$?
end_time=$(date +%s.%N)
elapsed=$(echo "$end_time - $start_time" | bc)

if [[ $weak_exit_code -eq 0 ]]; then
    weak_time=$(cat "$WEAK_TIME_FILE")
    echo "[$(date)] WEAK SCALING completed successfully in ${weak_time}s (${elapsed}s wall time)" >> "$LOG_FILE"
    echo "$NODES,$x,$y,$weak_time,$elapsed" >> "$WEAK_RESULTS"
else
    echo "[$(date)] ERROR in WEAK SCALING test (exit code: $weak_exit_code)" >> "$LOG_FILE"
    echo "$NODES,$x,$y,ERROR,$weak_exit_code" >> "$WEAK_RESULTS"
fi

# Final summary
echo "=============================================================" >> "$LOG_FILE"
echo "[$(date)] Experiment completed for $NODES nodes" >> "$LOG_FILE"
echo "Results summary:" >> "$LOG_FILE"
echo "  - Strong scaling: ${strong_time:-ERROR}s" >> "$LOG_FILE"
echo "  - Weak scaling: ${weak_time:-ERROR}s" >> "$LOG_FILE"
echo "=============================================================" >> "$LOG_FILE"
##!/bin/bash
#set -x 
#if [[ $ENV -eq "LEONARDO" ]]; then
#	echo "Set environment for LEONARDO" >> out_${NODES}nodes.txt 
#        module load openmpi/4.1.6--gcc--12.2.0
#elif [[ $ENV -eq "ORFEO" ]]; then
#	echo "Set environment for ORFEO" >> out_${NODES}nodes.txt 
#        module load openMPI
#fi
#EXEC=./stencil_parallel

#module load openmpi/4.1.6--gcc--12.2.0
#X_BASE=10000
#Y_BASE=10000

#x=$(($X_BASE * $NODES))
#y=$(($Y_BASE * $NODES))

#total_tasks=$(($NODES * $BEST_TASKS))

#export OMP_NUM_THREADS=$BEST_THREADS
#export OMP_PLACES=cores
#export OMP_PROC_BIND=close

#echo "Running multinode scaling with $NODES nodes ($total_tasks MPI tasks, $BEST_THREADS threads/task), input size ${x}x${y}" >> out_${NODES}nodes.txt 

#echo "Performing STRONG SCALING..." >> out_${NODES}nodes.txt 

#/usr/bin/time -f "%e" -o time_tmp_${NODES}_STRONG.txt mpirun -np $total_tasks $EXEC -x $x -y $y -o 1 \
#  >> out_${NODES}nodes.txt 2>> err_${NODES}nodes.txt

## mpirun -np $total_tasks $EXEC -x $x -y $y -o 1 \
##   >> out_${NODES}nodes.txt 2 >> err_${NODES}nodes.txt
## srun --ntasks=$total_tasks --ntasks-per-node=$BEST_TASKS \
##      --cpus-per-task=$BEST_THREADS \
##      $EXEC -x $x -y $y -o 1 \
##      >> out_${NODES}nodes.txt 2>> err_${NODES}nodes.txt

#if [[ $? -eq 0 ]]; then
#  echo "$NODES,$x,$y,$(cat time_tmp.txt)" >> multinode_results_STRONG.csv
#else
#  echo "$NODES,$x,$y,ERROR" >> multinode_results_STRONG.csv
#fi 
#echo "Finished STRONG SCALING logging"
#echo "-----------------------------------------------------------------------------"
#echo "Performing WEAK SCALING..." >> out_${NODES}nodes.txt 
#x=10000
#y=10000
#echo "Set x=$x, y=$y for every node combination." >> out_${NODES}nodes.txt 

#/usr/bin/time -f "%e" -o time_tmp_${NODES}_STRONG.txt mpirun -np $total_tasks $EXEC -x $x -y $y -o 1 \
#  >> out_${NODES}nodes.txt 2>> err_${NODES}nodes.txt

## mpirun -np $total_tasks $EXEC -x $x -y $y -o 1 \
##   >> out_${NODES}nodes.txt 2 >> err_${NODES}nodes.txt
## srun --ntasks=$total_tasks --ntasks-per-node=$BEST_TASKS \
##      --cpus-per-task=$BEST_THREADS \
##      $EXEC -x $x -y $y -o 1 \
##      >> out_${NODES}nodes.txt 2>> err_${NODES}nodes.txt

#if [[ $? -eq 0 ]]; then
#  echo "$NODES,$x,$y,$(cat time_tmp.txt)" >> multinode_results_STRONG.csv
#else
#  echo "$NODES,$x,$y,ERROR" >> multinode_results_STRONG.csv
#fi 

#echo "Finished WEAK SCALING logging"
