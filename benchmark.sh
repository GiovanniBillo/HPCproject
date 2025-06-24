#!/bin/bash

# Parameters to test
X_SIZE=1000
Y_SIZE=1000
ITERATIONS=500
NSOURCES=5
ENERGY_PER_SOURCE=1.0

# Threads to test for OpenMP
THREADS_LIST=(1 2 4 8)

# Check if periodic flag is set
PERIODIC_FLAG=""
if [[ "$1" == "-p" ]]; then
    PERIODIC_FLAG="-p 1"
    echo "Running with periodic boundary conditions enabled"
fi

# Clean and build
echo "Rebuilding executables..."
make clean > /dev/null
make openmp > /dev/null  # Only builds stencil_serial and stencil_serial_fr

echo "---------------------------------------------------"
echo "Benchmarking SERIAL version"
echo "---------------------------------------------------"
SERIAL_TIME=$(./stencil_serial_fr -x $X_SIZE -y $Y_SIZE -n $ITERATIONS -e $NSOURCES -E $ENERGY_PER_SOURCE $PERIODIC_FLAG | grep "Execution time" | awk '{print $3}')
echo "Serial time: $SERIAL_TIME s"

echo ""
echo "---------------------------------------------------"
echo "Benchmarking PARALLEL version (OpenMP)"
echo "---------------------------------------------------"
for threads in "${THREADS_LIST[@]}"; do
    export OMP_NUM_THREADS=$threads
    PAR_TIME=$(./stencil_serial -x $X_SIZE -y $Y_SIZE -n $ITERATIONS -e $NSOURCES -E $ENERGY_PER_SOURCE $PERIODIC_FLAG | grep "Execution time" | awk '{print $3}')
    SPEEDUP=$(echo "scale=2; $SERIAL_TIME / $PAR_TIME" | bc)
    EFFICIENCY=$(echo "scale=2; 100 * $SPEEDUP / $threads" | bc)

    echo "Threads: $threads | Time: $PAR_TIME s | Speedup: $SPEEDUP | Efficiency: $EFFICIENCY%"
done

