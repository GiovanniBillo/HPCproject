#!/bin/bash

# Load OpenMPI module
module load openMPI

# Compile the project
make clean && make
if [ $? -ne 0 ]; then
  echo "Compilation failed. Exiting."
  exit 1
fi

# Default values
NP=2
HOSTFILE="hostlist"
EXEC="./stencil_parallel"
X=1000
Y=1000
E=10
V=0
O=1
P=0  # periodic boundaries

# You can override these via CLI
while getopts "n:h:x:y:e:v:o:p:" opt; do
  case ${opt} in
    n) NP=${OPTARG} ;;
    x) X=${OPTARG} ;;
    y) Y=${OPTARG} ;;
    e) E=${OPTARG} ;;
    v) V=${OPTARG} ;;
    o) O=${OPTARG} ;;
    p) P=${OPTARG} ;;
    *) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
  esac
done

# Run the parallel program with chosen parameters
mpirun -np $NP $EXEC -x $X -y $Y -e $E -v $V -o $O -p $P

