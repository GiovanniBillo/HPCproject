#!/bin/bash
## original for just 1 MPI task per node
srun -A dssc -p EPYC --nodes=1 --ntasks-per-node=1  --cpus-per-task=112 --time=01:00:00 --pty bash

# now we try to allocate one big session for everything, in order to do a complete scalability study. 
# same number of tasks per node as cpus-per-task: we'll vary the quantity of each each time
# srun -A dssc -p EPYC --nodes=1 --ntasks-per-node=112  --cpus-per-task=112 --time=01:00:00 --pty bash
