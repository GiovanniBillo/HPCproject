# HPC project - 5-point stencil parallel implementation
In this repo I implemented a parallel version of the 5-point stencil algorithm to calculate the heat diffusion equation.
The code is organized as follows:
- `src/`: contains the main iteration loop, both the serial and parallel version (`stencil_template_parallel.c`).
- `include/`: contains helper functions for the algorithm per-se (in `stencil_template_parallel.h`) and helpers for time-tracking in the code(`timing.h).
- a Makefile to compile everything (with debug, profile and release options).
- scripts (`sunmit.sh` and `submit_improved.sh`) to submit the jobs directly.

A scalability study in 2 parts is then performed.
The result of the first part, scaling openMP threads, are contained in `thread_exp`.
The result of the second part, scaling nodes, are contained in `multinode_exp`.

All these results were obtained on the [DCGP partition of Cineca Leonardo](https://leonardo-supercomputer.cineca.eu/hpc-system/#jump-partition).

The approach was the following:
- Initially results are collected with "vanilla" settings, including an optimal Thread/Task ratio.
- Then, further results are collected and compared by modifying the openMP policy and parallelizing part of the initialization part (the initial factorization described in `stencil_template_parallel.h`, which compiles conditionally).   

Either python or shell scripts are present in the respective folders to obtain the tables and figures used in the presentation (also present on the repository and made with [Advanced Slides](https://mszturc.github.io/obsidian-advanced-slides/)).

All experiments can be repeated by running the scripts `submit.sh` (for "vanilla" configurations) and 'submit_improved.sh` for modified configurations.
