#ifndef TIMING_H
#define TIMING_H

#include <mpi.h>
#include <omp.h>

// Global time trackers (defined in main.c)
extern double comm_time;
extern double comp_time;

// Macro to wrap MPI calls and accumulate communication time
#define TIME_MPI_CALL(call, time_var)      \
    do {                                   \
        double t0 = MPI_Wtime();           \
        call;                              \
        time_var += MPI_Wtime() - t0;      \
    } while(0)

// Macro to wrap compute sections and accumulate compute time
#define TIME_OMP_BLOCK(code_block)                   \
    do {                                             \
        double _t0 = omp_get_wtime();                \
        code_block                                   \
        comp_time += omp_get_wtime() - _t0;          \
    } while (0)

#endif // TIMING_H
