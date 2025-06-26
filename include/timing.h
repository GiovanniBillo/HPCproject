#ifndef TIMING_H
#define TIMING_H

#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>

// ---- Global time trackers ---- //
extern double comm_time;
extern double comp_time;
extern double loop_start_time;
extern double loop_end_time;

#define MAX_THREADS 128  // Tune to your max expected thread count
extern double thread_times[MAX_THREADS]; // Thread-local timings

// ---- Macros ---- //
#define TIME_MPI_CALL(call, time_var)      \
    do {                                   \
        double t0 = MPI_Wtime();           \
        call;                              \
        time_var += MPI_Wtime() - t0;      \
    } while(0)

#define TIME_OMP_BLOCK(code_block)         \
    do {                                   \
        double _t0 = omp_get_wtime();      \
        code_block                         \
        comp_time += omp_get_wtime() - _t0;\
    } while (0)

#define START_LOOP_TIMER() (loop_start_time = MPI_Wtime())
#define STOP_LOOP_TIMER()  (loop_end_time = MPI_Wtime())

#define TIME_THREAD_BLOCK(code_block)                          \
    do {                                                       \
        int _tid = omp_get_thread_num();                       \
        double _t0 = omp_get_wtime();                          \
        code_block                                             \
        double _elapsed = omp_get_wtime() - _t0;               \
        thread_times[_tid] += _elapsed;                        \
    } while(0)

// ---- Timing Reporter ---- //
static inline void report_timing_stats(MPI_Comm comm, int Rank, int Ntasks, const char* label, int log_per_rank) {
    double total_time = loop_end_time - loop_start_time;
    double compute_time = total_time - comm_time;

    // Memory usage
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
    double mem_MB = rusage.ru_maxrss / 1024.0;

    // Reductions
    double global_comm_time, max_comm_time, avg_comm_time;
    MPI_Reduce(&comm_time, &global_comm_time, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    MPI_Reduce(&comm_time, &max_comm_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    avg_comm_time = global_comm_time / Ntasks;

    double global_comp_time, max_comp_time, min_comp_time;
    MPI_Reduce(&compute_time, &global_comp_time, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
    MPI_Reduce(&compute_time, &max_comp_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(&compute_time, &min_comp_time, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
    double avg_comp_time = global_comp_time / Ntasks;
    double imbalance_ratio = (max_comp_time - min_comp_time) / max_comp_time;

    double max_total_time;
    MPI_Reduce(&total_time, &max_total_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

    if (Rank == 0) {
        printf("\n=== %s Timing Report ===\n", label);
        printf("Wall clock time (max over ranks): %.6f s\n", max_total_time);
        printf("Total COMM time (sum): %.6f s\n", global_comm_time);
        printf("Max COMM time: %.6f s | Avg COMM time: %.6f s\n", max_comm_time, avg_comm_time);
        printf("Max COMP time: %.6f s | Min COMP time: %.6f s | Avg: %.6f s\n", max_comp_time, min_comp_time, avg_comp_time);
        printf("COMP imbalance ratio: %.2f%%\n", imbalance_ratio * 100);
    }

    // Thread-level timing summary (optional: only print if > 1 thread)
    int num_threads = omp_get_max_threads();
    if (num_threads > 1 && Rank == 0) {
        printf("--- Thread-level imbalance (per rank) ---\n");
    }

    double max_thread_time = 0.0, min_thread_time = 1e9;
    for (int i = 0; i < num_threads; ++i) {
        if (thread_times[i] > max_thread_time) max_thread_time = thread_times[i];
        if (thread_times[i] < min_thread_time) min_thread_time = thread_times[i];
    }

    if (num_threads > 1) {
        printf("Rank %d thread imbalance: max %.6f s, min %.6f s, ratio %.2f%%\n",
               Rank, max_thread_time, min_thread_time,
               100 * (max_thread_time - min_thread_time) / max_thread_time);
    }

    // Optional per-rank logging
    if (log_per_rank) {
        char fname[64];
        snprintf(fname, sizeof(fname), "timing_rank_%d.log", Rank);
        FILE *f = fopen(fname, "w");
        if (f) {
            fprintf(f, "RANK %d\n", Rank);
            fprintf(f, "Total time: %.6f s\n", total_time);
            fprintf(f, "Comm time: %.6f s\n", comm_time);
            fprintf(f, "Comp time: %.6f s\n", compute_time);
            fprintf(f, "Mem usage: %.2f MB\n", mem_MB);
            for (int i = 0; i < num_threads; i++) {
                fprintf(f, "Thread %d: %.6f s\n", i, thread_times[i]);
            }
            fclose(f);
        }
    }
}

#endif // TIMING_H

//#ifndef TIMING_H
//#define TIMING_H
//
//#include <mpi.h>
//#include <omp.h>
//
//// Global time trackers (defined in main.c)
//extern double comm_time;
//extern double comp_time;
//
//// Macro to wrap MPI calls and accumulate communication time
//#define TIME_MPI_CALL(call, time_var)      \
    //do {                                   \
        //double t0 = MPI_Wtime();           \
        //call;                              \
        //time_var += MPI_Wtime() - t0;      \
    //} while(0)
//
//// Macro to wrap compute sections and accumulate compute time
//#define TIME_OMP_BLOCK(code_block)                   \
    //do {                                             \
        //double _t0 = omp_get_wtime();                \
        //code_block                                   \
        //comp_time += omp_get_wtime() - _t0;          \
    //} while (0)
//
//#endif // TIMING_H
