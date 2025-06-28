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

#define NUM_TIMED_FUNCS 2  // 0 = get_total_energy, 1 = update_plane

extern double thread_times[NUM_TIMED_FUNCS][MAX_THREADS]; // Thread-local timings


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
    int num_threads = omp_get_max_threads();

	if (num_threads > 1 && Rank == 0) {
		printf("--- Thread-level imbalance (per rank) ---\n");

		for (int func_id = 0; func_id < NUM_TIMED_FUNCS; ++func_id) {
		    double max_time = 0.0, min_time = 1e9;
		    for (int t = 0; t < num_threads; ++t) {
			double tval = thread_times[func_id][t];
			if (tval > max_time) max_time = tval;
			if (tval < min_time) min_time = tval;
		    }

		    printf("Func %d | Max: %.6f s, Min: %.6f s, Ratio: %.2f%%\n",
			   func_id, max_time, min_time,
			   100.0 * (max_time - min_time) / (max_time > 0.0 ? max_time : 1.0));
		}
	    }

    // Optional per-rank log file
    if (log_per_rank) {
        char fname[64];
        snprintf(fname, sizeof(fname), "timing_rank_%d.log", Rank);
        FILE *f = fopen(fname, "w");
        if (f) {
            fprintf(f, "RANK %d\n", Rank);
            fprintf(f, "Label: %s\n", label);
            fprintf(f, "Total time: %.6f s\n", total_time);
            fprintf(f, "Comm time:  %.6f s\n", comm_time);
            fprintf(f, "Comp time:  %.6f s\n", compute_time);
            fprintf(f, "Mem usage:  %.2f MB\n", mem_MB);
            fprintf(f, "Threads:    %d\n\n", num_threads);

            for (int func_id = 0; func_id < NUM_TIMED_FUNCS; ++func_id) {
                fprintf(f, "Thread timings for function %d:\n", func_id);
                for (int t = 0; t < num_threads; ++t) {
                    fprintf(f, "  Thread %2d: %.6f s\n", t, thread_times[func_id][t]);
                }
                fprintf(f, "\n");
            }

            fclose(f);
        } else {
            fprintf(stderr, "Rank %d: Failed to open %s for writing.\n", Rank, fname);
        }
    }
}

#endif // TIMING_H

