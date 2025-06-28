/*

/*
 *
 *  mysizex   :   local x-extendion of your patch
 *  mysizey   :   local y-extension of your patch
 *
 */


#include "stencil_template_parallel.h"


// global time counters
double comm_time = 0.0, comp_time = 0.0;
double loop_start_time = 0.0, loop_end_time = 0.0;
double thread_times[NUM_TIMED_FUNCS][MAX_THREADS] = {1.0};
// ------------------------------------------------------------------
// ------------------------------------------------------------------

int main(int argc, char **argv)
{

  MPI_Comm myCOMM_WORLD;
  int  Rank, Ntasks;
  int verbose = 0;
  uint neighbours[4];
  int  Niterations;
  int  periodic;
  vec2_t S, N;
  
  int      Nsources;
  int      Nsources_local;
  vec2_t  *Sources_local;
  double   energy_per_source;

  plane_t   planes[2];  
  buffers_t buffers[2];
  
  int output_energy_stat_perstep;
  
  /* initialize MPI envrionment */
  {
    int level_obtained;
    
    // NOTE: change MPI_FUNNELED if appropriate
    MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &level_obtained );
    if ( level_obtained < MPI_THREAD_FUNNELED ) {
      printf("MPI_thread level obtained is %d instead of %d\n",
	     level_obtained, MPI_THREAD_FUNNELED );
      MPI_Finalize();
      exit(1); }
    
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
    MPI_Comm_size(MPI_COMM_WORLD, &Ntasks);
    MPI_Comm_dup (MPI_COMM_WORLD, &myCOMM_WORLD);

  

  }
  
  /* argument checking and setting */
  int ret = initialize ( &myCOMM_WORLD, Rank, Ntasks, argc, argv, &S, &N, &periodic, &output_energy_stat_perstep,
			 neighbours, &Niterations,
			 &Nsources, &Nsources_local, &Sources_local, &energy_per_source,
			 &planes[0], &buffers[0], &verbose);
  // Validate minimum ranks
  if (verbose > 0){
	  printf("Hello from rank %d of %d\n", Rank, Ntasks);
	  fflush(stdout);
  }

  if (Ntasks < 2) {
        //if (Rank == 0) {
            //fprintf(stderr, "Error: Need at least 2 MPI ranks\n");
        //}
        //MPI_Finalize();
        //return 1;
	  if (Rank == 0) {
	    fprintf(stderr, "Warning: Running in serial mode with a single MPI rank. Communications will be skipped.\n");
	  }
	  // set all neighbours to MPI_PROC_NULL to avoid sends/recvs
	  for (int i = 0; i < 4; ++i) neighbours[i] = MPI_PROC_NULL;	
  }
  else{
        printf("We have %d tasks to deal with \n", Ntasks);
  }

  if ( ret )
    {
      printf("task %d is opting out with termination code %d\n",
	     Rank, ret );
      fflush(stdout);
      
      MPI_Finalize();
      return 0;
    }
  else{
	if (verbose > 0){
	      printf("task %d is starting successfully with code %d\n",
		     Rank, ret );
	      fflush(stdout);
	
	}
  }
  
  
  int current = OLD;
  START_LOOP_TIMER(); 
  for (int iter = 0; iter < Niterations; ++iter)
    {
      
      /* new energy from sources */
      int inj_ret = inject_energy( periodic, Nsources_local, Sources_local, energy_per_source, N, &planes[current], verbose);
      if (inj_ret != 0) {
	    fprintf(stderr, "Rank %d: inject_energy failed with code %d\n", Rank, inj_ret);
	    MPI_Abort(MPI_COMM_WORLD, inj_ret);
      }
      else {
	      if (verbose > 0){
		    printf("Rank %d: inject_energy succeeded\n", Rank);
		    fflush(stdout);
	      }
      }

      /* -------------------------------------- */

      // [A] fill the buffers, and/or make the buffers' pointers pointing to the correct position

      const int register width = planes[current].size[_x_];
      const int register height = planes[current].size[_y_];
      const int register buffer_width = planes[current].size[_x_] + 2;
      const int register buffer_height = planes[current].size[_y_] + 2;
        
      pack_halos(&planes[current], buffers, width, height, neighbours, verbose, Rank);
	
      /* // [B] perform the halo communications */

      send_halos(buffers, neighbours, buffer_width, buffer_height, Rank, verbose, 1);      
      // [C] copy the haloes data

      update_halos(&planes[current], buffers, width, height, neighbours, verbose, Rank);

      /* --------------------------------------  */
      /* update grid points */
      
      update_plane( periodic, N, &planes[current], &planes[!current]);

      /* output if needed */
      if ( output_energy_stat_perstep )
	output_energy_stat ( iter, &planes[!current], (iter+1) * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
	
      /* swap plane indexes for the new iteration */
      current = !current;
      
    }
  STOP_LOOP_TIMER();  

  output_energy_stat ( -1, &planes[!current], Niterations * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
  if (Ntasks > 1){
  	TIME_MPI_CALL(MPI_Barrier(myCOMM_WORLD), comm_time);
  }

  memory_release(buffers, planes, Rank, verbose);

  report_timing_stats(myCOMM_WORLD, Rank, Ntasks, (Ntasks == 1) ? "Serial Mode" : "Main Loop", 0);
  
  MPI_Finalize();
  return 0;
}





















 



