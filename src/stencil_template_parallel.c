/*

/*
 *
 *  mysizex   :   local x-extendion of your patch
 *  mysizey   :   local y-extension of your patch
 *
 */


#include "stencil_template_parallel.h"

#define HALO_TAG 42

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
    //
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
        if (Rank == 0) {
            fprintf(stderr, "Error: Need at least 2 MPI ranks\n");
        }
        MPI_Finalize();
        return 1;
    } else{
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
  double t1 = MPI_Wtime();   /* take wall-clock time */
  
  unsigned int old_frame_size, new_frame_size;
  for (int iter = 0; iter < Niterations; ++iter)
    {
      
      MPI_Request reqs[8];
      
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
      //// const int register  old_frame_size = xsize+2;
//
      //const int register width_old = planes[current].size[_x_];
      //const int register height_old = planes[current].size[_y_];
      //// const int register  new_frame_size = xsize+2;
//
      old_frame_size = (planes[current].size[_x_]+2) * (planes[current].size[_y_]+2);
      new_frame_size = (planes[!current].size[_x_]+2) * (planes[!current].size[_y_]+2);
//
      //#define IDX(j) ((j + 1) * width_old - 1) 
      //for (int i = 0; i < 4; ++i) {
	  //for (int j = 0; j < height_old; ++j){
	   //// printf("The data currently in planes[current] is: %f \n", *(planes[current].data));
//
          //(buffers)[SEND][i][j] = planes[current].data[IDX(j)];
	//
	  //} 
      //}
      //
      
      pack_halos(&planes[current], buffers, width, height, neighbours, verbose, Rank);
	
      // [B] perform the halo communications
      //     (1) use Send / Recv
	// EAST-WEST Communication
	if (neighbours[EAST] != MPI_PROC_NULL) {
		if (verbose > 0){
		    printf("Rank %d: Sending EAST to %d. First 3 values: %f, %f, %f\n", 
		   Rank, neighbours[EAST],
		   ((double*)buffers[SEND][EAST])[0],  // Cast buffer to your data type (e.g., double*)
		   ((double*)buffers[SEND][EAST])[1],
		   ((double*)buffers[SEND][EAST])[2]);
		}
	    MPI_Send(buffers[SEND][EAST], old_frame_size, MPI_BYTE, neighbours[EAST], HALO_TAG, MPI_COMM_WORLD);
	}
	if (neighbours[WEST] != MPI_PROC_NULL) {
	    MPI_Recv(buffers[RECV][WEST], new_frame_size, MPI_BYTE, neighbours[WEST], HALO_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (verbose > 0){
		    printf("Rank %d: Received EAST  from %d. First 3 values: %f, %f, %f\n", 
		   Rank, neighbours[EAST],
		   ((double*)buffers[RECV][EAST])[0],  // Cast buffer to your data type (e.g., double*)
		   ((double*)buffers[RECV][EAST])[1],
		   ((double*)buffers[RECV][EAST])[2]);
		}

	}

	if (neighbours[WEST] != MPI_PROC_NULL) {
			if (verbose > 0){
		    printf("Rank %d: Sending WEST to %d. First 3 values: %f, %f, %f\n", 
		   Rank, neighbours[WEST],
		   ((double*)buffers[SEND][WEST])[0],  // Cast buffer to your data type (e.g., double*)
		   ((double*)buffers[SEND][WEST])[1],
		   ((double*)buffers[SEND][WEST])[2]);
		}
    MPI_Send(buffers[SEND][WEST], old_frame_size, MPI_BYTE, neighbours[WEST], HALO_TAG, MPI_COMM_WORLD);
	}
	if (neighbours[EAST] != MPI_PROC_NULL) {
	    MPI_Recv(buffers[RECV][EAST], new_frame_size, MPI_BYTE, neighbours[EAST], HALO_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (verbose > 0){
		    printf("Rank %d: Received WEST  from %d. First 3 values: %f, %f, %f\n", 
		   Rank, neighbours[WEST],
		   ((double*)buffers[RECV][WEST])[0],  // Cast buffer to your data type (e.g., double*)
		   ((double*)buffers[RECV][WEST])[1],
		   ((double*)buffers[RECV][WEST])[2]);
		}

	}

	// NORTH-SOUTH Communication
	if (neighbours[NORTH] != MPI_PROC_NULL) {
		if (verbose > 0){
		    printf("Rank %d: Sending NORTH to %d. First 3 values: %f, %f, %f\n", 
		   Rank, neighbours[NORTH],
		   ((double*)buffers[SEND][NORTH])[0],  // Cast buffer to your data type (e.g., double*)
		   ((double*)buffers[SEND][NORTH])[1],
		   ((double*)buffers[SEND][NORTH])[2]);
		}

	    MPI_Send(buffers[SEND][NORTH], old_frame_size, MPI_BYTE, neighbours[NORTH], HALO_TAG, MPI_COMM_WORLD);
	}
	if (neighbours[SOUTH] != MPI_PROC_NULL) {
	    MPI_Recv(buffers[RECV][SOUTH], new_frame_size, MPI_BYTE, neighbours[SOUTH], HALO_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (verbose > 0){
		    printf("Rank %d: Received NORTH  from %d. First 3 values: %f, %f, %f\n", 
		   Rank, neighbours[NORTH],
		   ((double*)buffers[RECV][NORTH])[0],  // Cast buffer to your data type (e.g., double*)
		   ((double*)buffers[RECV][NORTH])[1],
		   ((double*)buffers[RECV][NORTH])[2]);
		}

	}

	if (neighbours[SOUTH] != MPI_PROC_NULL) {
		if (verbose > 0){
		    printf("Rank %d: Sending SOUTH to %d. First 3 values: %f, %f, %f\n", 
		   Rank, neighbours[SOUTH],
		   ((double*)buffers[SEND][SOUTH])[0],  // Cast buffer to your data type (e.g., double*)
		   ((double*)buffers[SEND][SOUTH])[1],
		   ((double*)buffers[SEND][SOUTH])[2]);
		}

	    MPI_Send(buffers[SEND][SOUTH], old_frame_size, MPI_BYTE, neighbours[SOUTH], HALO_TAG, MPI_COMM_WORLD);
	}
	if (neighbours[NORTH] != MPI_PROC_NULL) {
	    MPI_Recv(buffers[RECV][NORTH], new_frame_size, MPI_BYTE, neighbours[NORTH], HALO_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (verbose > 0){
		    printf("Rank %d: Received SOUTH  from %d. First 3 values: %f, %f, %f\n", 
		   Rank, neighbours[SOUTH],
		   ((double*)buffers[RECV][SOUTH])[0],  // Cast buffer to your data type (e.g., double*)
		   ((double*)buffers[RECV][SOUTH])[1],
		   ((double*)buffers[RECV][SOUTH])[2]);
		}
	}
      //     (2) use Isend / Irecv
      //         --> can you overlap communication and compution in this way?
      
      // [C] copy the haloes data
      update_boundaries(&planes[current], buffers, width, height, neighbours, verbose, Rank);

      /* --------------------------------------  */
      /* update grid points */
      
      update_plane( periodic, N, &planes[current], &planes[!current]);

      /* output if needed */
      if ( output_energy_stat_perstep )
	output_energy_stat ( iter, &planes[!current], (iter+1) * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
	
      /* swap plane indexes for the new iteration */
      current = !current;
      
    }
  
  t1 = MPI_Wtime() - t1;

  output_energy_stat ( -1, &planes[!current], Niterations * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
  MPI_Barrier(myCOMM_WORLD);
  memory_release( buffers, planes, Rank, verbose );
  
  
  MPI_Finalize();
  return 0;
}





















 



