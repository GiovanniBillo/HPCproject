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

    // Validate minimum ranks
  printf("Hello from rank %d of %d\n", Rank, Ntasks);
  fflush(stdout);
    if (Ntasks < 2) {
        if (Rank == 0) {
            fprintf(stderr, "Error: Need at least 2 MPI ranks\n");
        }
        MPI_Finalize();
        return 1;
    } else{
        printf("We have %d tasks to deal with \n", Ntasks);
    }

  }
  
  
  /* argument checking and setting */
  int ret = initialize ( &myCOMM_WORLD, Rank, Ntasks, argc, argv, &S, &N, &periodic, &output_energy_stat_perstep,
			 neighbours, &Niterations,
			 &Nsources, &Nsources_local, &Sources_local, &energy_per_source,
			 &planes[0], &buffers[0] );

  if ( ret )
    {
      printf("task %d is opting out with termination code %d\n",
	     Rank, ret );
      fflush(stdout);
      
      MPI_Finalize();
      return 0;
    }
  else{
      printf("task %d is starting successfully with code %d\n",
	     Rank, ret );
      fflush(stdout);
  }
  
  
  int current = OLD;
  double t1 = MPI_Wtime();   /* take wall-clock time */
  
  unsigned int old_frame_size, new_frame_size;
  for (int iter = 0; iter < Niterations; ++iter)
    {
      
      MPI_Request reqs[8];
      
      /* new energy from sources */
      int inj_ret = inject_energy( periodic, Nsources_local, Sources_local, energy_per_source, N, &planes[current]);
      if (inj_ret != 0) {
	    fprintf(stderr, "Rank %d: inject_energy failed with code %d\n", Rank, inj_ret);
	    MPI_Abort(MPI_COMM_WORLD, inj_ret);
      }
      else {
	    printf("Rank %d: inject_energy succeeded\n", Rank);
	    fflush(stdout);
      }

      /* -------------------------------------- */

      // [A] fill the buffers, and/or make the buffers' pointers pointing to the correct position

      old_frame_size = (planes[current].size[_x_]+2) * (planes[current].size[_y_]+2);

      new_frame_size = (planes[!current].size[_x_]+2) * (planes[!current].size[_y_]+2);
      for (int i = 0; i < 4; ++i) {
	  // printf("The data currently in planes[current] is: %f \n", *(planes[current].data));

          (buffers)[SEND][i] = planes[current].data;
	
          (buffers)[RECV][i] = planes[!current].data;
	  
      // [B] perform the halo communications
      //     (1) use Send / Recv
    if (neighbours[EAST] != MPI_PROC_NULL){

	MPI_Send(
	    (buffers)[SEND][EAST],
	    old_frame_size,
	    MPI_BYTE, 
	    neighbours[EAST],
	    HALO_TAG,
	    MPI_COMM_WORLD);

	MPI_Recv(
	    (buffers)[RECV][EAST],
	    new_frame_size,
	    MPI_BYTE, 
	    neighbours[EAST],
	    HALO_TAG,
	    MPI_COMM_WORLD, 
	    MPI_STATUS_IGNORE);
    }
    else if (neighbours[WEST] != MPI_PROC_NULL){
	MPI_Send(
	    (buffers)[SEND][WEST],
	    old_frame_size,
	    MPI_BYTE, 
	    neighbours[WEST],
	    HALO_TAG,
	    MPI_COMM_WORLD);

	MPI_Recv(
	    (buffers)[RECV][WEST],
	    new_frame_size,
	    MPI_BYTE, 
	    neighbours[WEST],
	    HALO_TAG,
	    MPI_COMM_WORLD, 
	    MPI_STATUS_IGNORE);

    }else if (neighbours[NORTH] != MPI_PROC_NULL){
	MPI_Send(
	    (buffers)[SEND][NORTH],
	    old_frame_size,
	    MPI_BYTE, 
	    neighbours[NORTH],
	    HALO_TAG,
	    MPI_COMM_WORLD);

	MPI_Recv(
	    (buffers)[RECV][NORTH],
	    new_frame_size,
	    MPI_BYTE, 
	    neighbours[NORTH],
	    HALO_TAG,
	    MPI_COMM_WORLD, 
	    MPI_STATUS_IGNORE);

    }
    else if (neighbours[SOUTH] != MPI_PROC_NULL){
	MPI_Send(
	    (buffers)[SEND][SOUTH],
	    old_frame_size,
	    MPI_BYTE, 
	    neighbours[SOUTH],
	    HALO_TAG,
	    MPI_COMM_WORLD);

	MPI_Recv(
	    (buffers)[RECV][SOUTH],
	    new_frame_size,
	    MPI_BYTE, 
	    neighbours[SOUTH],
	    HALO_TAG,
	    MPI_COMM_WORLD, 
	    MPI_STATUS_IGNORE);
    }

      //     (2) use Isend / Irecv
      //         --> can you overlap communication and compution in this way?
      
      // [C] copy the haloes data
    }        
      /* --------------------------------------  */
      /* update grid points */
      
      update_plane( periodic, N, &planes[current], &planes[!current] );

      /* output if needed */
      if ( output_energy_stat_perstep )
	output_energy_stat ( iter, &planes[!current], (iter+1) * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
	
      /* swap plane indexes for the new iteration */
      current = !current;
      
    }
  
  t1 = MPI_Wtime() - t1;

  output_energy_stat ( -1, &planes[!current], Niterations * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
  MPI_Barrier(myCOMM_WORLD);
  memory_release( buffers, planes );
  
  
  MPI_Finalize();
  return 0;
}





















 



