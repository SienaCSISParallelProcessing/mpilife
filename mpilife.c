/* Conway's Game of Life - MPI implementation */
/* Jim Teresco, Tue Mar 18 22:42:46 EST 1997 */
/* Updated 20 years later for CSIS 400, Siena College, Fall 2017 */

/* 
   Preprocessor defines that will enable extra outputs:

   OUTPUT_ALL: print out the world on every iteration
   OUTPUT_END: print out the world at the end
   USE_FILE: used in combination with the above, if defined, print to a file
     (or files, when OUTPUT_ALL) with basename as defined in USE_FILE
*/

#define OUTPUT_END
//#define USE_FILE "lifeworld"

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>

int main(int argc, char *argv[]) {
  int gridsize, myrows;
  double init_pct;
  int num_iters, iter;
  int **grid[2], curr, prev;
  int i, j, neigh_count;
  long live_count, birth_count, death_count;
  long global_live, global_birth, global_death;
  int numprocs, mypid;
  int msgcount;
  MPI_Request requests[4];
  MPI_Status status[4];
  FILE *fp;
  int proc_turn;
#ifdef USE_FILE
  char filename[FILENAME_MAX];
#endif
  
  MPI_Init(&argc,&argv);

  if (argc != 4) {
    fprintf(stderr,"Usage: %s gridsize init_pct num_iters\n",argv[0]);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  MPI_Comm_rank(MPI_COMM_WORLD,&mypid);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

  srand48(time(NULL)*(mypid+1));

  /* read parameters from the command line */
  gridsize=atoi(argv[1]);
  init_pct=atof(argv[2]);
  num_iters=atoi(argv[3]);

  if (gridsize%numprocs) {
    fprintf(stderr,"%s: grid size must be a multiple of number of procs\n", argv[0]);
    MPI_Abort(MPI_COMM_WORLD,1);
  }

  myrows=gridsize/numprocs;
  if (mypid == 0)
    printf("Using grid size %d (%d rows on each of %d procs)\n",
	   gridsize,myrows,numprocs);

  /* allocate the grids */
  grid[0]=(int **)malloc((myrows+2)*sizeof(int *));
  for (i=0;i<=myrows+1;i++)
    grid[0][i]=(int *)malloc((gridsize+2)*sizeof(int));
  grid[1]=(int **)malloc((myrows+2)*sizeof(int *));
  for (i=0;i<=myrows+1;i++)
    grid[1][i]=(int *)malloc((gridsize+2)*sizeof(int));

  /* initialize the grids (incl boundary buffer all 0's) */
  for (i=0;i<=myrows+1;i++)
    for (j=0;j<=gridsize+1;j++) {
      grid[0][i][j]=0;
      grid[1][i][j]=0;
    }

  /* start current grid as 0 */
  curr=0; prev=1;

  /* initialize the current grid based on the desired percentage of
     living cells specified on the command line */
  live_count=0;
  for (i=1;i<=myrows;i++)
    for (j=1;j<=gridsize;j++) {
      if (drand48()<init_pct) {
	grid[curr][i][j]=1;
	live_count++;
      }
      else grid[curr][i][j]=0;
    }

  printf("Proc %d: Initial grid has %ld live cells out of %ld\n",mypid,
	 live_count,(long)myrows*gridsize);
  MPI_Reduce(&live_count, &global_live, 1, MPI_LONG, MPI_SUM, 0,
	     MPI_COMM_WORLD);
  if (mypid == 0)
    printf("Global:  Initial grid has %ld live cells out of %ld\n",global_live,
	   (long)gridsize*gridsize);

#ifdef OUTPUT_ALL
  for (proc_turn = 0; proc_turn < numprocs; proc_turn++) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (mypid == proc_turn) {
#ifdef USE_FILE
      sprintf(filename, "%s.000000.txt", USE_FILE);
      fp = fopen(filename, "a");
#else
      fp = stdout;
#endif
      if (mypid == 0) {
	fprintf(fp, "Initial grid:\n");
      }
      for (i=1; i<=myrows; i++) {
	fprintf(fp, "[%3d] ", mypid);
	for (j=1; j<=gridsize; j++) {
	  fprintf(fp, "%c", (grid[curr][i][j] ? '*' : '-'));
	}
	fprintf(fp, "\n");
      }
#ifdef USE_FILE
      fclose(fp);
#else
      fflush(stdout);
#endif
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  /* we can now start iterating */
  for (iter=1; iter<=num_iters; iter++) {

    /* swap the grids */
    curr=1-curr; prev=1-prev;

    if (mypid == 0)
      printf("Iteration %d...\n",iter);

    /* update the local buffers of off-processor neighbors in prev grid */
    /* send from proc to next highest and next lowest neighbor */
    msgcount=0;
    if (mypid < numprocs-1) {
      MPI_Isend(grid[prev][myrows],gridsize+2,MPI_INT,mypid+1,1,
		MPI_COMM_WORLD,&requests[msgcount]);
      msgcount++;

      MPI_Irecv(grid[prev][myrows+1],gridsize+2,MPI_INT,mypid+1,2,
		MPI_COMM_WORLD,&requests[msgcount]);
      msgcount++;
    }

    if (mypid > 0) {
      MPI_Isend(grid[prev][1],gridsize+2,MPI_INT,mypid-1,2,
		MPI_COMM_WORLD,&requests[msgcount]);
      msgcount++;

      MPI_Irecv(grid[prev][0],gridsize+2,MPI_INT,mypid-1,1,
		MPI_COMM_WORLD,&requests[msgcount]);
      msgcount++;
    }

    /* wait for communication to complete */
    MPI_Waitall(msgcount,requests,status);

    /* perform the actual iteration */
    live_count=0; birth_count=0; death_count=0;

    /* visit each grid cell */
    for (i=1;i<=myrows;i++)
      for (j=1;j<=gridsize;j++) {
	neigh_count=
	  (grid[prev][i-1][j-1]+grid[prev][i-1][j]+grid[prev][i-1][j+1]+
	   grid[prev][i][j-1]+grid[prev][i][j+1]+
	   grid[prev][i+1][j-1]+grid[prev][i+1][j]+grid[prev][i+1][j+1]);
	switch (neigh_count) {
	case 2:
	  /* no change */
	  grid[curr][i][j]=grid[prev][i][j];
	  break;
	case 3:
	  /* birth */
	  if (!grid[prev][i][j]) birth_count++;
	  grid[curr][i][j]=1;
	  break;
	default:
	  /* death of loneliness or overcrowding */
	  if (grid[prev][i][j]) death_count++;
	  grid[curr][i][j]=0;
	  break;
	}
	live_count+=grid[curr][i][j];
      }

    /* print the stats */
    printf("Proc %d Counters- living: %ld, died: %ld, born: %ld\n",mypid,
	   live_count, death_count, birth_count);
    MPI_Reduce(&live_count, &global_live, 1, MPI_LONG, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&death_count, &global_death, 1, MPI_LONG, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&birth_count, &global_birth, 1, MPI_LONG, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    if (mypid == 0)
      printf("Global Counters- living: %ld, died: %ld, born: %ld\n",
	   global_live, global_death, global_birth);

#ifdef OUTPUT_ALL
    for (proc_turn = 0; proc_turn < numprocs; proc_turn++) {
      MPI_Barrier(MPI_COMM_WORLD);
      if (mypid == proc_turn) {
#ifdef USE_FILE
	sprintf(filename, "%s.%06d.txt", USE_FILE, iter);
	fp = fopen(filename, "a");
#else
	fp = stdout;
#endif
	if (mypid == 0) {
	  fprintf(fp, "Grid at iter %d:\n", iter);
	}
	for (i=1; i<=myrows; i++) {
	  fprintf(fp, "[%3d] ", mypid);
	  for (j=1; j<=gridsize; j++) {
	    fprintf(fp, "%c", (grid[curr][i][j] ? '*' : '-'));
	  }
	  fprintf(fp, "\n");
	}
#ifdef USE_FILE
	fclose(fp);
#else
	fflush(stdout);
#endif
      }
    }
#endif
    
    /* end of iteration - synchronize */
    /* this is not strictly necessary, as the reduce collective
       communication operations will synchronize the processes */
    MPI_Barrier(MPI_COMM_WORLD);
  }

  
#ifdef OUTPUT_END
  for (proc_turn = 0; proc_turn < numprocs; proc_turn++) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (mypid == proc_turn) {
#ifdef USE_FILE
      sprintf(filename, "%s.final.txt", USE_FILE);
      fp = fopen(filename, "a");
#else
      fp = stdout;
#endif
      if (mypid == 0) {
	fprintf(fp, "Final grid:\n");
      }
      for (i=1; i<=myrows; i++) {
	fprintf(fp, "[%3d] ", mypid);
	for (j=1; j<=gridsize; j++) {
	  fprintf(fp, "%c", (grid[curr][i][j] ? '*' : '-'));
	}
	fprintf(fp, "\n");
      }
#ifdef USE_FILE
      fclose(fp);
#else
      fflush(stdout);
#endif
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  /* free the grids */
  for (i=0;i<=myrows+1;i++)
    free(grid[0][i]);
  free(grid[0]);
  for (i=0;i<=myrows+1;i++)
    free(grid[1][i]);
  free(grid[1]);

  MPI_Finalize();
  return 0;
}
