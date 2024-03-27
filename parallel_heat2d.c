# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <omp.h>
# include <string.h>
int main ( int argc, char *argv[] );
  double u[10000][10000];
  double w[10000][10000];
/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for HEATED_PLATE_OPENMP.
    
  Parameters:

    Commandline argument 1, double EPSILON, the error tolerance.  

    Commandline argument 2, char *OUTPUT_FILE, the name of the file into which
    the steady state solution is written when the program has completed.

    Commandline argument 3, Integer N , the matrix NxN
  Local parameters:

    Local, double DIFF, the norm of the change in the solution from one iteration
    to the next.

    Local, double MEAN, the average of the boundary values, used to initialize
    the values of the solution in the interior.

    Local, double U[M][N], the solution at the previous iteration.

    Local, double W[M][N], the solution computed at the latest iteration.
    OUTPUT:
    output.txt          :Time is appended to it   
    parallel_{N}.gp     :used for gnuplot
    parallel_{N}.txt    :used to store final value
*/
{
# define M atoi(argv[3])
# define N atoi(argv[3])

  double diff;
  double epsilon = 0.001;
  sscanf ( argv[1], "%lf", &epsilon );
  int i;
  int iterations;
  int iterations_print;
  int j;
  FILE *fp;
  double mean;
  char output_file[80];
  sscanf ( argv[2], "%s", output_file );
  double my_diff;

  double wtime;

  printf ( "\n" );
  printf ( "HEATED_PLATE_OPENMP\n" );
  printf ( "  C/OpenMP version\n" );
  printf ( "  A program to solve for the steady state temperature distribution\n" );
  printf ( "  over a rectangular plate.\n" );
  printf ( "\n" );
  printf ( "\n  The steady state solution will be written to '%s'.\n", output_file );
  printf ( "  Spatial grid of %d by %d points.\n", M, N );
  printf ( "  The iteration will be repeated until the change is <= %e\n", epsilon ); 
  printf ( "  Number of processors available = %d\n", omp_get_num_procs ( ) );
  printf ( "  Number of threads =              %d\n", omp_get_max_threads ( ) );
/*
  Set the boundary values, which don't change. 
*/
  mean = 0.0;

#pragma omp parallel shared ( w ) private ( i, j )
  {
#pragma omp for
    for ( i = 1; i < M - 1; i++ )
    {
      w[i][0] = 100.0;
    }
#pragma omp for
    for ( i = 1; i < M - 1; i++ )
    {
      w[i][N-1] = 100.0;
    }
#pragma omp for
    for ( j = 0; j < N; j++ )
    {
      w[M-1][j] = 100.0;
    }
#pragma omp for
    for ( j = 0; j < N; j++ )
    {
      w[0][j] = 0.0;
    }
/*
  Average the boundary values, to come up with a reasonable
  initial value for the interior.
*/
#pragma omp for reduction ( + : mean )
    for ( i = 1; i < M - 1; i++ )
    {
      mean = mean + w[i][0] + w[i][N-1];
    }
#pragma omp for reduction ( + : mean )
    for ( j = 0; j < N; j++ )
    {
      mean = mean + w[M-1][j] + w[0][j];
    }
  }
/*
  OpenMP note:
  You cannot normalize MEAN inside the parallel region.  It
  only gets its correct value once you leave the parallel region.
  So we interrupt the parallel region, set MEAN, and go back in.
*/
  mean = mean / ( double ) ( 2 * M + 2 * N - 4 );
  printf ( "\n" );
  printf ( "  MEAN = %f\n", mean );
/* 
  Initialize the interior solution to the mean value.
*/
#pragma omp parallel shared ( mean, w ) private ( i, j )
  {
#pragma omp for
    for ( i = 1; i < M - 1; i++ )
    {
      for ( j = 1; j < N - 1; j++ )
      {
        w[i][j] = mean;
      }
    }
  }
/*
  iterate until the  new solution W differs from the old solution U
  by no more than EPSILON.
*/
  iterations = 0;
  iterations_print = 1;
  printf ( "\n" );
  printf ( " Iteration  Change\n" );
  printf ( "\n" );
  wtime = omp_get_wtime ( );

  diff = epsilon;

  while ( epsilon <= diff )
  {
# pragma omp parallel shared ( u, w ) private ( i, j )
    {
/*
  Save the old solution in U.
*/
# pragma omp for
      for ( i = 0; i < M; i++ ) 
      {
        for ( j = 0; j < N; j++ )
        {
          u[i][j] = w[i][j];
        }
      }
/*
  Determine the new estimate of the solution at the interior points.
  The new solution W is the average of north, south, east and west neighbors.
*/
# pragma omp for
      for ( i = 1; i < M - 1; i++ )
      {
        for ( j = 1; j < N - 1; j++ )
        {
          w[i][j] = ( u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1] ) / 4.0;
        }
      }
    }
/*
  C and C++ cannot compute a maximum as a reduction operation.

  Therefore, we define a private variable MY_DIFF for each thread.
  Once they have all computed their values, we use a CRITICAL section
  to update DIFF.
*/
    diff = 0.0;
# pragma omp parallel shared ( diff, u, w ) private ( i, j, my_diff )
    {
      my_diff = 0.0;
# pragma omp for
      for ( i = 1; i < M - 1; i++ )
      {
        for ( j = 1; j < N - 1; j++ )
        {
          if ( my_diff < fabs ( w[i][j] - u[i][j] ) )
          {
            my_diff = fabs ( w[i][j] - u[i][j] );
          }
        }
      }
# pragma omp critical
      {
        if ( diff < my_diff )
        {
          diff = my_diff;
        }
      }
    }

    iterations++;
    if ( iterations == iterations_print )
    {
      printf ( "  %8d  %f\n", iterations, diff );
      iterations_print = 2 * iterations_print;
    }
  } 
  wtime = omp_get_wtime ( ) - wtime;

  printf ( "\n" );
  printf ( "%dX%d PARALLEL %8d  %f\n",N,N, iterations, diff );
  printf ( "\n" );
  // printf ( "  Error tolerance achieved.\n" );
  printf ( "%dX%d PARALLEL Wallclock time = %f\n",N,N, wtime );
  FILE *out_file = fopen("output.txt","a");
  fprintf ( out_file,"MATRIX: %d X %d\t PARALLEL CPU time = %f\n",N,N, wtime );
  fclose(out_file);
  printf("===========================================================\n");
  /* 
  Write the solution to the output file.
  */
 
    fp = fopen ( output_file, "w" );

    // fprintf ( fp, "%d\n", M );
    // fprintf ( fp, "%d\n", N );

    for ( i = 0; i < M; i++ )
    {
      for ( j = 0; j < N; j++)
      {
        fprintf ( fp, "%6.2f ", w[i][j] );
      }
      fputc ( '\n', fp);
    }
    fclose ( fp );

    printf ( "\n" );
    printf ("  Solution written to the output file '%s'\n", output_file );
    // Image Visualization
    char plt_name[30];
    snprintf(plt_name, sizeof(plt_name), "parallel_%d.gp", N);
    FILE* plt_file = fopen(plt_name,"w");
    fprintf(plt_file,"set pm3d map\n");
    fprintf(plt_file,"splot \'%s\' matrix with image\n",output_file);
    fclose(plt_file);
  /* 
    All done! 
  */
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "HEATED_PLATE_OPENMP:\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;

# undef M
# undef N
}