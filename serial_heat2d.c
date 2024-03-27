# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

int main ( int argc, char *argv[] );
double cpu_time ( );
  double u[10000][10000];
  double w[10000][10000];
/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*

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
    serial_{N}.gp       :used for gnuplot
    serial_{N}.txt      :used to store final value
*/
{

  int M=atoi(argv[3]);
  int N=atoi(argv[3]);
    // printf("%d %d\n",atoi(argv[3]),atoi(argv[4]));
  double ctime;
  double ctime1;
  double ctime2;
  double diff;
  double epsilon;
  FILE *fp;
  int i;
  int iterations;
  int iterations_print;
  int j;
  double mean;
  char output_file[80];
  sscanf ( argv[2], "%s", output_file );
  int success;


  printf ( "\n" );
  printf ( "HEATED_PLATE\n" );
  printf ( "  C version\n" );
  printf ( "  A program to solve for the steady state temperature distribution\n" );
  printf ( "  over a rectangular plate.\n" );
  printf ( "\n" );
  printf ( "  Spatial grid of %d by %d points.\n", M, N );
/* 
  Read EPSILON from the command line or the user.
*/
  if ( argc < 2 ) 
  {
    // printf ( "\n" );
    // printf ( "  Enter EPSILON, the error tolerance:\n" );
    success = scanf ( "%lf", &epsilon );
  }
  else
  {
    success = sscanf ( argv[1], "%lf", &epsilon );
  }

  if ( success != 1 )
  {
    // printf ( "\n" );
    // printf ( "HEATED_PLATE\n" );
    // printf ( "  Error reading in the value of EPSILON.\n");
    return 1;
  }

  // printf ( "\n" );
  // printf ( "  The iteration will be repeated until the change is <= %f\n", epsilon );
  diff = epsilon;
/* 
  Read OUTPUT_FILE from the command line or the user.
*/
/*
  if ( argc < 3 ) 
  {
    printf ( "\n" );
    printf ( "  Enter OUTPUT_FILE, the name of the output file:\n" );
    success = scanf ( "%s", output_file );
  }
  else
  {
    success = sscanf ( argv[2], "%s", output_file );
  }

  if ( success != 1 )
  {
    printf ( "\n" );
    printf ( "HEATED_PLATE\n" );
    printf ( "  Error reading in the value of OUTPUT_FILE.\n");
    return 1;
  }

  printf ( "\n" );
  printf ( "  The steady state solution will be written to '%s'.\n", output_file );
*/

/* 
  Set the boundary values, which don't change. 
*/
  for ( i = 1; i < M - 1; i++ )
  {
    w[i][0] = 100.0;
  }
  for ( i = 1; i < M - 1; i++ )
  {
    w[i][N-1] = 100.0;
  }
  for ( j = 0; j < N; j++ )
  {
    w[M-1][j] = 100.0;
  }
  for ( j = 0; j < N; j++ )
  {
    w[0][j] = 0.0;
  }
/*
  Average the boundary values, to come up with a reasonable
  initial value for the interior.
*/
  mean = 0.0;
  for ( i = 1; i < M - 1; i++ )
  {
    mean = mean + w[i][0];
  }
  for ( i = 1; i < M - 1; i++ )
  {
    mean = mean + w[i][N-1];
  }
  for ( j = 0; j < N; j++ )
  {
    mean = mean + w[M-1][j];
  }
  for ( j = 0; j < N; j++ )
  {
    mean = mean + w[0][j];
  }
  mean = mean / ( double ) ( 2 * M + 2 * N - 4 );
/* 
  Initialize the interior solution to the mean value.
*/
  for ( i = 1; i < M - 1; i++ )
  {
    for ( j = 1; j < N - 1; j++ )
    {
      w[i][j] = mean;
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
  ctime1 = cpu_time ( );

  while ( epsilon <= diff )
  {
/*
  Save the old solution in U.
*/
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
    diff = 0.0;
    for ( i = 1; i < M - 1; i++ )
    {
      for ( j = 1; j < N - 1; j++ )
      {
        w[i][j] = ( u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1] ) / 4.0;

        if ( diff < fabs ( w[i][j] - u[i][j] ) )
        {
          diff = fabs ( w[i][j] - u[i][j] );
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
  ctime2 = cpu_time ( );
  ctime = ctime2 - ctime1;

  printf ( "\n" );
  printf ( "%dX%d SERIAL  %8d  %f\n",N,N, iterations, diff );
  printf ( "\n" );
  printf ( "  Error tolerance achieved.\n" );
  printf ( "%dX%d SERIAL CPU time = %f\n",N,N, ctime );
  FILE *out_file = fopen("output.txt","a");
  fprintf ( out_file,"MATRIX: %d X %d\t SERIAL CPU time = %f\n",N,N, ctime );
  fclose(out_file);
/* 
  Write the solution to the output file.
*/

  fp = fopen ( output_file, "w" );

//   fprintf ( fp, "%d\n", M );
//   fprintf ( fp, "%d\n", N );

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
  /*IMAGE VISUALIZATION
  */
  printf ("  Solution written to the output file '%s'\n", output_file );
  char plt_name[30];
  snprintf(plt_name, sizeof(plt_name), "serial_%d.gp", N);
  FILE* plt_file = fopen(plt_name,"w");
  fprintf(plt_file,"set pm3d map\n");
  fprintf(plt_file,"splot \'%s\' matrix with image\n",output_file);
  fclose(plt_file);
  
/* 
  All done! 
*/
  // printf ( "\n" );
  // printf ( "HEATED_PLATE:\n" );
  // printf ( "  Normal end of execution.\n" );
  printf("=====================================================================\n");

  return 0;

# undef M
# undef N
}
/******************************************************************************/

double cpu_time ( void )

/******************************************************************************/
/*
  Purpose:

    CPU_TIME returns the current reading on the CPU clock.
  Parameters:

    Output, double CPU_TIME, the current reading of the CPU clock, in seconds.
*/
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}