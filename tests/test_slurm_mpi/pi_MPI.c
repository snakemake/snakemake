#include <stdio.h>
#include <math.h>
#include "stdlib.h"
#include "mpi.h"


double f(double x){
  return ( 4.0/(1.0 + x*x) );
}


int main(int argc,char *argv[])
{

  int myrank, size;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  long long int i, n;
  double PI25DT = 3.141592653589793238462643;
  double mypi, pi, h, mysum, x;
  
  

  /* Argument handling from command line */
  if( 2 == argc ){
    n = atoll(argv[1]);
  }else{
    if( 0 == myrank ){
      n = 10000000;
      printf("Too many or no argument given; using n = %d instead.\n", n);
    }
  }
  

  MPI_Bcast(&n, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

  h = 1.0/( (double)n );
  mysum = 0.0;
  
  for(i = myrank+1 ; i <= n ; i += size){
    x = h*((double)i - 0.5);
    mysum = mysum + f(x);
  }
  
  mypi = h*mysum;

  MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  /* Output result if myrank==0 */
  if( 0 == myrank ){
    char* filename = argv[1];   
    FILE *fp;
    char buf[0x100];    
    snprintf(buf, sizeof(buf), "%s", filename);
    fp = fopen(buf, "w+");
    fprintf(fp, "\nUsing %d processes and the value n = %d.\n",size,n);
    fprintf(fp, "Calculated pi: %.16f, with error %.16f\n\n", pi, fabs(pi - PI25DT));
    fclose(fp);    
  }

  MPI_Finalize();
  return 0;
}

