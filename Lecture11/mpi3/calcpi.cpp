#include <iostream>
#include <mpi.h>
#include <math.h>

#define  N 10000000
#ifndef MARK
#define MARK fprintf(stdout,"%s @ %u\n",__FILE__,__LINE__)
#endif

using namespace std;
int main(int argc, char ** argv){
  int n, rank, size, i;
  double mypi, pi, h, sum, x;
  double startime, endtime;
     
  MPI::Init(argc, argv);
  size = MPI::COMM_WORLD.Get_size();
  rank = MPI::COMM_WORLD.Get_rank();

  cout << "rank = " << rank << endl;
  cout << "size = " << size << endl;

  while (1){
    if (rank == 0){
      cout << "enter n " << endl;
      cin >> n;
      cout << "number of intervals = " << n; 
      startime = MPI::Wtime();
    }


    MPI::COMM_WORLD.Bcast(&n, 1, MPI::INT, 0);
    if (n==0){ cout << endl; break;}
    else{
      h = 1.0 / (double) n;
      sum = 0.0;
      for (i = rank + 1; i <= n; i+= size){
	x = h* ( (double) i - 0.5);
	sum += (4.0 / (1.0 + x*x) );
      }
      mypi = h * sum;

      MPI::COMM_WORLD.Reduce(&mypi, &pi, 1, MPI::DOUBLE, MPI::SUM, 0);
	       
      if (rank == 0){
	endtime = MPI::Wtime();
	cout << "time = " << endtime - startime << endl;
	cout << "pi is around " << pi << endl;
      }
    }
  }
	  
  MPI::Finalize();
  return 0;

}


	  
     

