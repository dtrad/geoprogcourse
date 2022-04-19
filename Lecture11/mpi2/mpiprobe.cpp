#include <iostream>
#include <mpi.h>
using namespace std;
#define MARK fprintf(stdout,"%s @ %u\n",__FILE__,__LINE__)
int main(int argc, char ** argv){
  int size;
  int rank;
  float x=0;
  int i=0;     
  MPI::Intracomm world;
  MPI::Group world_group;
  MPI::Status status;     
  MPI::Init(argc,argv);
  world = MPI::COMM_WORLD;
  size  = world.Get_size();
  rank  = world.Get_rank();
  if (rank == 0){
    i = 10;
    cout << "rank=" << rank << ": size=" << size << endl;
    world.Send(&i,1,MPI::INT,2,0);
  }
  else if (rank == 1){
    x = 1.1;
    cout << "rank=" << rank << endl;
    world.Send(&x,1,MPI::FLOAT,2,0);
  }
  else if (rank == 2){
    for (int ix=1;ix<=2;ix++){
      world.Probe(MPI::ANY_SOURCE,0,status);
      if (status.Get_source() == 0){
	world.Recv(&i,1,MPI::INT,0,0);
      }
      else if (status.Get_source() == 1){
	world.Recv(&x,1,MPI::FLOAT,1,0);
      }
    }
  }     
  if (rank == 2){
    cout << "rank = " << rank << endl;
    cout << " x = " << x << endl;
    cout << " i = " << i << endl;
  }  
  MPI::Finalize();
  return 0;

}


	  
     

