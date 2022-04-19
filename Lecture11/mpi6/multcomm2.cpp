#include <iostream>
#include <mpi.h>
#include <math.h>
#include "MyAlloc.hpp"
#include <stdlib.h>
#include <limits.h>
//******************************************************************************
// mpi example 6: Geophysical programming University of Calgary
// Create communicator example 
// credits: Daniel Trad (based on Using MPI by Gropp, Lusk, Skjellum)
//*****************************************************************************/

#define MARK fprintf(stdout,"%s @ %u\n",__FILE__,__LINE__)
#define CHUNKSIZE 1000
#define REQUEST 1
#define REPLY   2

using namespace std;
int main(int argc, char ** argv){
  int iter;
  int in, out, i, max, ranks[1], done;
  double x, y, Pi, error, epsilon;
  int numprocs, myid, server, totalin, totalout, workerid;
  int rands[CHUNKSIZE], request;

  MPI::Intracomm world, workers;
  MPI::Group world_group, worker_group;
  MPI::Status status;

  MPI::Init(argc,argv);
  world = MPI::COMM_WORLD;

  numprocs = world.Get_size();
  myid    = world.Get_rank();
     
  server = numprocs - 1;
  if (myid == 0) epsilon = 0.001;
  
  world.Bcast(&epsilon,1,MPI::DOUBLE,0);
  world_group = world.Get_group();
  ranks[0]=server;
  worker_group = world_group.Excl(1, ranks);
  workers = world.Create(worker_group);
     
  if (myid == server){
    cout << "server node = " << server << endl;
    do {
      world.Recv(&request,1,MPI::INT,MPI::ANY_SOURCE,REQUEST,status);
      if (request){
	cout << "request from " << status.Get_source() << endl;
	for (i=0; i< CHUNKSIZE;i++) rands[i]=random();
	world.Send(rands,CHUNKSIZE,MPI_INT,status.Get_source(),REPLY);

      }
    }while (request > 0);
  }
  else{
    request = 1;
    done = in = out = 0;
    max = INT_MAX;
    world.Send(&request,1,MPI::INT,server,REQUEST);
    workerid = workers.Get_rank();	  
    cout << "node = " << workerid << endl;
    iter = 0;
    while(!done){
      iter++;
      request = 1;
      world.Recv(rands,CHUNKSIZE,MPI::INT,server,REPLY,status);
      for (i=0;i<CHUNKSIZE;) {		    
	x = (((double) rands[i++])/max) * 2 - 1;
	y = (((double) rands[i++])/max) * 2 - 1;
	if ( x*x + y*y < 1.0 ) in++;
	else out++;		    
      }
      workers.Allreduce(&in,&totalin,1,MPI::INT,MPI::SUM);
      workers.Allreduce(&out,&totalout,1,MPI::INT,MPI::SUM);
      Pi = (4.0*totalin)/(totalin + totalout);
      error = fabs(Pi - 3.1415926);
      done = ((error < epsilon) || (totalin + totalout) > 1000000);
      request = (done) ? 0 : 1;
      if (myid == 0){
	cout << "pi= " << Pi << endl;
	world.Send(&request,1,MPI::INT,server,REQUEST);
      }
      else{
	if (request) world.Send(&request,1,MPI::INT,server,REQUEST);
      }
    }	  
  }

  if (myid == server) cout << "server done" << endl;
  else cout << "before free node= " << workerid << endl;  
 
  if (myid == 0){
    cout << "totalin = " << totalin << endl;
    cout << "totalout = " << totalout << endl;
  }
  //worker_group.Free();
  //workers.Free();
  if (myid != server) cout << "after free node= " << workerid << endl;
  MPI::Finalize();

  return 0;

}


	  
     

