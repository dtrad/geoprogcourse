#include <iostream>
#include <mpi.h>
using namespace std;
// creating a new datatype for MPI
#define MARK fprintf(stdout,"%s @ %u\n",__FILE__,__LINE__)
void printA(double *A, int n, int m){
  cout << "A = " << endl;
  for (int i=0;i<n;i++){
    for (int j=0;j<m;j++){
      cout << A[i*n+j] << " ";
    }
    cout << endl;
  }
}

int main(int argc, char ** argv){
     float x=0;
     int i=0;
     MPI::Status status;

     MPI::Init(argc,argv);
     MPI::Intracomm world = MPI::COMM_WORLD;
     int size = world.Get_size();
     int rank    = world.Get_rank();
     if (size < 3){
       cerr << " needs at least 3 processes " << endl;
     }
     int ntest=10;
     int sizeofa= ntest*ntest;
     double *a = new double[sizeofa];
     memset(a,0,sizeof(double)*sizeofa);
     int disp[ntest];
     int blocklen[ntest];
     for (i=0;i<ntest;++i){
	  disp[i]=ntest*i+i;
	  blocklen[i]=ntest-i;
	  if (!rank){
	    a[i*ntest+i]=9;
	    if (i< ntest-1) a[i*ntest+i+1]=8;
	  }
     }
     // print contents for different processes for test (remove later)
     bool test = 1;
     if (test){
       if (!rank){ 
	 cout << " in rank = " << rank << endl;
	 printA(a,ntest,ntest);
       }
       world.Barrier();
       if (rank==1){ 
	 cout << " in rank = " << rank << endl;
	 printA(a,ntest,ntest);
       }
       world.Barrier();
     }

     MPI::Datatype upper;
     upper = MPI::DOUBLE.Create_indexed(ntest,blocklen,disp);
     upper.Commit();
	  
     if (rank == 0){
	  i = 10;
	  cout << "rank " << rank << ":sending a"<< endl;
	  world.Send(a,1,upper,2,0);
	  //world.Send(a,10000,MPI::DOUBLE,2,0);
	  cout << "rank " << rank << ":sent a"<< endl;
     }
     else if (rank == 1){
	  x = 1.1;
	  cout << "rank=" << rank << endl;
	  world.Send(&x,1,MPI::FLOAT,2,0);
	  cout << "rank " << rank << ":sent x"<< endl;
     }
     else if (rank == 2){
	  cout << "rank=" << rank << endl;
	  for (int ix=1;ix<=2;ix++){
	       world.Probe(MPI::ANY_SOURCE,0,status);
	       if (status.Get_source() == 0){
		    cout << "receiving a " << endl;
		    world.Recv(a,1,upper,0,0);
		    //world.Recv(a,10000,MPI::DOUBLE,0,0);
		    cout << "received a" << endl;
	       }
	       else if (status.Get_source() == 1){
		    cout << "receiving x " << endl;
		    world.Recv(&x,1,MPI::FLOAT,1,0);
		    cout << "received x " << endl;
	       }
	  }
     }
     world.Barrier();
     if (rank == 2){ 
       cout << " in rank = " << rank << endl;
       cout << " x = " << x << endl;
       printA(a,ntest,ntest);
     }

     MPI::Finalize();
     return 0;

}


	  
     

