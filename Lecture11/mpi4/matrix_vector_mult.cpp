#include <iostream>
#include <mpi.h>
#include <math.h>
#include "MyAlloc.hpp"
//******************************************************************************
// mpi example 4: Geophysical programming University of Calgary
// Matrix vector multiplication
// credits: Daniel Trad (based on Using MPI by Gropp, Lusk, Skjellum
//*****************************************************************************/
#ifndef MARK
#define MARK fprintf(stdout,"%s @ %u\n",__FILE__,__LINE__)
#endif

int main(int argc, char ** argv){
     int rank, size, i, j, ip, master, anstype, sender;
     int nrows, ncols;
     double **A;
     double  *b;
     double  *c;
     double  *buffer, ans;     

     int numsent = 0;
     int rec = 0;
     // assign parameters
     nrows = 5000;
     ncols = 50000;
     master = 0;
	  
     MPI::Init(argc, argv);
     size = MPI::COMM_WORLD.Get_size();
     rank = MPI::COMM_WORLD.Get_rank();
     MPI::Status status;
     
     cout << "rank = " << rank << endl;
     cout << "node = " << rank << ": size = " << size << endl;
     
     new1d(b,ncols);
     new1d(buffer,ncols);	  
     if (rank == master){
	  new2d(A,nrows,ncols);
	  new1d(c,nrows);
	  for (j=0; j<ncols;j++){
	       b[j]=j;
	       for (i=0; i<nrows;i++){
		    A[i][j]=j;
	       }
	  }
	  MPI::COMM_WORLD.Bcast(b,ncols,MPI::DOUBLE,master);
	  int np = nrows;
	  if ((size-1) < nrows) np = size - 1;
	  for (ip=1;ip<= np;ip++){
	       for (j=0;j<ncols;j++){
		    buffer[j]=A[ip-1][j];
	       }
	       MPI::COMM_WORLD.Send(buffer,ncols,MPI::DOUBLE,ip,ip);
	       numsent++;
	       cout << "master: numsent = " << numsent << " ip = " << ip << endl;
	  }
	  cout << "master: send first finish" << endl;
	  for (ip=1;ip<=nrows;ip++){
	       MPI::COMM_WORLD.Recv(&ans,1,MPI::DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,status);
	       sender = status.Get_source();
	       anstype = status.Get_tag();
	       c[anstype]= ans;
	       if (numsent < nrows){
		    for (j=0;j<ncols;j++){
			 buffer[j]=A[numsent][j];
		    }
		    MPI::COMM_WORLD.Send(buffer,ncols,MPI::DOUBLE,sender,numsent+1);
		    numsent++;
		    cout << "master: numsent = " << numsent << " ip = " << ip << endl;
	       }
	       else{
		    MPI::COMM_WORLD.Send(MPI_BOTTOM,0,MPI::DOUBLE,sender,0);
	       }
	  }
	  cout << "master: send second finish" << endl;
     }
     else{
	  MPI::COMM_WORLD.Bcast(b,ncols,MPI::DOUBLE,master);
	  cout << "node = " << rank << " received b" << endl;
	  if (rank > nrows){
	       cout << "nothing to do" << endl;
	  }
	  else{
	       int tag =1;
	       while (tag != 0){
		    MPI::COMM_WORLD.Recv(buffer,ncols,MPI::DOUBLE,master,MPI_ANY_TAG,status);
		    cout << "node = " << rank << ":received " << ++rec << endl;
		    tag = status.Get_tag();
		    if (tag == 0){
			 cout << "done" << endl;
			 break;
		    }
		    else{
			 ans = 0;
			 for (j=0;j<ncols;j++) ans+= (buffer[j]*b[j]);
			 MPI::COMM_WORLD.Send(&ans,1,MPI::DOUBLE,master,tag);
			 cout << "node = " << rank << ":sent back " << endl;
		    }
	       }
	  }
     }
     
     cout << "node =" << rank << " is done " << endl;
     
     
     if (rank == master){
	  del1d(c);
	  del2d(A);
     }
     MPI::Finalize();
     del1d(b);
     del1d(buffer);	  

     return 0;

}


	  
     

