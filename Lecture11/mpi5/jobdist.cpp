#include <iostream>
#include <mpi.h>
#include <vector>
//******************************************************************************
// mpi example 5: Geophysical programming University of Calgary
// job distribution example
// credits: Daniel Trad 
//*****************************************************************************/
using namespace std;

#define  N 10000000000
#ifndef MARK
#define MARK fprintf(stdout,"%s @ %u\n",__FILE__,__LINE__)
#endif

#define COMM MPI_COMM_WORLD

int doSomething(int inode, int ijob){
  cout << "node " << inode << " working on job " << ijob << endl;
  for (unsigned long int i=0;i< N;i++);
  return 0;
}


int main(int argc, char ** argv){
  int myNode;
  int totalnodes;
  int data = 0;
  int busy=0;
  int done=0;
  int inode;
  int njobs =9;
	  
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(COMM, &totalnodes);
  MPI_Comm_rank(COMM, &myNode);
     
  if (!busy) cout << "Node number = " << myNode << " is active and waiting" << endl;
     
  if (myNode == 0){
    // master send information to each node and wait for result
    int ijob = 1;
    while (ijob < njobs){
      cout << "total nodes = " << totalnodes << endl;
      for (inode = 1; inode < totalnodes; inode++){
	data = ijob;
	cout << "sending job " << ijob << " to node " << inode << endl;
	MPI_Send(&data,1,MPI_INT,inode,inode,COMM);
	ijob++;
      }
	       
      for (inode = 1; inode < totalnodes; inode++){
	MPI_Recv(&data,1,MPI_INT,inode,inode,COMM,&status);
	cout << "receiving job " << data << "from node " << inode << endl;
      }
    }
    done = 1;
  }
     
  if (myNode != 0){
    // slaves receive data, process it and send result back
    while(data < njobs - (totalnodes - myNode)){
      MPI_Recv(&data,1,MPI_INT,0,myNode,COMM,&status);
      busy = 1;
      busy = doSomething(myNode, data);
      MPI_Send(&data,1,MPI_INT,0,myNode,COMM);
      cout << "myNode = " << myNode << "--> data = " << data << endl;
    }
  }
  // everybody: broad cast signal done (synchronization)
  // Wait till all  processes have finished.
  MPI_Bcast(&done,1,MPI_INT,0,MPI_COMM_WORLD);     
  // send info to log that all nodes are done. 
  for (inode = 0; inode < totalnodes; inode++){
    if (myNode == inode){
      cout << "finalize node " << myNode << "with done = " << done << endl;
	       
    }
  }

  MPI_Finalize();		      

}


	  
     

