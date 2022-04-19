//******************************************************************************
// mpi example 5: Geophysical programming University of Calgary
// helloworld, based on an example from Acceleware.
// credits: Acceleware, Daniel Trad 
//*****************************************************************************/

#include <stdio.h>
#include <mpi.h>
#include <iostream>
using namespace std;
#define BUFSIZE 128
#define TAG 0

int main(int argc, char* argv[])
{
    int rank, size;
    char processorName[MPI_MAX_PROCESSOR_NAME];
    int version, subversion;
    int length;
    // Initialize MPI
    MPI_Init(&argc, &argv);
    // Get the rank of this process, and
    // the total number of processes in this communicator
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    
    // Get the MPI version number (MPI_Get_version) from rank 0 and print all
    if (rank == 0){
      MPI_Get_version(&version,&subversion);
      printf("master rank %d size =%d version=%d subversion=%d \n", 
	     rank,size,version, subversion);
    }

    // Get the processor name for this process    
      MPI_Get_processor_name(processorName,&length);
    
    // Print processor name, rank, and total number of processes
    processorName[length]='\n';
    printf("rank = %d name = %s",rank,processorName);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 1) cout << "rank=" << rank << "processor" << processorName << endl;
    // Finalize MPI
    MPI_Finalize();
    return 0;
}
