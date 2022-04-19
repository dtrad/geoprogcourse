#include "su.h"
#include "segy.h"
#include "header.h"
#include <mpi.h>
#define MARK fprintf(stderr,"rank=%d %s @ %u\n",rank,__FILE__,__LINE__)
char *sdoc[] = {
" 								",
" SUMPITEST - Seismic Unix MPI TEST                             ",
" Example: sumpitest < input.su > output.su                     ",
NULL};
/**************** end self doc ***********************************/
segy tr;
int main(int argc, char** argv){
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  if (rank==0){
    // read data fron stdin and send it to node 1
    initargs(argc,argv);
    for (int i=0;i<argc;i++) fprintf(stderr,"argc=%d argv[%d]=%s\n",argc,i,argv[i]);
    FILE* fp1 = efopen(argv[1],"r");
    while(fgettr(fp1,&tr)){
      fprintf(stderr,"rank=%d inputing trace tr.tracl=%d\n",rank,tr.tracl);
      MPI_Send(&tr.ns,1,MPI_INT,1,0,MPI_COMM_WORLD);
      MPI_Send(tr.data,tr.ns,MPI_FLOAT,1,1,MPI_COMM_WORLD);
    }
    int zero = 0;
    MPI_Send(&zero,1,MPI_INT,1,0,MPI_COMM_WORLD);
    fclose(fp1);
  }
  if (rank==1){
    FILE* fp2 = efopen(argv[2],"w");
    tr.tracl = 0;
    MPI_Recv(&tr.ns,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    while(tr.ns){
      MPI_Recv(tr.data,tr.ns,MPI_FLOAT,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      tr.tracl++;
      fprintf(stderr,"rank=%d outputting trace tr.tracl=%d\n",rank, tr.tracl);
      fputtr(fp2,&tr);
      MPI_Recv(&tr.ns,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    fclose(fp2);
  }

  MPI_Finalize();
  return EXIT_SUCCESS;
}
