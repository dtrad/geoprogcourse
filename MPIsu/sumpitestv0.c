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
    MARK;
    initargs(argc,argv);
    for (int i=0;i<argc;i++) fprintf(stderr,"argc=%d argv[%d]=%s\n",argc,i,argv[i]);
    FILE* fp1 = efopen(argv[1],"r");
    while(fgettr(fp1,&tr)){
      fprintf(stderr,"rank=%d inputing trace tr.tracl=%d\n",rank,tr.tracl);
    }
    fclose(fp1);
  }
  if (rank) MARK;
  if (rank==1){
    MARK;
    FILE* fp1 = efopen(argv[1],"r");
    FILE* fp2 = efopen(argv[2],"w");

    while(fgettr(fp1,&tr)){
      fprintf(stderr,"rank=%d outputting trace tr.tracl=%d\n",rank, tr.tracl);
      fputtr(fp2,&tr);
    }
    fclose(fp1);
    fclose(fp2);
  }
  if (rank) MARK;

  MPI_Finalize();
  return EXIT_SUCCESS;
}
