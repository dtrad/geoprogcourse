#include <iostream>
#include <mpi.h>
// pack an integer and a char array into a buffer for mpi
using namespace std;
#define MARK fprintf(stdout,"%s @ %u\n",__FILE__,__LINE__)
int main(int argc, char ** argv){
     float i=0;
     int sizechar = 100;
     int sizebuffer = sizechar + sizeof(int);

     char c[sizechar];
     char buffer[sizebuffer];
     int position = 0; // starting position in buffer
     MPI::Status status;

     MPI::Init(argc,argv);
     MPI::Intracomm world = MPI::COMM_WORLD;
     int size = world.Get_size();
     if (size < 2){
       cerr << " need at least two nodes" << endl;
       return 1;
     }
     int rank = world.Get_rank();
     if (rank == 0){
       memset(c,0,sizeof(char)*sizechar);
       c[0]='m';
       c[1]='e';
       c[2]='s';
       c[3]='s';
       i = 10;
       MPI::INT.Pack(&i,1,buffer,sizebuffer,position,world);
       cout << "position=" << position << endl;
       MPI::CHAR.Pack(c,sizechar,buffer,sizebuffer,position,world);
       cout << "position=" << position << endl;
       world.Send(buffer,position,MPI::PACKED,1,0);
     }
     else if (rank == 1){
       world.Recv(buffer,sizebuffer,MPI::PACKED,0,0,status);
       MPI::INT.Unpack(buffer,sizebuffer,&i,1,position,world);
       MPI::CHAR.Unpack(buffer,sizebuffer,c,sizechar,position,world);
       for (int j=0;j<sizechar;j++) cout << c[j];
       cout << endl;
       cout << "i=" << i << endl;
     }
     
     MPI::Finalize();
     return 0;

}


	  
     

