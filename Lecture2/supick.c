/* SUPICK: GOPH 699 UofC	*/

#include "su.h"
#include "segy.h"
#include <math.h>
#include <stdbool.h>
/*********************** self documentation **********************/
char *sdoc[] = {
  " 								",
  " SUPICK - example of simple SU program.                      ",
  " supick < stdin > stdout					",
  "                                                             ",
  " Required parameters:					",
  "        none			  			        ",
  " Optional parameters:                                        ",
  " float small=0         value considered equivalent to zero   ",
  " int   verbose=0       level of detail to print out          ",
  " Example:                                                    ",
  "       supick < input.su > output.su small=0.0001 verbose=1  ",
  NULL};

/* Credits:
   Daniel Trad, UofC
 */
/**************** end self doc ***********************************/

// Usually is not good idea to have global variables but in small programs is reasonable.
segy tr; // define a segy trace to use for input/output
int verbose; // level of detail for debugging info.

// example of a function to find, store and print first sample on each trace.
void firstbreak(unsigned short dt, float** data, int nx, int nt, float *time, float small){
  for (int ix=0;ix<nx;ix++){
    int it=0;
    while ((fabs(data[ix][it]) <= small)&&(it<nt)) it++;
    time[ix]=it*dt;
    if (verbose>1) fprintf(stderr,"time=%f --> data[%d][%d]=%f\n",time[ix],ix,it,data[ix][it]);
  }
  return;
}

int main(int argc, char **argv){
  /* Initialize */
  initargs(argc, argv);
  requestdoc(1); /* two file args required */
  float small;   // threshold to find zero value

  if (!getparfloat("small",&small))   small=0.;
  if (!getparint("verbose",&verbose)) verbose=1;

  // SU function to check input parameters
  checkpars();

  // open a file to write maximum for each trace an location.
  FILE *fp1 = efopen("list", "w"); 

  // loop over traces. Read first trace to get info and rewind 
  fgettr(stdin,&tr); // read first trace for stdin, 
  erewind(stdin);    // rewind to the begining of stdin

  // get info from trace
  if (!tr.dt) err("error  tr.dt not set ");
  if (!tr.ns) err("error  tr.ns not set");
  unsigned short dt = tr.dt*0.001; // convert from microsec to milisec
  int nt=tr.ns;  // number of samples per trace
  int nx=0; // number of traces


  // check we have required information in input data
  if (!tr.ntr){
    // ntr is not set, read all data set and count traces.
    if (verbose) fprintf(stderr,"ntr is not set, counting traces...\n");
    int itr=0;
    while (gettr(&tr)) itr++;
    nx = itr;
    erewind(stdin);
  }
  else nx=tr.ntr;

  // create a matrix for the data set. 
  float** data=ealloc2float(nt,nx); // reverse convention in SU, data[nx][nt] 
  float*  offset=ealloc1float(nx);  // vector to store offsets
  float*  time=ealloc1float(nx);    // vector to store first times
  int itr = 0;	// number of the trace being processed at each loop
  while (gettr(&tr)){ // do while there are input traces
    for (int it=0;it< nt;it++)  data[itr][it]=tr.data[it];
    offset[itr]=(float) tr.offset;
    ++itr;
  }
  
  // call function to do some matrix operation in data.
  firstbreak(dt,data,nx,nt,time,small);


  // put back data into traces after some processing. 
  int count = 0;
  for (int ix=0;ix<nx;ix++){
    if (verbose > 1) fprintf(stderr,"offset=%f time=%f \n",offset[ix],time[ix]);    
    // write info to file
    fprintf(fp1,"offset=%f time=%f \n",offset[ix],time[ix]);
        
    for (int it=0;it< nt;it++) tr.data[it]=data[ix][it];
    tr.data[NINT(time[ix]/dt)]=-2;
    puttr(&tr); count++;
  }

  fprintf(stderr,"End of run, %d traces in, %d traces out \n",nx,count);

  fclose(fp1);
  free2float(data);
  free1float(offset);
  free1float(time);
	
  return(CWP_Exit());
}
