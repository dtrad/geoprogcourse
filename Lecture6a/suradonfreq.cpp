#include "su.h"
#include "segy.h"
#include <signal.h>
#include <math.h>
#include <time.h>
#include "header.h"
void save_gather(float **d, int nh, int nt, float dt, char* name);
/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SURADONFREQ - Frequency domain Radon Transform                      ",
  " 	   								",
  " suradonfreq < stdin > stdout [optional parameters]       		",
  "                                                                     ",
  " sequencial processingfor a line of cdps. The data input has to      ",
  " be sort by cdp and offset, for example:                             ",
  "  susort cdp offset < input | suradonfreq .....                      ",
  "                                                                     ",
  NULL};
/* Credits:
 *	Daniel Trad, UBC, 2000. Last update June, 2016.
 * Trace header fields accessed: ns, cdp, dt, offset
 **************** end self doc ***********************************/

static void closefiles(void);
/* Globals (so can trap signal) defining temporary disk files */
char tracefile[BUFSIZ];	/* filename for the file of traces	*/
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *tracefp;		/* fp for trace storage file		*/
FILE *headerfp;		/* fp for header storage file		*/
FILE *modelfilep;       /* fp for model file output  		*/
///////////////////////////////////////////////////////////////////

int main(int argc, char **argv){
  segy tr, tr2;
  cwp_String modelfile=" "; // output model
  time_t start, finish;
  double elapsed_time;
  
  // Data space
  int nt;   // number of time samples 
  int nhmax; // max number of offsets.
  float dt; // Time sampling
  float dx;  // midpoint interval
  long ntr;  // Total number of input traces

  // Radon space /////////////////////////  
  int nq;  
  float qmin; 
  float quantile;
  int ngettr; // Number of bytes read by gettr (0 after)
  float depth;   // depth for pseudo-hyperbolic RT 
  int rtmethod;
  char *solver;
  // define Radon axis by minimum/maximum moveout at maximum offset
  float moveoutmin; // minimum moveout (ms) at far offset (usually negative)
  float moveoutmax; // maximum moveout (ms) at far offset 
  float maxoffset; // max offset used as a reference for moveout axis
  ///////////////////////////////////////////////////////////////////////

  
  int  cdpmin;          //cdp limits for output
  int  cdpmax;          //cdp limits for output
  int  dcdp;            // interval
  float  dxcdp;         // output cdp interval 
  int verbose;		// flag for echoing info		
  float fmax;

  // Initialize 
  initargs(argc, argv);
  requestdoc(1);

  // Register starting time
  start=time(0);

  // Get info from first trace 
  if (!gettr(&tr)) err("can't read first trace");
  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");

  //if (!tr.offset) err("offset header field must be set");
  dt = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;
  tr2=tr; /* this trace has to be used later in the cdp loop */ 
  fprintf(stderr,"nt=%d,dt=%f\n",nt,dt);  
  if (fmax==0) fmax=1./(2*dt); // default Nyquist

  // Get parameters 
  if (!getparint("cdpmin", &cdpmin))  cdpmin = 1;
  if (!getparint("cdpmax", &cdpmax))  cdpmax = 1;
  if (!getparint("dcdp",&dcdp))         dcdp = 1;

  if (!getparfloat("moveoutmin",&moveoutmin)) moveoutmin=-50;
  if (!getparfloat("moveoutmax",&moveoutmax)) moveoutmax=200;
  if (!getparfloat("maxoffset",&maxoffset))   maxoffset = 1000;
  if (!getparstring("modelfile",&modelfile))  modelfile="model.su";


  float conversion = maxoffset;
  if (rtmethod == 2) conversion = maxoffset*maxoffset;


  // open permanent and temporary files
  modelfilep=efopen(modelfile,"w");
  tracefp = etmpfile();
  headerfp = etmpfile();
  
  // axes
  //for (int i=0;i<nq;i++) q[i]=qmin+i*dq;
  fprintf(stderr,"Starting cdp loop \n");
  int itr = 0; // total number of traces;
  int icdp = 0; // total number of cdps;
  int ix = 0; // number of traces for current cdp;
  int oldcdp=tr.cdp;
  int nx=1000; // initial value;
  float** data=ealloc2float(nt,nx);  
  float* offset=ealloc1float(nx);
  memset(data[0],0,nx*nt*FSIZE);
  bool haveMore=true;
  do{
    // if new trace is new cdp
    fprintf(stderr,"cdp=%d\n",tr.cdp);   
    icdp++;
    erewind(headerfp);
    // store current cdp to temporary file
    ix=0;
    do{
      efwrite(&tr,HDRBYTES,1,headerfp);
      memcpy(data[ix],tr.data,nt*FSIZE);
      offset[ix]=tr.offset;
      ix++;
    } while ((gettr(&tr))&&(tr.cdp == oldcdp));
    if (tr.cdp == oldcdp) haveMore=false; // end of data break;
    oldcdp=tr.cdp; // keep track of last trace;
    nx=ix++; // ntraces for this cdp;
    // Add Radon here
    save_gather(data,nx,nt,dt,"data");
    // Output filtered data;
    erewind(headerfp);
    for (ix=0;ix<nx;ix++){
      efread(&tr2,HDRBYTES,1,headerfp);
      memcpy(tr2.data,data[ix],nt*FSIZE);
      //fprintf (stderr,"ix=%d\n",ix);
      puttr(&tr2);
    }
    itr+=nx;
    fprintf(stderr,"itr=%d, current cdp=%d nx=%d\n",itr,icdp,nx);

  }while (haveMore);
  


  finish=time(0);
  elapsed_time=difftime(finish,start);
  fprintf(stderr,"Total time required: %f \n", elapsed_time);
  return EXIT_SUCCESS;
}

/* for graceful interrupt termination */
static void closefiles(void)
{
	efclose(headerfp);
	efclose(tracefp);
	eremove(headerfile);
	eremove(tracefile);
	exit(EXIT_FAILURE);
}
