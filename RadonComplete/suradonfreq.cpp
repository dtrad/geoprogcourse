/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* SUEOMIG:  $Date: March 1999  */

#include "su.h"
#include "segy.h"
#include "header.h"
#include "radonfrequency.h"
#include <signal.h>
#include <math.h>
#include <time.h>
#include "Mute.hpp"


/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SURADONFREQ - Frequency domain Radon Transform                      ",
  " 	   								",
  " suradonfreq < stdin > stdout [optional parameters]       		",
  "                                                                     ",
  " Multiple attenuation for a line of cdps. The data input has to      ",
  " be sort by cdp and offset, for example:                             ",
  "  susort cdp offset < input | suradonfreq .....                      ",
  "                                                                     ",
  " Notes:                                                              ",
  " Because NMO is applied internally the module requires the velocities",
  " Velocities in SU are supplied with par files with format like this  ",
  " cdp=1,20                                                            ",
  " vnmo=1700,2000,2900,3000,3500                                       ",
  " tnmo=0.5,0.7,1,1.3,1.5                                              ",
  " vnmo=1700,2200,2900,3300,3500                                       ",
  " tnmo=0.5,0.7,1,1.3,1.6                                              ",
  " the module interpolates velocities between cdp numbers.             ",
  " Also the parameter nmofactor allows one to control the correction.  ",
  " i.e nmofactor=0 will do no correction, 1 will do full correction    ",
  "                                                                     ",  
  " solver defines whether to have high resolution RT (cgfft, wtcgls)   ",
  " or standard RT (toep). Option adj skips the solver.                 ",
  " rtmethod defines if working with linear (1) , parabolic (2) o       ",
  " pseudohyperbolic RT (3). In this last case need to set also depth   ",
  " (hyperbolic for one target depth)                                   ",
  "                                                                     ",
  " cdpmin=1        First CDP                                           ",
  " cdpmax=1        Last CDP 						",
  " dcdp=1          Interval in cdps                                    ",
  " dxcdp=50        cdp interval in feet or meters of cdps to process   ",
  " par=            parfile with cdp number, nmo velocities and nmotime ",
  "                 as obtained  from Velan (standard PARFILE in SU)    ",
  " rtmethod=2      shape of integration surface                        ",
  "                 =1 linear                                           ",
  "                 =2 parabolic                                        ",
  "                 =3 Pseudo-hyperbolic (Foster and Mosher)            ",
  " depth=1000      only used for rtmethod=3                            ",
  " solver=cgfft    Method to solve the system of equations             ",
  "                 =cgfft very fast cg with fft                        ",
  "                 =wtcgls standard cg with weight functions           ",
  "                 =toep  Toeplitz solver (no high resolution)         ",
  "                 =adj Simple adjoint                                 ",
  " modelfile=model.su  optional filename to output RT model            ",
  " nhmax=200       maximum number of traces for cdp                    ",
  " moveoutmin=-50  minimum moveout at maxoffset                        ", 
  " moveoutmax=200  maximum moveout at maxoffset                        ", 
  " maxoffset=1550  maximum offset for conversions to Radon parameter   ",
  " itercg=50       Internal iterations for CG                          ",
  " iter_end=3      External iterations for IRLS                        ",
  " eps1=5e-1       Numerator hyperparameter for Wd (quantil of data)   ",
  " eps2=5e-1       Denominator  hyperparameter for Wm (quantil of model)",
  " fmax=0.8/(2*dt) Maximum frequency to preserve in the transform      ",
  "                 =0 Nyquist (1/2*dt)                                 ", 
  " norm=0          model weights derivated from Cauchy norm            ",
  "                 =1 model weights derivated from L1 norm             ", 
  " mute=0          =1 mute multiples                                   ",
  " verbose=0       =1  verbose output                                  ",
  " t0=0            First useful time                                   ",
  " smute=2         stretch greater than smute x 100 % is muted		",
  " nmofactor=1     nmofactor * offset is used for NMO	                ",
  " quantile=1       filter for large data outliers                     ",
  "                                                                     ",
  " mute parameters :                                                   ",
  " fbmute=0        =1 apply first break mute after inverse             ",
  " tmin_m=0        min time for pass mute window                       ",
  " tmax_m=nt*dt    max time for pass mute window                       ",
  " moveoutmin_m=moveoutmin minimum moveout for mute window             ",
  " moveoutmax_m=moveoutmax maximum moveout for mute window             ",
  " thres_m=0.2    values less than threshold are removed inside window ",
  " slope_m=3      slope in pass mute window (ih top / ih bottom)       ",
  "                                                                    ",
  NULL};
/* Credits:
 *	Daniel Trad, UBC, 2000. Last update June, 2016.
 * Trace header fields accessed: ns, cdp, dt, offset
 **************** end self doc ***********************************/
static void closefiles(void);
int getFirstBreak(float* data, int nt, float threshold, int minIt);
/* Globals (so can trap signal) defining temporary disk files */
char tracefile[BUFSIZ];	/* filename for the file of traces	*/
char headerfile[BUFSIZ];/* filename for the file of headers	*/
FILE *tracefp;		/* fp for trace storage file		*/
FILE *headerfp;		/* fp for header storage file		*/
FILE *modelfilep;       /* fp for model file output  		*/
///////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  segy tr,tr2;
  cwp_String modelfile=""; /* output sufile for the model */ 
  ///////////////////////////////////////////////////////////////
  time_t start,finish;
  double elapsed_time;
  
  // Data space
  int nt;   // number of time samples 
  int nx;    // number of midpoints
  int nhmax; // max number of offsets.
  float dt; // Time sampling
  float dx;  // midpoint interval

  
  ///////////////////////////////////////////////////////////////
  int ih, ix; // General counters 
  register int it;
  
  float t0=0;      // First useful time
  unsigned int ntrmax;  // Numebr of traces to process from data file
  /* Velocity */
  int ncdp;	/* number of cdps specified */
  

  char *tmpdir;		/* directory path for tmp files		*/
  cwp_Bool istmpdir=cwp_false;/* true for user-given path		*/ 
  long ntr;     // Total number of input traces
  int itr;
  int  cdpmin;          //limit for the output
  int  cdpmax;          //limit for the output
  int  dcdp;
  float  dxcdp;           // output cdp interval 
  int verbose;		/* flag for echoing info		*/
  float fmax;
  // Radon space  
  int nq;  
  float qmin;
  
  //  LS mig  
  float step;
  float eps1;
  float eps2;
  int itercg;
  int iter_end;
  int norm;
  
  float factor;
  float smute;
  float nmofactor;
  
  float quantile;
  int ngettr; // Number of bytes read by gettr (0 after)
  float depth;   // depth for pseudo-hyperbolic RT 
  int rtmethod;
  char *solver;
  
  // define Radon axis by minimum/maximum moveout at maximum offset
  float moveoutmin; // minimum moveout (ms) at far offset (usually negative)
  float moveoutmax; // maximum moveout (ms) at far offset 
  float maxoffset; // max offset used as a reference for moveout axis
  int fbmute; // if 1 apply first break mute;
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
  
  // Get parameters 
  
  if (!getparint("cdpmin", &cdpmin))  cdpmin = 1;
  if (!getparint("cdpmax", &cdpmax))  cdpmax = 1;
  if (!getparint("dcdp",&dcdp))         dcdp = 1;
  if (!getparfloat("dxcdp", &dxcdp))  dxcdp = 50;
  if (!getparint("nhmax", &nhmax))  nhmax = 300;
  if (!getparint("verbose", &verbose)) verbose = 0;    
  if (!getparfloat("eps1", &eps1))  eps1 = 5e-1;
  if (!getparfloat("eps2", &eps2))  eps2 = 5e-1;
  if (!getparfloat("step", &step))  step = 0.98;
  if (!getparint("itercg", &itercg))  itercg = 50;
  if (!getparint("iter_end", &iter_end))  iter_end = 3;
  if (!getparint("norm",&norm)) norm=0;
  if (!getparuint("ntrmax",&ntrmax)) ntrmax = 1000000;
  if (!getparfloat("t0",&t0)) t0=0;
  if (!getparfloat("factor",&factor)) factor=0.8;
  if (!getparfloat("smute",&smute)) smute=2;
  if (!getparfloat("nmofactor",&nmofactor)) nmofactor=1;
  if (!getparstring("modelfile",&modelfile)) modelfile="model.su";
  if (!getparfloat("quantile",&quantile)) quantile=1;
  if (!getparint("rtmethod",&rtmethod)) rtmethod=2;
  if (!getparfloat("depth",&depth)) depth=2000;
  if (!getparfloat("fmax",&fmax)) fmax=0.8/(2*dt);
  if (!getparstring("solver",&solver)) solver="cgfft";
  if (!getparfloat("moveoutmin",&moveoutmin)) moveoutmin=-50;
  if (!getparfloat("moveoutmax",&moveoutmax)) moveoutmax=200;
  if (!getparfloat("maxoffset",&maxoffset))   maxoffset = 1550;
  if (!getparint("fbmute",&fbmute)) fbmute=0;
  /* Format the string solver to have the same length */ 
  if      (STREQ(solver,"toep"))   solver="toep__";
  else if (STREQ(solver,"wtcgls")) solver="wtcgls";         
  else if (STREQ(solver,"cgfft"))  solver="cgfft_";
  else if (STREQ(solver,"adj"))    solver="adj___";
  if (fmax==0) fmax=1./(2*dt);
  
    /* Specify mute     ********************************************************/  
  /* for mute in the migrated space we need to set the geometry of the mask */
  int mute; 
  mutemask_par mutepar;
  int moveoutmin_m;
  int moveoutmax_m;
  if (!getparint("mute",&mute)) mute = 0;
  if (!getparfloat("tmin_m",&mutepar.tmin)) mutepar.tmin = 0.; 
  if (!getparfloat("tmax_m",&mutepar.tmax)) mutepar.tmax = nt*dt;
  // for moveout, give in ms, need to convert to trace number once dq is known.
  if (!getparint("moveoutmin_m",&moveoutmin_m)) moveoutmin_m = moveoutmin;
  if (!getparint("moveoutmax_m",&moveoutmax_m)) moveoutmax_m = moveoutmax;
  if (!getparfloat("slope_m",&mutepar.slope)) mutepar.slope = 0;     
  if (!getparfloat("thres_m",&mutepar.threshold)) mutepar.threshold = 0.2;
  
  float conversion = maxoffset;
  if (rtmethod == 2) conversion = maxoffset*maxoffset;
  
  Mute* myMute = 0;
  if (mute) myMute = new Mute(mutepar,moveoutmin_m,moveoutmax_m);
  
  /************************************************************************/ 
  
  modelfilep=efopen(modelfile,"w");
  nx=(cdpmax-cdpmin)/dcdp + 1;

  if (verbose) fprintf(stderr,"number of ouput cdps = %d\n",nx);

  /* get velocity functions, linearly interpolated in time */
  // Velint will contain velocity law for one csp
  ncdp = countparval("cdp");


  float* cdp = ealloc1float(ncdp+1); // get velocities crashes if not added +1
  float** ovv = ealloc2float(nt,ncdp+1); // array of slowness (1/velocity^2) 
  float* velint=ealloc1float(nt); // array of velocities for each cdp
  getvelocities(dt,nt,ncdp,cdp,ovv);

  /* Look for user-supplied tmpdir */

  if (!getparstring("tmpdir",&tmpdir) &&
      !(tmpdir = getenv("CWP_TMPDIR"))) tmpdir="";
  if (!STREQ(tmpdir, "") && access(tmpdir, WRITE_OK))
    err("you can't write in %s (or it doesn't exist)", tmpdir);  
 
  // Store CMP traces in the disk for later use
  if (STREQ(tmpdir,"")){
    tracefp = etmpfile();
    headerfp = etmpfile();
    if (verbose) warn("using tmpfile() call");
  }else{ /* user-supplied tmpdir */
    char directory[BUFSIZ];
    strcpy(directory, tmpdir);
    strcpy(tracefile, temporary_filename(directory));
    strcpy(headerfile, temporary_filename(directory));
    /* Trap signals so can remove temp files */
    signal(SIGINT,  (void (*) (int)) closefiles);
    signal(SIGQUIT, (void (*) (int)) closefiles);
    signal(SIGHUP,  (void (*) (int)) closefiles);
    signal(SIGTERM, (void (*) (int)) closefiles);
    tracefp = efopen(tracefile, "w+");
    headerfp = efopen(headerfile, "w+");
    istmpdir=cwp_true;		
    if (verbose) warn("putting temporary files in %s", directory);
  }

  // before allocating need to calculate nq.
  // user gives min/max moveout in ms, convert to q parameter;
  moveoutmin*=1e-3; // convert to secs
  moveoutmax*=1e-3;
  float aperture=2*maxoffset;
  float dh=dxcdp;
  float dq;
  float qmax;
  radon_param_init(rtmethod,fmax,factor,aperture,dh, 
                   moveoutmin,moveoutmax,maxoffset,dq,nq,qmin,qmax);
  fprintf(stderr,"nq=%i dq=%g qmin=%g qmax=%g maxoffset=%f \n",nq,dq,qmin,qmax,maxoffset);
  // Allocate memory for data and model
  // SU alloc has reversed dimension parameters (nfast first)
  float** mm=ealloc2float(nt,nq);
  float** cmp=ealloc2float(nt,nhmax);
  float* d=ealloc1float(nt);
  float* x=ealloc1float(nx);
  float* h=ealloc1float(nhmax);
  float* t=ealloc1float(nt);
  float* q=ealloc1float(nq);
  float* Wm=ealloc1float(nt*nq);
  int* nhcmp=ealloc1int(nx);
  int* hmute=alloc1int(nhmax); // store the first non-zero sample
  memset(hmute,0,nhmax*sizeof(int));
  
  if (myMute){
      myMute->init(maxoffset,rtmethod,moveoutmin,dq,nq,dt,nt);
      myMute->write_curve_mute();
  }
  
 
  /***************************************/
  for (int i=0;i<nq;i++) q[i]=qmin+i*dq;

  /* Create axis for output cpds, equivalent offset and time */
  dx=dxcdp;
  for(ix=0;ix<nx;ix++) x[ix]=cdpmin+dcdp*ix;
  for(it=0;it<nt;it++) t[it]=t0+it*dt;

  
  

  ///////////////////////////////////////////////////////////////////////

  fprintf(stderr,"Starting cdp loop ****\n");
  fprintf(stderr,"conversion before =%f \n",conversion); 
  /*************************************************************/
  /* Loop on the CMP gathers  */
  for (int ix=0;ix<nx;ix++){
    /* put back last trace */
    tr=tr2;    
    // search for the first cdp  
    while (tr.cdp<x[ix]) if (!(ngettr=gettr(&tr))) break;
    
    fprintf(stderr,"tr.cdp=%d  -> cdp=%f \n",tr.cdp,x[ix]);

    /* Initialize to zero the ouput trace and the cmp gather corresp to ix */  
    for (ih=0;ih<nhmax;ih++) memset(cmp[ih],(int) '\0',nt*FSIZE);

    erewind(tracefp);
    erewind(headerfp);   
    /**************************************************/    
    // Save the current cdp only (x[ix]) into a temporal file
    ntr = 0;
    
    do {
      ntr++;
      if (verbose) fprintf(stderr,"ntr=%d,tr.cdp=%d,tr.offset=%d\n",
			   ntr,(int) tr.cdp,(int) tr.offset);
      efwrite(&tr,HDRBYTES,1,headerfp);
      efwrite(tr.data,FSIZE, nt, tracefp);
    } while ( (ngettr=gettr(&tr)) && (tr.cdp == x[ix]));
    fprintf(stderr,"x=%f, ntr=%d\n",x[ix],ntr);
    /* Save last trace for next cdp */
    tr2=tr;    
    /**************************************************/
    // Read the current cdp
    erewind(tracefp);
    erewind(headerfp);  
    
    if (ntr<3) continue; // empty cdp
    float threshold = 1;
    int minIt=0;
    for(ih=0;ih<ntr;ih++){
      efread(tr.data, FSIZE, nt, tracefp);
      efread(&tr,HDRBYTES,1,headerfp);
      h[ih]=tr.offset;      
      if (fbmute) hmute[ih]=getFirstBreak(tr.data,nt,threshold,minIt);
      //minIt=hmute[ih];
      memcpy(cmp[ih],tr.data,nt*FSIZE);
    }
    nhcmp[ix]=ntr;


    ///////////////////////////////////////////////////
    /* compute new square slowness  */
    if (ncdp) interpovv(nt,ncdp,cdp,ovv,x[ix],velint);
    float** mask=0;
    if (myMute) mask=myMute->getMask();
    radonsolver0_subtract(cmp,h,nhcmp[ix],t,nt,dt,mm,q,nq,velint,itercg,iter_end,step,
		      eps2,eps1, quantile,norm,factor,smute,nmofactor,rtmethod,depth,
		      fmax,solver,mask,mute); 

    // Output trace
    //erewind(tracefp);
    erewind(headerfp);  
    
    for (ih=0;ih<nhcmp[ix];ih++){ 
      efread(&tr,HDRBYTES,1,headerfp);
      for (int it=0;it<hmute[ih];it++) cmp[ih][it]=0;
      memcpy((void *) tr.data,(const void *) cmp[ih],nt*sizeof(float));      
      if (verbose) fprintf(stderr,"Original data ===> tr.offset=%d\n",tr.offset);
      puttr(&tr);
    }
    
    fprintf(stderr,"nhcmp[%d]=%d  conversion=%f\n",ix,nhcmp[ix],conversion);    
    
    for (itr=0;itr<nq;itr++){
      for (it=0;it<nt;it++) tr.data[it]=mm[itr][it];
      tr.cdp=(int) x[ix]; // front of the CMP
      tr.dt=(int) (dt*1e6);       
      tr.ntr=nq;
      tr.ns=nt;
      tr.tracl=itr+1;
      tr.tracr=itr+1;
      tr.f2=q[itr]; // radon transform axis (sec/offset  or sec/offset**2)
      //float tmp = (1e3*q[itr]*conversion);
      tr.stas=(short) (1e3*q[itr]*conversion); // moveout in ms
      //fprintf(stderr,"stas=%d, tmp=%f\n",tr.stas,tmp);
      tr.sx=(int) x[ix]; 
      tr.gx=(int) x[ix];
      fputtr(modelfilep,&tr);    
    }
  }

  

  free1float(cdp);
  free2float(ovv);
  free1float(Wm);
  free2float(mm);
  free1float(q);
  free1float(velint);
  free1float(t);
  free1float(h);
  free1float(x);
  free1int(nhcmp);
  free2float(cmp);
  free1float(d);
  free1int(hmute);
  efclose(headerfp);
  efclose(tracefp);
  efclose(modelfilep);


  if (istmpdir) eremove(headerfile);
  if (istmpdir) eremove(tracefile);  

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

int getFirstBreak(float* data, int nt, float threshold, int minIt){
    int it=minIt;while ((it<nt)&&(fabs(data[it])< threshold)) it++;
    return it;
}
