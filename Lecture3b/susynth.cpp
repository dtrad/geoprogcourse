/* Copyright (c) Daniel Trad, University of British Columbia, 1999.*/
/* SUSYNTH:  $Date: December 1999  */

#include "su.h"
#include "segy.h"

#define TRACE { \
  fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__); \
  fflush(stderr); \
}
void savefile(float*d, int n1, const char* name);
int count_ascii_file_lines(const char * name);
int read_ascii_file(const char *name,float *x, int n);
void kolmogoroff(int nt, float *x); 
void conv(int lx, int ifx, float *x, int ly, int ify, float *y,
	  int lz, int ifz, float *z);
void ricker1_wavelet (int nt, float dt, float fpeak, float *wavelet);
/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUSYNTH    > output                                                ",
  " Required parameters: none                                           ",
  "                                                                     ",
  " Generates synthetic data based on the double square root equation   ",
  " The events are first generated as impulses and then convolved with  ",
  " a Ricker wavelet. A minimum phase Ricker is also available by using ",
  " the Kolmogoroff algorithm.                                          ",
  "                                                                     ",
  " Events can be linear, parabolic and hyperbolic.                     ",
  " The curvature or slope of these events are defined by arrays        ",
  " (parameters  vel and tau). Their amplitude are defined by the       ",
  " array coeff.                                                        ",
  "                                                                     ",
  " AVO effects can be included by using the option avo. For example    ",
  " avo=2,3 means that the second and third event will have AVO efects  ",
  " This avo effect is purely geometric and is defined in terms of      ",
  " frequency and phase. avof=0 means maximum aplitude at zero offset   ",
  "                                                                     ",
  " The output offset can be generated automatically or read             ",
  " from an ascii  file                                                 ",
  "                                                                     ",
  "                                                                     ",
  "                                                                     ",
  " Optional parameters                                                 ",
  "  nh      offset traces                                              ", 
  "  nx      Number of midpoints                                        ",
  "  nt          number of time samples                                 ",
  "  dh          offset interval                                        ",
  "  dx          midpoint interval                                      ",
  "  dt          sampling interval                                      ",
  "  hnear   first offset                                               ",
  "  offsetfile  if irregular offset is desired an ascii file should    ",
  "              be given.                                              ",
  "  tau      array of two way travel time                              ",
  "  coeff    array of reflection coefficients                          ",
  "  vel      array RMS velocities (for example use surmsvel)           ",
  "  avo=0    array that defines which events have AVO           	",
  "  avof=0   frequency for AVO  (one value only)                       ", 
  "  avop=0   phase for AVO  (one value only)                           ",
  "  wave=0      0 Ricker 1 MP Ricker                                   ",
  "  fpeak=25    peak frequency for the wavelet                         ", 
  "  nw          number of points for the ricker wavelet                ",
  "  shape=3     1 Line 2 Parabola 3 Hyperbola 				",
  "                                                                     ",
  "  Example 1                                                          ",
  "    susynth nh=30 nx=3 dh=100 dx=250 fpeak=20 nt=1024  \\            ",
  "     fx=0 vel=1500,2000,3000  tau=0.5,1,2 coef=1,-0.5,1 \\           ",
  "     wave=0 avof=1.1 avop=1 avo=2,1   | suxwigb  &                   ",
  "  Example 2                                                          ",
  "     susynth dh=50 nh=63 nx=1 fpeak=20 nt=512 \\                     ",
  "      vel=1700,2000,2500,2900,3000,3500 tau=0.5,0.7,1,1,1.3,1.5 \\   ",
  "      hnear=-1550 | suaddnoise  > data.su &                          ", 
  NULL};

/* Credits:
 *	Daniel Trad. University of British Columbia. 1998
/**************** end self doc ***********************************/

int main(int argc, char **argv) {
    int nt; // trace length
    int ntr; // number of traces;
    int nvel; // number of velocities (layers)
    int ntau; // number of two way travel times (layers)
    int ncoef; // number of reflectivities (layers)
    // notice nvel == ntau == ncoef 
    int nx; // number of midpoints
    int nh; // number of offsets    
    float *coef; // seismic coefficients array 
    float *wavelet = 0; // array to contain the wavelet
    float *vel = 0; // velocity array 
    float *data = 0;
    float *tau = 0; // two way travel time array 
    
    float fpeak; // peak frequency for the Ricker wavelet
    float dt;
    float offset;
    float dh; // offset sampling 
    float dx; // midpoint sampling 
    float tau0;
    float hnear; // minimum offset 
    float fx;
    float time;
    float veloc;
    float hoff;    
    
    int wave; // choice for wavelet
    float avof; // frequency of AVO 
    float avop; // AVO phase
    int *avo; // array containing which events have avo effect 
    int navo;
    float amplitude;
    float avoscale;
    int nw;
    int shape; // 1 Line 2 Parabola 3  Hyperbola
    float midxs;
    float *h;
    cwp_String offsetfile = NULL; /*input ascii file for offset if interpolation is desired */
    segy tr;

    // Initialize 

    initargs(argc, argv);
    requestdoc(0);
    fprintf(stderr, "argc=%d\n", argc);

    if (!getparfloat("fpeak", &fpeak)) fpeak = 25;
    if (!getparint("nt", &nt)) nt = 512;
    if (!getparint("nh", &nh)) nh = 50;
    if (!getparint("nx", &nx)) nx = 1;
    if (!getparfloat("dt", &dt)) dt = 0.004;
    if (!getparfloat("dh", &dh)) dh = 25;
    if (!getparfloat("dx", &dx)) dx = 250;
    if (!getparfloat("hnear", &hnear)) hnear = -((nh - 1) * dh / 2.);
    if (!getparfloat("fx", &fx)) fx = 0;
    if (!getparint("wave", &wave)) wave = 1;
    if (!getparint("nw", &nw)) nw = 50;
    if (!getparfloat("avof", &avof)) avof = 0;
    if (!getparfloat("avop", &avop)) avop = 0;
    if (!getparint("shape", &shape)) shape = 3;
    if (!getparstring("offsetfile", &offsetfile)) offsetfile = NULL;

    /* If offsetfile name is given read it */
    if (offsetfile) {
        nh = count_ascii_file_lines(offsetfile);
        h = ealloc1float(nh);
        read_ascii_file(offsetfile, h,nh);
    }
    else{
        h = ealloc1float(nh);
        if (!h) perror("can't allocate memory for offset");
    }
    ntr = nx*nh;

    //////////////////////////////////////////////////////////////////
    nvel = countnparval(1, "vel");
    if (nvel == 0) nvel = 1;
    if ((vel = ealloc1float(nvel)) == NULL)
        fprintf(stderr, "***Space for vel could not be allocated\n");
    if (!getnparfloat(1, "vel", vel)) vel[0] = 1500;
    /////////////////////////////////////////////////////////////////        
    ntau = countnparval(1, "tau");
    if (ntau == 0) ntau = 1;
    if ((tau = ealloc1float(ntau)) == NULL)
        fprintf(stderr, "***Space for vel could not be allocated\n");
    if (!getnparfloat(1, "tau", tau)) tau[0] = 0.5;
    if (ntau != nvel) err("nvel and ntau must be equal\n");
    //////////////////////////////////////////////////////////////////
    ncoef = countnparval(1, "coef");
    if (ncoef == 0) ncoef = 1;
    if (ncoef != nvel) warn("***, ncoef differs from nvel. used coeff=1 instead");
    coef = ealloc1float(nvel);

    if (!getnparfloat(1, "coef", coef) || ncoef != nvel)
        for (int iv = 0; iv < nvel; iv++) coef[iv] = 1;


    //////////////////////////////////////////////////////////////////
    navo = countnparval(1, "avo");
    avo = ealloc1int(navo);

    if (!getnparint(1, "avo", avo))
        for (int iv = 0; iv < navo; iv++) avo[iv] = 0;
    //////////////////////////////////////////////////////////////////

    for (int it = 0; it < nvel; it++) fprintf(stderr, "vel[%d]=%f,tau[%d]=%f\n", it, vel[it], it, tau[it]);

    for (int iavo = 0; iavo < navo; iavo++) fprintf(stderr, "avo[%d]=%d\n", iavo, avo[iavo]);

    fprintf(stderr, "nt=%d,dt=%f;nw=%d,ntr=%d,fx=%f,navo=%d,vel=%d\n", nt, dt, nw, ntr, fx, navo, nvel);
    wavelet = ealloc1float(nw);
    data = ealloc1float(nt);


    if (wave == 0 || wave == 1) ricker1_wavelet(nw, dt, fpeak, wavelet);
    if (wave == 1) kolmogoroff(nw, wavelet);
    //find wavelet peak and index
    int iwmax = 0;
    if (1) {
        float amax = 0;
        for (int it = 0; it < nw; it++) {
            if (fabs(wavelet[it]) > amax) {
                amax = wavelet[it];
                iwmax = it;
            }
        }
    }


    savefile(wavelet, nw, "wavelet.dat");   
    fprintf(stderr, "**************************\n");
    for (int ix = 0; ix < nx; ix++) { // loop over midpoints
        midxs = fx + ix*dx; // regular midpoint axis 
        for (int ih = 0; ih < nh; ih++) { // Loop over offsets           
            if (offsetfile) offset = h[ih]; // irregular offset from file.
            else offset = ih * dh + hnear;  // regular offset created on the fly
            hoff = offset / 2.;
            memset(data,0,sizeof(float)*nt);
            for (int iv = 0; iv < nvel; iv++) { // loop over layers                
                veloc = vel[iv];
                tau0 = tau[iv];
                /* Use the DSR to compute travel time */
                if (shape == 3)
                    time = sqrt(tau0 * tau0 / 4 + ((midxs + hoff)*(midxs + hoff))
                        / (veloc * veloc)) + sqrt(tau0 * tau0 / 4 +
                        ((midxs - hoff)*(midxs - hoff)) / (veloc * veloc));
                else if (shape == 2)
                    time = (tau0 / 2 + ((midxs + hoff)*(midxs + hoff))
                        / (veloc * veloc))+(tau0 / 2 +
                        ((midxs - hoff)*(midxs - hoff)) / (veloc * veloc));
                else if (shape == 1)
                    time = tau0 + (fabs(2 * hoff) / veloc);
                // assign sample by linear interpolation.
                // for more precision this can be replaced by ints8r (sinc interpolation)
                float itfloat = time / dt;
                int itint = (int) floor(itfloat + 0.5);
                float a1 = 1 - (itfloat - itint);
                float a2 = itfloat - itint;
                // arbitrary AVO function
                amplitude = coef[iv];
                for (int iavo = 0; iavo < navo; iavo++){
                    if (avo[iavo] == (iv + 1)) {
                        avoscale = cos((avof / 1000.) * offset - (avop / 10.));
                        amplitude = coef[iv] * avoscale;
                    }
                }
                if (itint < nt) data[itint] = data[itint] + a1 * amplitude;
                if ((itint + 1) < nt) data[itint + 1] = data[itint + 1] + a2 * amplitude;
            }
            

            if (wave == -1) // output spikes without convolution
                for (int it = 0; it < nt; it++) tr.data[it] = data[it];
            else
                conv(nt, -iwmax, data, nw, 0, wavelet, nt, 0, tr.data);

            // assign header words
            tr.cdp = (int) (midxs + 1);
            tr.offset = (int) offset;
            tr.ns = (unsigned short) nt;
            tr.dt = (unsigned short) (dt * 1e6); //microsec
            tr.ntr = (int) ntr;
            puttr(&tr); // output trace
        }
    }

    free1float(data);
    free1float(wavelet);
    free1int(avo);
    free1float(coef);
    free1float(tau);
    free1float(vel);
    free1float(h);
    return EXIT_SUCCESS;
}
// auxiliary functions
int count_ascii_file_lines(const char * name){
    int nlines=0;
    float temp;
    FILE* fp = efopen(name,"r");
    int nn;
    do{
        nn = fscanf(fp,"%f",&temp);
        fprintf(stderr,"%f",temp);
        nlines++;
    } while (nn==1);
    nlines=nlines-1;
    fclose(fp);
    fprintf(stderr,"nlines = %d\n",nlines);
    return nlines;
}
int read_ascii_file(const char *name,float *x, int nx)
{
  FILE *fp;
  fp=efopen(name,"r");
  for (int ix=0;ix<nx;ix++){
    fscanf(fp,"%f",&x[ix]); 
    fprintf(stderr,"h[%d]=%f\n",ix,x[ix]);
  }     
  efclose(fp);
  return(nx);
}
void kolmogoroff(int nt, float *x) 
{
  /* 
     Given a non minimum phase wavelet it returns its minimum phase version 
     The program kolmogoroff() first takes the logarithm of the spectrum, 
     then returns to the time domain and sets to zero the noncausal part. 
     It returns to frequency, exponentiates, and returns to the time domain
     with a wavelet that will be proven to be minimum-phase.  
     See Claerbout, 1992. Processing versus inversion for the description. 
  */
    
  int it;
  int nfft;
  complex czero;
  complex *cx;

  czero.r=czero.i=0;
  nfft=npfa(nt);
  cx=ealloc1complex(nfft);
  
  /* store the wavelet in a complex array */
  for (it=0;it<nt;it++) {cx[it].r=x[it];cx[it].i=0;}
  for (it=nt;it<nfft;it++) cx[it]=czero;
  /* Compute spectrum */
  pfacc(1,nfft,cx);  

  for (it=0;it<nfft;it++) cx[it]*=sqrt(1./nfft);
  for (it=0;it<nfft;it++) cx[it]=cx[it]*conjg(cx[it]); 
  /* Take log */
  for (it=0;it<nfft;it++) cx[it]=log(cx[it]);
  pfacc(-1,nfft,cx);  
  for (it=0;it<nfft;it++) cx[it]=cx[it]*sqrt(1./nfft);
  cx[0]/=2.;
  cx[nfft/2]/=2.;

  for (it=nfft/2+1;it<nfft;it++) cx[it]=czero;
  pfacc(1,nfft,cx);  for (it=0;it<nfft;it++) cx[it]=cx[it]*sqrt(1./nfft);
  for (it=0;it<nfft;it++)  cx[it]=exp(cx[it]);
  pfacc(-1,nfft,cx);  for (it=0;it<nfft;it++) cx[it]=cx[it]*sqrt(1./nfft);

  for (it=0;it<nfft/2;it++) x[it]=cx[nfft/2+it].r; 

  free1complex(cx);

  return;
}

void savefile(float*d, int n1, const char* name){
  FILE* fp =0;
      
  if ((fp=fopen(name,"w"))==NULL){
    fprintf(stderr,"cant open file  \n");         
  }
  for (int i=0;i<n1;i++) fprintf(fp,"%f\n",d[i]); // ascii
  //fwrite(d,sizeof(float),n1,fp);                // binary
  fclose(fp);
}












