#include "su.h"
#include "radonsolver.h"
#define SIDE 1    // side 1 --> right mute ; side 1 left mute

void radon_param(float fmax, float *x, int nh, float dx, float qmin, float& qmaxt, 
                 float& qmax, float& dq, int nq, int  rtmethod, float factor);
void interval(float *x,int lx,float& mx, float& ax);
void xplotgather(float **d, int nh, int nt, float dt, char *s, char *s2);
void xplotgather(float **d, int nh, float* h, int nt, float dt, char *s, char *s2);
void adaptive_subtract(float **data1, float **data2, int nh, int nt);
void mute(float **model, int nq, int nt, float *q, float *t, float parmute, float t0, int side);

void radonsolver0_subtract(float **data, float *h, int nh,float *t, int nt, float dt, 
                           float **model, float *q, int nq, 
                           float *vel, int itercg, int iter_end, float step, 
                           float eps2, float eps1, float quantil, int norm, 
                           float factor, float smute, float nmofactor, 
                           int rtmethod, float depth, float fmax, char *solver, 
                           float **M, int muteflag){

  // Radon driver function. It takes data with dimensions and calculates model (Radon).
  // it requires in input all user parameters for inversion. 
  // Also performs nmo (before), mute (in model), prediction (inverse RT) and subtraction. 
  // Daniel Trad, UBC, June 2000. 
  
  int i, it, ih;
  float  qmax, qmaxt;
  float  dx_av, dx_max;
  float *dtemp;
  float dq;
  float qmin=q[0];
  float *Wd; 
  float eps=1e-7;
  int testadj=0;
  float **datapred;
  int side=SIDE;   // side 1 --> right mute ; side 1 left mute
  dtemp=ealloc1float(nt);
  Wd=ealloc1float(nh);
  datapred=ealloc2float(nt,nh);
   
  for (ih=0;ih<nh;ih++) Wd[ih]=1;
  if (0){ // recalculate q axis for this cdp offset
      interval(h,nh,dx_max,dx_av);
      fprintf(stderr,"dx_max=%f, dx_av=%f\n", dx_max, dx_av);  
      radon_param(fmax,h,nh,dx_av,qmin,qmaxt,qmax,dq,nq,rtmethod,factor);
      for (i=0;i<nq;i++)  q[i]=qmin+i*dq; 
  }
  else{ // use same q for all cdps, precalculated in initialization
      dq=q[1]-q[0];
      qmax=q[nq-1];
  }

  fprintf(stderr,"freq max=%f,nq=%d dq=%e, nh=%d, nt=%d\n", fmax,nq,dq,nh,nt);
  
  
  if (nmofactor)
  for (ih=0;ih<nh;ih++){
    nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,0,nt,dt,smute);  
    for (it=0;it<nt;it++) data[ih][it]=dtemp[it];
  }

  //  plotgather(data,nh,nt,dt,"suxwigb");
  if (0) plotgather_pipe(data,nh,nt,"After_NMO");
  memset( (void *) model[0], (int) '\0', nq * nt *FSIZE);  
  
  // forward RT using different algorithms (specified in solver)
  radonsolver(h,nh,data,t,nt,dt,model,q,nq,dq,eps1,eps2,eps,fmax,Wd,itercg,iter_end,
	       norm,step,testadj,rtmethod,depth,solver);
  
  /* Any automatic filtering process can be done here on model[nq][nt] */
  if (0){ // set to 1 to plot the model during processing
    xplotgather(model,nq,q,nt,dt,"unmutedmodel","xbox=0 legend=1 perc=98&");
    if (M) xplotgather(M,nq,q,nt,dt,"Mask","xbox=600 legend=1&");
  }
  if (muteflag) AtimesB(model,M,nq,nt);
  if ((0)&&(muteflag)){ // set to 1 to plot model after muting 
    xplotgather(model,nq,q,nt,dt,"mutedmodel","xbox=0 legend=1 perc=98 &");
  }  
  /* compute inverse RT or predicted data */
  hrrti(datapred,h,dt,model,q,fmax,nt,nh,nq,rtmethod,depth);
  /* subtract predicted data */
  /* To get predicted rather than difference set this */
  #define subtract 0
  if ((muteflag)&&(subtract)) adaptive_subtract(data,datapred,nh,nt);
  else{
    for (ih=0;ih<nh;ih++)
      for (it=0;it<nt;it++)
	data[ih][it]=datapred[ih][it];
  }  
  

   
  if (nmofactor)
  for (ih=0;ih<nh;ih++){
    for (it=0;it<nt;it++) dtemp[it]=data[ih][it];    
    nmo(data[ih],dtemp,t,nmofactor*h[ih],vel,1,nt,dt,smute);  
  }    

  free1float(Wd);
  free1float(dtemp);
  free2float(datapred);

  return;

}


void xplotgather(float **d, int nh, int nt, float dt, char *s, char *s2)
{
  char buf[120];
  save_gather(d,nh,nt,dt,s);
  sprintf(buf,"suximage < %s title=%s curve=curve1 npair=5 hbox=900 wbox=700 %s\n",s,s,s2);
  system(buf);
  return;
}

void xplotgather(float **d, int nh, float* h, int nt, float dt, char *s, char *s2)
{
  char buf[120];
  save_gather(d,nh,h,nt,dt,s);
  sprintf(buf,"suximage < %s title=%s curve=curve1 npair=5 hbox=900 wbox=700 %s\n",s,s,s2);
  system(buf);
  return;
}
void interval(float *x,int lx,float& mx, float& ax){
    int i;
    float dx;		
    mx=0;
    ax=0;	
    for (i=1;i<lx;i++){
	dx=fabs(x[i]-x[i-1]);
	if (dx>mx) mx=dx;
	ax=ax+dx;
    }
    ax=ax/(lx-1);
    return;
}

void radon_param_init(float rtmethod, float fmax, float factor,
                      float aperture, float dx, 
                      float moveoutmin, float moveoutmax, float maxoffset,
                      float& dq, int& nq, float& qmin, float& qmax){

    // initial estimates for radon parameters.
    // these parameters will be re-estimated as aperture changes but nq will be fixed     
    if (rtmethod==2 || rtmethod==5) { //PRT
        aperture=aperture*aperture;
        maxoffset=maxoffset*maxoffset;
    }


    dq= factor*fabs(1/(fmax*aperture));          
    qmin=moveoutmin/(maxoffset);
    qmax=moveoutmax/(maxoffset);
    nq=(qmax-qmin)/dq + 1;
    return;
    
}  
    



void radon_param(float fmax, float *x, int nh, float dx, float qmin, float& qmaxt, 
                 float& qmax, float& dq, int nq, int  rtmethod, float factor)
/*

   Given the field geometry (dx, xmin, xmax and fmax) and the chosen qmin,    
   it computes dq, maximum allowable qmax, nq.
    Hence, the NMO must be adjust such that q < qmax.	
    The dx can be the average or maximum.
    rtmethod=1 LRT 
    rtmethod=2 PRT 
    This version allows to define a different factor to undersample the q space
    For example when we want to go beyond the qmax defined by Nyquist
    Daniel Trad- UBC- 16-2-99

*/
{  
   dq=0, qmax=0, qmaxt=0;
   float xmin=x[0];
   float xmax=x[nh-1];
   
   if (rtmethod==2 || rtmethod==5) { //PRT
	  dq= 1/(fmax*(xmax-xmin)*(xmax-xmin));
          dq=factor*fabs(dq); 
          qmax=qmin+dq*(nq-1);
	  qmaxt = 1/(2*fmax*(fabs(xmax-xmin))*dx);
    }
   else if(rtmethod==1 || rtmethod==3) { //LRT
	  dq= 1/(fmax*fabs(xmax-xmin));
          dq=factor*fabs(dq); 
          qmax=qmin+dq*(nq-1);
	  qmaxt = 1/(fmax*dx);  
    }
    return;
}
	
void adaptive_subtract(float **data, float **datapred, int nh, int nt){

  /* Calculate a scale factor that minimize the difference between
     data and datapred in a least squares sense */
  
  float dpd=0;
  float dpdp=0;
  float scale=0;
  int ih,it;

  fprintf(stderr,"substracting multiples \n");
  dpd=dot(nh*nt,data[0],datapred[0]);
  dpdp=dot(nh*nt,datapred[0],datapred[0]);
  scale=dpd/dpdp;
  fprintf(stderr,"scale ===>%f\n",scale);
  for (ih=0;ih<nh;ih++) 
    for (it=0;it<nt;it++)
      data[ih][it]=data[ih][it]-scale*datapred[ih][it];

  return;
}

















