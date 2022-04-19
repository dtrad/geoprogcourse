#include "su.h"
#include "radontd_win2.h"

#define NWIN 1 // Number of data space windows


void radontd_sparse_win(float *t, float *q, float *h, float **m,float **d, int nt, 
                        int nh, int nq, float dt, float **vgrid, inv_par inv, 
                        int dataprec, int nw, float fpeak, int typewav, int LI, 
                        float parmute, int mute, int itm, int plot);

void radontd_sparse(float *t, float *q, float *h, float **m,float **d, int nt, 
                    int nh, int nq, float dt, float *vel, float dperv, float pervmin, 
                    float t0,  inv_par inv, int centralq, int dataprec, 
                    int nw, float fpeak, int typewav, int LI, 
                    float parmute, int mute, int itm, int plot)

  /*
    RADONTD_SPARSE
    RADON TRANSFORM IN TIME DOMAIN: with linear interpolation and wavelet convolution
    This function is called by suradontd1.cpp and suradonlinetd.cpp 
    This routine is just an interface to radontd_sparse_win, that performs the RT 
    The interface takes windows in the  data space and calls the RT

    Daniel Trad - UBC October 2000
    
    E-mail: dtrad@geop.ubc.ca
  */
{
  register int it;
  int  ih, iq, it0;
  /****************************Model window********************************/ 
  int nwin=NWIN;
  int it0_total=(int) (t0/dt+0.5);
  int ntwin=(nt-it0_total)/nwin;
  float **mw;     // window for the model 
  float **dw;     // window for the model
  float **vgridw; // window for the model
  int iwin;
  itm=itm-it0_total; if (itm<0) itm=0;
  /***********************************************************************************/
  //   Velocity axis //////////////////////////////////////////////////
  /* The grid of square slowness is stored in the array vgrid and this is written 
     into  a binary file vgrid */
  float **vgrid;

  vgrid=alloc2float(nt,nq);
  vgridw=alloc2float(ntwin,nq);
  mw=alloc2float(ntwin,nq);
  dw=alloc2float(ntwin,nh);

  fprintf(stderr,"t0=%f,it0_total=%d,dt=%f,t[%d]=%f\n",t0,it0_total,dt,it0_total,t[it0_total]);
  fprintf(stderr,"nq=%d,nt=%d,pervmin=%f,dperv=%f,centralq=%d\n",
	  nq,nt,pervmin,dperv,centralq);

 
  irreg_slowness_axis(nq,nt,pervmin,dperv,t,vel,q,vgrid,centralq);
  save2dfilepower(vgrid,nq,nt,dt,-0.5,"vgrid");

  // Loop for the data windows
  for (iwin=0;iwin<nwin;iwin++){  
    
    it0=iwin*ntwin+it0_total; 
    for (ih=0;ih<nh;ih++) for (it=0;it< ntwin; it++) dw[ih][it]=d[ih][it+it0];

    if (0){
      save_gather(dw,nh,ntwin,dt,"dw");
      system("suxwigb < dw  title=\"plot dw\" ");
    }

    for (iq=0;iq<nq;iq++) for (it=0;it< ntwin; it++) vgridw[iq][it]=vgrid[iq][it+it0];  
    //memcpy((void *) dw[ih],(const void *) &d[ih][it],ntwin*sizeof(float));

    radontd_sparse_win(&t[it0],q,h,mw,dw,ntwin,nh,nq,dt,vgridw,inv,dataprec,nw,fpeak,typewav,LI,parmute,mute,itm,plot);   
    
    if ((plot==1)||(plot==3)){
      save_gather(mw,nq,q,ntwin,dt,"mw");
      system("suxwigb < mw perc=97 key=f2 title=\"RT_in_radontd_win\" ");
    }
    if ((plot==1)||(plot==3)){
      save_gather(dw,nh,h,ntwin,dt,"dw");
      system("suxwigb < dw perc=97 key=offset  title=\"primaries\" ");
    }
    for (iq=0;iq<nq;iq++) for (it=0;it< ntwin; it++) m[iq][it+it0]=mw[iq][it];
    for (ih=0;ih<nh;ih++) for (it=0;it< ntwin; it++) d[ih][it+it0]=dw[ih][it];
  }

  free2float(mw);
  free2float(dw);
  free2float(vgrid);
  free2float(vgridw);
  return;
  
}

void radontd_sparse_win(float *t, float *q, float *h, float **m,float **d, int nt, int nh, int nq, float dt, float **vgrid, inv_par inv, int dataprec, int nw, float fpeak, int typewav, int LI, float parmute, int mute, int itm, int plot)

  /*
    RADONTD_SPARSE_WIN
    RADON TRANSFORM INSIDE A WINDOW IN THE TIME DOMAIN

    Daniel Trad - UBC October 2000
    E-mail: dtrad@geop.ubc.ca
  */
{
  register int it;
  int  ih, iq;
  float **Wm=0; // Model weights
  float **Wd=0;// Data weights
  float **dtemp=0;// Data weights
  unsigned int **index=0; // sparse matrix L
  int nsparse;
  int ntaper=5; // Number of lateral traces to taper 
  float test;
  float *wavelet=0;
  int testadj=0; // set to 1 to test operators. 
  int root_nqxnh=(int) sqrt(nq*nh);
  //if (inv.iter_end==1) inv.restart=0;
  float sigmam;
  float sigmad; 
  int iter;
  //////////////////////////////////////////////////////////////////////////
  // Define pointers to functions to be used for operators
  // hyperbolic RT 
  void (*radon3s) (float *m, float *d, unsigned int **index, int adj, int nt, 
		   int nh, int nq, int nsparse);
  
  // hyperbolic RT with wavelet convolution
  void (*radon3sw) (float *m, float *d, unsigned int **index, int adj, int nt, 
		   int nh, int nq, int nsparse, float *wavelet, int nw);
  
  // convolution operator
  void (*convop)  (int, int, int, float *,int, float *,float *);

  // Define the actual functions to be used 
  if (!LI){ // nearest neighbor (binning)
    radon3sw=radonhyp;
    radon3s=radonhyp;
  }
  else{     // linear interpolation
    radon3sw=radonhyp_li;
    radon3s=radonhyp_li;
  }
  convop=contruc_2;
  // Generate a wavelet for the forward and adjoint operator and plot it
  if (nw){    
    //typewav=1==>   Read a wavelet generated by Matlab (stored in an ascii file)
    //typewav=2==>   Generate a Ricker with SU
    //typewav=3==>   Design a simple wavelet by hand    
    wavelet=ealloc1float(nw);
    nw=get_wavelet(wavelet,"wavelet.txt",nw,typewav,dt,fpeak);
    if (!nw) abort();
    wavelet=erealloc1float(wavelet,nw);
    fprintf(stderr,"Test dot product for contran\n");
    test=testadjop_conv(convop,nw,wavelet,nt,nt+nw-1);
  }    

  Wm=ealloc2float(nt,nq);
  Wd=ealloc2float(nt,nh);
  dtemp=ealloc2float(nt,nh);

  if (0) if (inv.taperflag==1) taper(d,nt,nh,ntaper,0);
  /***************************************************************************/  

  for (iq=0;iq<nq;iq++) memset((void *)m[iq],(int)'\0',nt*FSIZE);
  
  nsparse=nt*nq*nh; 
  fprintf(stderr,"nsparse=%d,t[0]=%f,nt=%d\n",nsparse,t[0],nt);
  // allocate the big monster
  size_t size=sizeof(unsigned int);
  if ((index=(unsigned int **) alloc2(nsparse,2,size))==NULL) 
    err("Cannot allocate index\n");

  // Assign the elements of index
  if (!LI) build_index_slowness(t,h,q,vgrid,nt,nh,nq,index);
  else build_index_slowness_li(t,h,q,vgrid,nt,nh,nq,index);

  // Test the adjoint
  if (testadj && nw ) test=testadjop(radon3sw,index,nt,nh,nq,nsparse,wavelet,nw);
  else if (testadj) test=testadjop(radon3s,index,nt,nh,nq,nsparse);
  /***************************************************************************/

  fprintf(stderr,"Wd = 1 ............\n");    
  for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) Wd[ih][it]=1.0;

  // This is a filter for outliers or dominant bad data
  if (dataprec) filterOutlier(dataprec,h,nh,t,nt,dt,d,Wd); 

  /***************************************************************************/
  // Compute the adjoint before the first iteration to obtain an initial model   
  for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) dtemp[ih][it]=d[ih][it]*Wd[ih][it];

  if (nw) radon3sw(m[0],dtemp[0],index,1,nt,nh,nq,nsparse,wavelet,nw);
  else if (inv.restart==1) radon3s(m[0],dtemp[0],index,1,nt,nh,nq,nsparse);
  else radonhyp_sinc(m[0],t,h,q,dtemp[0],vgrid,1,nt,nh,nq); 

  // Scale the adjoint to get an initial model 
  for (iq=0;iq<nq;iq++) for (it=0;it<nt;it++) m[iq][it]/=root_nqxnh;  

  /***************************************************************************/

  fprintf(stderr,"inv.norm=%d,inv.eps1=%f,inv.eps2=%f,inv.itercg=%d,inv.iter_end=%d,"
            "inv.eps=%f,inv.restart=%d,inv.step=%f\n",inv.norm,inv.eps1,inv.eps2,
            inv.itercg,inv.iter_end,inv.eps,inv.restart,inv.step);    

  // IRLS loop. Modelweight computes the model preconditioner
  for (iter=1;iter<=inv.iter_end;iter++){
    if (iter>1) inv.restart=1;
    
    // norm==1 ==> L1 , ==0  Cauchy else l2
    if (iter>1) 
      deviations(m[0],nq*nt,d[0],nh*nt,inv.norm,inv.eps1,inv.eps2,&sigmam,&sigmad);
    
    weights_td(m[0],nq*nt,inv.norm,sigmam,Wm[0],iter);
    if (iter==inv.iter_end) inv.itercg*=1;  // Last external iteration is 3 times longer 
    if (nw){ 
      fprintf(stderr,"iteration with wavelet convolution\n");
      wpcgnr(radon3sw,nt,nh,nq,nsparse,m[0],d[0],Wd[0],Wm[0],index,inv,wavelet,nw);
    }
    else wpcgnr(radon3s,nt,nh,nq,nsparse,m[0],d[0],Wd[0],Wm[0],index,inv);
  }

  // either predict data directly or subtract 
   
  if (mute){  
    if (nw) subtract_multiples(radon3sw,m,d,index,0,nt,nh,nq,nsparse,wavelet,nw,parmute,q,h,itm,plot);
    else subtract_multiples(radon3s,m,d,index,0,nt,nh,nq,nsparse,parmute,q,h,itm,plot);
  }
  else{
    if (nw) radon3sw(m[0],d[0],index,0,nt,nh,nq,nsparse,wavelet,nw);
    else radon3s(m[0],d[0],index,0,nt,nh,nq,nsparse);
  }
  //float scalefactor=1;//((float)nh)/nq;
  //for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) d[ih][it]*=scalefactor;
  
  free2((void **) index);
  free2float(Wm);
  free2float(Wd);  
  free2float(dtemp);
  if (nw) free1float(wavelet);

  return;
  
}
void irreg_slowness_axis(int nq, int nt, float pervmin, float dperv, float *t, float *vel,float *q, float **vgrid, int centralq)
{
  /*   Velocity Grid 
   Irregular velocity grid
   Define irregular spacing along the horizontal
   Here q  is perturbartion around the central velocity law
   dperq is a parameter used to generate the increasing space
   perqmin defines the minimum distance to the perturbation */
     
  float *perv;
  int it;
  int iq;
  int nqh=centralq;
  float vaux;
  fprintf(stderr,"pervmin=%g,dperv=%g\n",pervmin,dperv);
  if ((perv=alloc1float(nq+1))==NULL) err(" Cannot allocate perv");
  perv[nqh]=0;
  for (iq=nqh;iq<nq;iq++) perv[iq+1]=perv[iq]*(1+dperv)+pervmin;
  for (iq=nqh-1;iq>=0;iq--) perv[iq]=perv[iq+1]*(1+dperv)-pervmin;
  for (iq=0;iq<nq;iq++){
    //fprintf(stderr,"perv[%d]=%f\n",iq,perv[iq]);
    for (it=0;it<nt;it++){
      vaux=perv[iq];
      vgrid[iq][it]=1./(vel[it]*vel[it])+vaux;
      //if (vgrid[iq][it]<0) vgrid[iq][it]=0;//vgrid[iq][MAX(it-1,0)];
    }
    q[iq]=(vgrid[iq][0]);
    fprintf(stderr,"vtop[%d]=%6.0f<======>,vbot[%d]=%6.0f\n",iq,sqrt(1./q[iq]),iq,sqrt(1./vgrid[iq][nt-1]));
  }
  //  for(iq=0;iq<nq;iq++) fprintf(stderr,"q[%d]=%f\n",iq,q[iq]);
  free1float(perv); 
  return;
}
void build_index_slowness(float *t, float *h, float *q, float **vel, int nt, int nh, int nq,unsigned int **index)
{
  register int it;
  int ih,iq;
  float time,hxh,pxhxh;
  int iqxnt,ihxnt;
  int itime;
  int nsparse=nt*nq*nh;
  float dt=t[1]-t[0];

  for (it=0;it<nsparse;it++) index[0][it]=index[1][it]=0;
  
  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){    
      iqxnt=iq*nt;
      for (it=0;it<nt;it++){
	if (vel[iq][it]>0){
	  pxhxh=hxh*vel[iq][it];
	  time=sqrt(t[it]*t[it]+pxhxh);
	  itime=(int) ((time-t[0])/dt+0.5);
	  if (itime<nt){
	    index[0][ih*nq*nt+iq*nt+it]=ihxnt+itime;
	    index[1][ih*nq*nt+iq*nt+it]=iqxnt+it;
	  }
	  else{
	    index[0][ih*nq*nt+iq*nt+it]=0;
	    index[1][ih*nq*nt+iq*nt+it]=0;
	  }
	}
      }            
    }
  }
  return;
}
void build_index_slowness_li(float *t, float *h, float *q, float **vel, int nt, int nh, int nq,unsigned int **index)
{
  register int it;
  int ih,iq;
  float time,hxh,pxhxh;
  int iqxnt,ihxnt;
  int itime;
  int nsparse=nt*nq*nh;
  float dt=t[1]-t[0];


  for (it=0;it<nsparse;it++) index[0][it]=index[1][it]=0;
  
  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){    
      iqxnt=iq*nt;
      for (it=0;it<nt;it++){
	if (vel[iq][it]>0){
	  pxhxh=hxh*vel[iq][it];
	  // In this function we multiply by 10 before truncation to save
	  // the first decimal digit. Later the times are multiply by 0.1
	  time=sqrt(t[it]*t[it]+pxhxh);
	  itime=(int) ((10*(time-t[0]))/dt+0.5);

	  if (itime<((nt-1)*10)){
	    index[0][ih*nq*nt+iq*nt+it]=10*ihxnt+itime;
	    index[1][ih*nq*nt+iq*nt+it]=iqxnt+it;
	  }
	  else{
	    index[0][ih*nq*nt+iq*nt+it]=0;
	    index[1][ih*nq*nt+iq*nt+it]=0;
	  }
	}
      }            
    }
  }
  return;
}

/****************************************************************************/


void radonhyp_sinc(float *m,float *t, float *h, float *q, float *d, float **vel,int adj,int nt, int nh, int nq)
{
  register  int it;
  int ih,iq;
  float *ttn,*dint,*tnt,hxh,pxhxh;
  int iqxnt,ihxnt;
  int nx=nt*nq;
  int ny=nt*nh;
  float dt=t[1]-t[0];
  float dt2=dt*dt;



  if ((ttn=alloc1float(nt))==NULL)
    err("cannot allocate memory for ttn \n");
  if ((tnt=alloc1float(nt))==NULL)
    err("cannot allocate memory for tnt \n");
  if ((dint=alloc1float(nt))==NULL)
    err("cannot allocate memory for dint \n");
  fprintf(stderr,"Sinc interpolation for the adjoint +++++++++\n");

  if (adj) memset((void *) m,(int)'\0',nx*FSIZE);// for (it=0;it<nx;it++) m[it]=0;
  else memset((void *) d,(int)'\0',ny*FSIZE);// for (it=0;it<ny;it++) d[it]=0;   

  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){    
      iqxnt=iq*nt;
      for (it=0;it<nt;it++){
	pxhxh=hxh*vel[iq][it];
	ttn[it]=sqrt(t[it]*t[it]/dt2+pxhxh/dt2);
      }	
      if (adj) ints8r(nt,1.0,0,&d[ihxnt],0.0,0.0,nt,ttn,dint);
      else{
	yxtoxy(nt,1.0,0.0,ttn,nt,1.0,0.0,-nt,nt,tnt);
	ints8r(nt,1.0,0,&m[iqxnt],0.0,0.0,nt,tnt,dint);
      }
      if (adj) for (it=0;it<nt;it++) m[iqxnt+it]+=dint[it];
      else  for (it=0;it<nt;it++) d[ihxnt+it]+=dint[it];
    }            
  }
  free1float(dint);
  free1float(ttn);
  free1float(tnt);  
  return;
}

void radonhyp(float *m, float *d, unsigned int **index, int adj, int nt, int nh, int nq, int nsparse)
{
  int j;
  int ny=nh*nt;
  int nx=nq*nt;

  d[0]=0;
  m[0]=0;
  
  if (!adj){
    memset((void *) d,(int)'\0',ny*FSIZE);
    for (j=0;j<nsparse;j++) d[index[0][j]]+=m[index[1][j]];
  }
  else{
    memset((void *) m,(int)'\0',nx*FSIZE);
    for (j=0;j<nsparse;j++) m[index[1][j]]+=d[index[0][j]];
  }
  
  /* 
    A problem appears if some of the values of index are never computed
    because the zero index of d and m are mapped each other for index=0
    I make these two elements equal to zero just to prevent this problem, 
    It does not affect the data or model significantly.  
  */ 
			      
  d[0]=0;
  m[0]=0;

  return;
}


void radonhyp(float *m, float *d, unsigned int **index, int adj, int nt, int nh, int nq, int nsparse, float *wavelet, int nw)
{
  int j,ih;
  int ny=nh*nt;
  int nx=nq*nt;
  float *dtemp;
  float *dtemp2;

  //fprintf(stderr,"nw=%d,nt=%d,nh=%d,nq=%d\n",nw,nt,nh,nq);
  d[0]=0;
  m[0]=0;
  
  dtemp2=ealloc1float(ny);
  dtemp=ealloc1float(nt+nw);
  memset((void *) dtemp,(int)'\0',(nt+nw)*FSIZE); 
  memset((void *) dtemp2,(int)'\0',ny*FSIZE); 

  //for (j=0;j<nw;j++) fprintf(stderr,"nw=%d, wavelet[%d]=%f\n",nw,j,wavelet[j]);
  
  if (!adj){
    memset((void *) d,(int)'\0',ny*FSIZE);
    for (j=0;j<nsparse;j++) d[index[0][j]]+=m[index[1][j]];
    d[0]=0;
    for (ih=0;ih<nh;ih++){
      memset((void *) dtemp,(int)'\0',(nt+nw)*FSIZE); 
      contruc_2(0,0,nw,wavelet,nt,&d[ih*nt],dtemp);
      memcpy((void *) &d[ih*nt],(const void *) dtemp,nt*sizeof(float));
    }
  }
  if (adj){
    for (ih=0;ih<nh;ih++){
      memset((void *) dtemp,(int)'\0',(nt+nw)*FSIZE); 
      contruc_2(1,0,nw,wavelet,nt,dtemp,&d[ih*nt]);
      memcpy((void *) &dtemp2[ih*nt],(const void *) dtemp,nt*sizeof(float));
    }
    memset((void *) m,(int)'\0',nx*FSIZE);
    for (j=0;j<nsparse;j++) m[index[1][j]]+=dtemp2[index[0][j]];
    m[0]=0;
  }  
  
  /* 
    A problem appears if some of the values of index are never computed
    because the zero index of d and m are mapped each other for index=0
    I make these two elements equal to zero just to prevent this problem, 
    It does not affect the data or model significantly.  
  */
  d[0]=0;
  m[0]=0;

  free1float(dtemp2);
  free1float(dtemp);
  
  return;
}
void radonhyp_li(float *m, float *d, unsigned int **index, int adj, int nt, int nh, int nq, int nsparse)
{
  // This function is similar to radonhyp but it performs linear interpolation
  // The first column of index has been multiplied by 10 to keep the first decimal
  // place after truncation. 
  /// D. Trad - October 2000
  int j;
  int ny=nh*nt;
  int nx=nq*nt;
  int id1;
  int im1;
  float a;
  float b;

  d[0]=0;
  m[0]=0;
  
  if (!adj){
    memset((void *) d,(int)'\0',ny*FSIZE);
    for (j=0;j<nsparse;j++){
      b=0.1*index[0][j];
      id1=(int) b;
      im1=index[1][j];
      if (id1 && im1 ){
	a=b-id1;
	d[id1]+=(1.0-a)*m[im1];
	d[id1+1]+=a*m[im1];
      }
    }
  }
  else{
    memset((void *) m,(int)'\0',nx*FSIZE);
    for (j=0;j<nsparse;j++){
      b=0.1*index[0][j];
      id1=(int) b;
      im1=index[1][j];
      if (id1 && im1 ){
	a=b-id1;
	m[im1]+=(1-a)*d[id1]+a*d[id1+1];
      }
    }
  }

  /* 
    A problem appears if some of the values of index are never computed
    because the zero index of d and m are mapped each other for index=0
    I make these two elements equal to zero just to prevent this problem, 
    It does not affect the data or model significantly.  
  */ 
			      
  d[0]=0;
  m[0]=0;

  return;
}



void radonhyp_li(float *m, float *d, unsigned int **index, int adj, int nt, int nh, int nq, int nsparse, float *wavelet, int nw)
{
  int j,ih;
  int ny=nh*nt;
  int nx=nq*nt;
  float *dtemp;
  float *dtemp2;
  int im1;
  int id1;
  float a;
  float b;
  //fprintf(stderr,"nw=%d,nt=%d,nh=%d,nq=%d\n",nw,nt,nh,nq);
  d[0]=0;
  m[0]=0;
  
  dtemp2=ealloc1float(ny);
  dtemp=ealloc1float(nt+nw);
  memset((void *) dtemp,(int)'\0',(nt+nw)*FSIZE); 
  memset((void *) dtemp2,(int)'\0',ny*FSIZE); 

  //for (j=0;j<nw;j++) fprintf(stderr,"nw=%d, wavelet[%d]=%f\n",nw,j,wavelet[j]);
  
  if (!adj){
    memset((void *) d,(int)'\0',ny*FSIZE);
    for (j=0;j<nsparse;j++){
      b= ( 0.1 * index[0][j] );
      id1=(int) b;
      im1=index[1][j];
      if (id1 && im1 ){
	a=b-id1;
	d[id1]+=(1-a)*m[im1];
	d[id1+1]+=a*m[im1];
      }
    }
    d[0]=0;
    for (ih=0;ih<nh;ih++){
      memset((void *) dtemp,(int)'\0',(nt+nw)*FSIZE); 
      contruc_2(0,0,nw,wavelet,nt,&d[ih*nt],dtemp);
      memcpy((void *) &d[ih*nt],(const void *) dtemp,nt*sizeof(float));
    }
  }
  if (adj){
    for (ih=0;ih<nh;ih++){
      memset((void *) dtemp,(int)'\0',(nt+nw)*FSIZE); 
      contruc_2(1,0,nw,wavelet,nt,dtemp,&d[ih*nt]);
      memcpy((void *) &dtemp2[ih*nt],(const void *) dtemp,nt*sizeof(float));
    }
    memset((void *) m,(int)'\0',nx*FSIZE);
    for (j=0;j<nsparse;j++){
      b=( 0.1 *index[0][j] );
      id1=(int) b;
      im1=index[1][j];
      if (id1 && im1 ){
	a=b-id1;
	m[im1]+=(1-a)*dtemp2[id1]+a*dtemp2[id1+1];
      }
    }
    m[0]=0;
  }  
  
  /* 
     A problem appears if some of the values of index are never computed
     because the zero index of d and m are mapped each other for index=0
     I make these two elements equal to zero just to prevent this problem, 
     It does not affect the data or model significantly.  
  */
  d[0]=0;
  m[0]=0;
  
  free1float(dtemp2);
  free1float(dtemp);
  
  return;
}





void weights_td(float *m, int nx, int norm, float sigmam, float *Wm, int iter)
  /*
  The right Wm from Cauchy is 
  Wm[i]=sqrt(eps1*eps1+m[i]*m[i]/(maxm*maxm));
  But if M^-1 ATA x = M^-1 AT b is solved instead
  of the standard form M=WmT *Wm 
  Actually it works even better with (I don't know why)
  Wm[i]=Wm[i]*Wm[i];
  if (Wm[i]>2) Wm[i]=2; 
  */
{ 
      int i;
      
      if (iter==1){ 
	for (i=0;i<nx;i++) Wm[i]=1;
	return;
      }
      
      if (norm==1) for (i=0;i<nx;i++) Wm[i]=(fabs(m[i])*sigmam);
      else if(norm==0)
	for (i=0;i<nx;i++) Wm[i]=(sigmam*sigmam+fabs(m[i]*m[i]));
      else if(norm==2){ // Mask
	for (i=0;i<nx;i++)
	  Wm[i]=1+200./(1+exp(1*(fabs(m[i])-sigmam)+0.5));
      }
      return;
}





float testadjop(void (*oper) (float *,float *,unsigned int **,int ,int ,int, int, int, 
			      float *wavelet, int nw),unsigned int **index,int nt, int nh, 
		int nq, int nsparse, float *wavelet, int nw)
{
  float *dr1;
  float *mr1;
  float *dr2;
  float *mr2;
  float dp1;
  float dp2;
  int it;
  int iq;
  int ih;
  float test;
  int ny=nt*nh;
  int nx=nt*nq;
  //////////////  


  if ((dr1=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr1=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");
  if ((dr2=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr2=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");

 
  for (it=0;it<nt;it++) for (ih=0;ih<nh;ih++) dr1[ih*nt+it]=frannor();
  for (it=0;it<nt;it++) for (iq=0;iq<nq;iq++) mr1[iq*nt+it]=frannor();

  oper(mr2,dr1,index,1,nt,nh,nq,nsparse,wavelet,nw);
  oper(mr1,dr2,index,0,nt,nh,nq,nsparse,wavelet,nw);

  dp1=dot(ny,dr1,dr2);
  dp2=dot(nx,mr1,mr2);

  if (dp2!=0) test=dp1/dp2;
  else test=0;

  fprintf(stderr,"Test adjoint = %f dp1=%f, dp2=%f \n",test,dp1,dp2);
  return(test);
  
  free1float(mr2);
  free1float(dr2);
  free1float(mr1);
  free1float(dr1);

}

float testadjop(void (*oper) (float *,float *,unsigned int **,int ,int ,int, int, int),unsigned int **index,int nt, int nh, int nq, int nsparse)
{
  float *dr1;
  float *mr1;
  float *dr2;
  float *mr2;
  float dp1;
  float dp2;
  int it;
  int iq;
  int ih;
  float test;
  int ny=nt*nh;
  int nx=nt*nq;
  //////////////  

  if ((dr1=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr1=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");
  if ((dr2=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr2=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");

 
  for (it=0;it<nt;it++) for (ih=0;ih<nh;ih++) dr1[ih*nt+it]=frannor();
  for (it=0;it<nt;it++) for (iq=0;iq<nq;iq++) mr1[iq*nt+it]=frannor();

  oper(mr2,dr1,index,1,nt,nh,nq,nsparse);
  oper(mr1,dr2,index,0,nt,nh,nq,nsparse);

  dp1=dot(ny,dr1,dr2);
  dp2=dot(nx,mr1,mr2);

  if (dp2!=0) test=dp1/dp2;
  else test=0;

  fprintf(stderr,"Test adjoint = %f \n",test);
  return(test);
  
  free1float(mr2);
  free1float(dr2);
  free1float(mr1);
  free1float(dr1);

}


/****************************************************************************/

void Lumley_precond(float **M,float *t, float *h, int nt, int nh, float a, float b, float c)
{
  unsigned short it,ih;

  fprintf(stderr,"Inside Lumley nt=%d,nh=%d,a=%f,b=%f,c=%f\n",nt,nh,a,b,c);

  for (ih=0;ih<nh;ih++)
    for(it=0;it<nt;it++)
      M[ih][it]=1+c*(1+pow((fabs(h[ih]/1000.)),a))/(pow((1+t[it]),b));
    
  return;
}

void filterOutlier(int dataprec, float* h, int nh, float* t, int nt, float dt, float** d, float** Wd) {

    if (dataprec == 1 || (dataprec == 3 && 0)) {
        //radonhyp(m[0],Wd[0],index,0,nt,nh,nq,nsparse);
        fprintf(stderr, "Filter outliers with quantile............\n");
        float qup = quest(0.999, nh*nt, d[0]);
        float qmean = quest(0.50, nh*nt, d[0]);
        for (int ih = 0; ih < nh; ih++)
            for (int it = 0; it < nt; it++)
                if (fabs(d[ih][it]) > qup) Wd[ih][it] = qmean / qup;
    }
    if (0)
        if (dataprec == 2) { // Data preconditioner  
            fprintf(stderr, "Lumley Precond............\n");
            int c = 1;
            float a = 0.5;
            int b = 1;
            for (int ih = 0; ih < nh; ih++)
                for (int it = 0; it < nt; it++)
                    Wd[ih][it] = ((1 + 0.4 * sqrt(fabs(h[ih] / 1000.))) / (1 + 0 * t[it]));
            if (0) {
                save_gather(Wd, nh, nt, dt, "Wd");
                system("suxwigb < Wd  title=\"plot Wd\" &");
            }
        } else if (dataprec == 3) {
            float wdtemp;
            fprintf(stderr, "Wd = 1 ............\n");
            for (int ih = 1; ih < nh - 1; ih++) {
                wdtemp = fabs(h[ih - 1] - h[ih + 1]) / 2;
                for (int it = 0; it < nt; it++) Wd[ih][it] *= wdtemp;
            }
            for (int it = 0; it < nt; it++) {
                Wd[0][it] *= fabs(h[0] - h[1]);
                Wd[nh - 1][it] *= fabs(0 - h[nh - 2]) / 2;
            }
            if (1) {
                save_gather(Wd, nh, nt, dt, "Wd");
                system("suxwigb < Wd  title=\"plot Wd\" &");
            }
        } else if (dataprec == 4) {
            fprintf(stderr, "Wd = 1 ............\n");
            for (int ih = 0; ih < nh; ih++)
                for (int it = 0; it < nt; it++)
                    if (d[ih][it]) Wd[ih][it] = 1. / MAX(fabs(d[ih][it]), 1);
                    else Wd[ih][it] = 0;
            if (1) {
                save_gather(Wd, nh, nt, dt, "Wd");
                system("suxwigb < Wd  title=\"plot Wd\" &");
            }
        }
    return;

    
}













































