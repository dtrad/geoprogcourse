#include <iostream>
#include <cstdlib>
#include <cstring> // need memcpy
#include <cmath>
#include "rtmclass.hh"
#include "MyAlloc.hpp"
#include <time.h>
#include <stdexcept>
#include <sys/time.h>


#ifdef _OPENMP
#include <omp.h>
#endif

// #ifdef _OPENMP
// #undef _OPENMP
// #endif

#ifndef MARK
#define MARK fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif

void RtmClass::stepForward(float** p0, float** p1, bool adj)
/*< original forward modeling step fourth order>*/
{
  // implement Forward and adjoint here.

}

// class for forward modeling with second, fourth and eight order staggered grid.
RtmClass::RtmClass(int nx, int nz, int nt, int nb, int nbup, int ng, int ns, float dx, float dz, float dt, float* vel, float fm, float amp){
  cerr << " RTM class constructor " << endl;
  myDz = dz;
  myDx = dx;
  myDt = dt;

  myNx = nx;
  myNz = nz;
  myNt = nt;

  myNg = ng;
  myNs = ns;
  
  myFm = fm;
  myAmp = amp;

  myPI = acos(-1.);
  myNb = nb;
  myNbUp = nbup; // nb at surface. =0 to produce multiples.

  mySP0 = 0;
  mySP1 = 0;
  myGP0 = 0;
  myGP1 = 0;


  myDObs = 0;
  myTrans = 0;
  myVV    = 0;
  myV0    = 0;

  myBndr   = 0;
  myBndrUp = 0;
  myWlt    = 0;
  mySxz    = 0;
  myGxz    = 0;
  myModTmp = 0;
  myVerbose = true;
  myNthreads = 1;
  myChunksize = 1;
  myFromBoundary = 1;


#ifdef _OPENMP
  myNthreads = omp_get_max_threads();
  cerr << "Using " << myNthreads << " threads " << endl;
  myChunksize=(myNx+myNb+myNb)/myNthreads/8;
  cerr << "chunksize" << myChunksize << endl;
#endif

  
  init(vel);
  if (!isStable(vel)){
    cerr << "WARNING********UNSTABLE*********" << endl;
    //throw runtime_error("UNSTABLE");
  }
  else cerr << "STABLE!" << endl;
  cerr << " ********* RTM constructor done ******** " << endl;

  
}
RtmClass::~RtmClass(){
  cerr << " calling RTM destructor " << endl;

  if (myV0)     del2d(myV0);
  if (myVV)     del2d(myVV);
  if (mySP0)    del2d(mySP0);
  if (mySP1)    del2d(mySP1);
  if (myGP0)    del2d(myGP0);
  if (myGP1)    del2d(myGP1);
  if (myDObs)   del1d(myDObs);

  if (myTrans)  del1d(myTrans);
  if (myBndr)   del1d(myBndr);
  if (myBndrUp) del1d(myBndrUp);
  if (myWlt)    del1d(myWlt);
  if (mySxz)    del1d(mySxz);
  if (myGxz)    del1d(myGxz);

  if (myModTmp) del1d(myModTmp);

  if (myRWBndr) del1d(myRWBndr);
  if (!myFromBoundary) if (mySp0array) del3d(mySp0array);

  cerr << " RTM destructor done " << endl;
}

void RtmClass::setBndr(float*& bndr, int nb){
  // initialize sponge ABC coefficients 
  new1d(bndr,nb);
  float t;
  for(int ib=0;ib<nb;ib++){
    t=0.015*(nb-1-ib);
    bndr[ib]=expf(-t*t);
  }

}

void RtmClass::init(float* vel){

  myNzpad=myNz+myNb+myNbUp;
  myNxpad=myNx+myNb+myNb;


  myC11 = 4./3./(myDz*myDz);
  myC21 = 4./3./(myDx*myDx);
  myC12 = -1./12./(myDz*myDz);
  myC22 = -1./12./(myDx*myDx);
  myC0  = -2.*(myC11+myC12+myC21+myC22);
  new2d(myV0,myNx,myNz);
  new2d(myVV,myNxpad,myNzpad);
  new2d(mySP0,myNxpad,myNzpad);
  new2d(mySP1,myNxpad,myNzpad);
  new2d(myGP0,myNxpad,myNzpad);
  new2d(myGP1,myNxpad,myNzpad);
  if (!myFromBoundary) new3d(mySp0array,myNt,myNxpad,myNzpad);

  new1d(myDObs,myNg*myNt);
  new1d(myTrans,myNg*myNt);
  new1d(myModTmp,myNx*myNz);
  new1d(myWlt,myNt);
  new1d(mySxz,myNs);
  new1d(myGxz,myNg);
  new1d(myRWBndr,myNt*4*(myNx+myNz));
  
  memcpy(myV0[0],vel,myNx*myNz*sizeof(float));
  memset(myDObs,0,myNg*myNt*sizeof(float));
  memset(myTrans,0,myNg*myNt*sizeof(float));

  setBndr(myBndr,myNb);
  setBndr(myBndrUp,myNbUp);
  
  expand2d(myVV,myV0);
  float tmp;
  for (int ix=0;ix < myNxpad;ix++){
    for (int iz=0;iz < myNzpad;iz++){
      tmp = myVV[ix][iz]*myDt;
      myVV[ix][iz] = tmp*tmp;
    }
  }

  for(int it=0; it < myNt; it++){
    tmp = myPI*myFm*(it*myDt - 1.0/myFm);
    tmp = tmp*tmp;
    myWlt[it] = myAmp*(1.-2.*tmp)*expf(-tmp);
    //if (myWlt[it]) cerr << "source[" << it << "]="<< myWlt[it] << endl;
  }
  
  

}

float* RtmClass::getWavelet(){
  return myWlt;
}

void RtmClass::expand2d(float** b, float** a){
#ifdef _OPENMP
#pragma omp parallel for 
#endif
  for (int ix=0;ix<myNx;ix++) {
    for (int iz=0;iz<myNz;iz++) {
      b[myNb+ix][myNbUp+iz] = a[ix][iz];
    }
  }
  
  
  for (int ix=0; ix<myNxpad; ix++){ 
    for (int iz=0; iz<myNbUp;    iz++){  // top
      b[ix][          iz] = b[ix][myNbUp        ];
    }

    for (int iz=0; iz<myNb;    iz++){ 
      b[ix][myNzpad-iz-1] = b[ix][myNzpad-myNb-1];/* bottom*/
    }
  }
  for (int ix=0; ix<myNb;    ix++) {
    for (int iz=0; iz<myNzpad; iz++) {
      b[ix 	     ][iz] = b[myNb  		][iz];/* left */
      b[myNxpad-ix-1 ][iz] = b[myNxpad-myNb-1	][iz];/* right */
    }
  }
  
}


void RtmClass::window2d(float **a, float **b)
/*< window 'b' to 'a': source(b)-->destination(a) >*/
{
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for (int ix=0;ix<myNx;ix++) {
      for (int iz=0;iz<myNz;iz++) {
	a[ix][iz]=b[myNb+ix][myNbUp+iz] ;
      }
    }
}

void RtmClass::rw_snapshot(float** p, int it, bool read){
  // read/write snapshot completely 
  // if read=true read, else write
  
  if (!read){
    for (int ix=0;ix< myNxpad;ix++){
      for (int iz=0;iz< myNzpad;iz++){
	mySp0array[it][ix][iz]=p[ix][iz];
      }
    }
  }
  else{
    for (int ix=0;ix< myNxpad;ix++){
      for (int iz=0;iz< myNzpad;iz++){
	p[ix][iz]=mySp0array[it][ix][iz];
      }
    }
  }

}


void RtmClass::boundaryRW(float** p, int it, bool read){

  float* spo = &myRWBndr[it*4*(myNx+myNz)];
  // read/write using effective boundary saving strategy
  // if read = true, read the boundary out, else save/write boundary
  if (read){
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int ix=0;ix<myNx;ix++){
      for (int iz=0; iz<2; iz++){
	p[ix+myNb][iz-   2+myNb]=spo[iz+4*ix];
	p[ix+myNb][iz+myNz+myNb]=spo[iz+2+4*ix];
      }
    }
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int iz=0;iz<myNz;iz++){
      for (int ix=0;ix<2;ix++){
	p[ix-   2+myNb][iz+myNb]=spo[4*myNx+iz+myNz*ix];
	p[ix+myNx+myNb][iz+myNb]=spo[4*myNx+iz+myNz*(ix+2)];
      }
    }
  }
  else{
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int ix=0;ix<myNx;ix++){
      for (int iz=0; iz<2; iz++){
	spo[iz+  4*ix]=p[ix+myNb][iz-   2+myNb];
	spo[iz+2+4*ix]=p[ix+myNb][iz+myNz+myNb];
      }
    }
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int iz=0;iz<myNz;iz++){
      for (int ix=0;ix<2;ix++){
	spo[4*myNx+iz+myNz*ix]    =p[ix-   2+myNb][iz+myNb];
	spo[4*myNx+iz+myNz*(ix+2)]=p[ix+myNx+myNb][iz+myNb];
      }
    }
  }

}

void RtmClass::initGeometry(int xbeg, int zbeg, int jx, int jz, int n, string type){
    // set geometry ////////////////////////////////////
    if (!(xbeg>=0 && zbeg>=0 && xbeg+(n-1)*jx<myNx && zbeg+(n-1)*jz<myNz)){
      cerr << type << " exceeds the computing zone" << endl;
      exit(1);
    }
    sg_init(zbeg, xbeg, jz, jx, type);
}



void RtmClass::sg_init(int szbeg, int sxbeg, int jsz, int jsx, string type){
  if      (type == "shots") sg_init(mySxz,szbeg,sxbeg,jsz,jsx,myNs);
  else if (type == "rcvrs") sg_init(myGxz,szbeg,sxbeg,jsz,jsx,myNg);
}

void RtmClass::sg_init(int* sxz, int szbeg, int sxbeg, int jsz, int jsx, int ns){
  int sx, sz;
  for (int is=0;is < ns;is++){
    sz=szbeg+is*jsz;
    sx=sxbeg+is*jsx;
    sxz[is]=myNz*sx+sz;
  }
}


void RtmClass::applySponge(float**p0)
/*< apply absorbing boundary condition >*/
{
  int ix,iz,ib,ibx,ibz;
  float w;

#ifdef _OPENMP
#pragma omp parallel for  private(ib,iz,ix,ibz,ibx,w)		
#endif
  for(ib=0; ib<myNb; ib++) {
    w = myBndr[ib];

    ibz = myNzpad-ib-1;
    for(ix=0; ix<myNxpad; ix++) {
      p0[ix][ibz] *= w; /* bottom sponge */
    }

    ibx = myNxpad-ib-1;
    for(iz=0; iz<myNzpad; iz++) {
      p0[ib ][iz] *= w; /*   left sponge */
      p0[ibx][iz] *= w; /*  right sponge */
    }
  }

  for(ib=0; ib<myNbUp; ib++) {
    for(ix=0; ix<myNxpad; ix++) {
      p0[ix][ib ] *= myBndrUp[ib]; // upper sponge
    }
  }

}

void RtmClass::addSource(int is, int it, bool add){
  addSource(&mySxz[is],mySP1,1,&myWlt[it],add);
}


void RtmClass::addSource(int *sxz, float **p, int ns, float *source, bool add)
/*< add source term >*/
{
  int is, sx, sz;
  if(add){
    for(is=0;is<ns; is++){
      sx=sxz[is]/myNz+myNb;
      sz=sxz[is]%myNz+myNbUp;
      p[sx][sz]+=source[is];
    }
  }else{
    for(is=0;is<ns; is++){
      sx=sxz[is]/myNz+myNb;
      sz=sxz[is]%myNz+myNbUp;
      p[sx][sz]-=source[is];
    }
  }
}



void RtmClass::transpose(float* trans){
  matrixTranspose(myDObs,trans,myNg,myNt);
}


void RtmClass::matrixTranspose(float *matrix, float *trans, int n1, int n2)
/*< matrix transpose: matrix tansposed to be trans >*/
{
  int i1, i2;
  
  for(i2=0; i2<n2; i2++)
    for(i1=0; i1<n1; i1++)
      trans[i2+n2*i1]=matrix[i1+n1*i2];
}


void RtmClass::setRecording(sf_file& wave, int ft, int jt){
  myWavefile = wave;
  myFTrecord = ft;
  myJTrecord = jt;
}

void RtmClass::saveWave(int it, float** p){
  if (!((it-myFTrecord)%myJTrecord)){
    window2d(myV0,p);
    sf_floatwrite(myV0[0],myNx*myNz,myWavefile);
  }
  
}

bool RtmClass::isStable(float* velp){
  float vmax= 0;
  size_t n = myNx*myNz;
  for (size_t i=0;i<n;i++) vmax= max(vmax,velp[i]);

  bool stable = false;
  float factor = (vmax*myDt*sqrt(1./(myDx*myDx)+1./(myDz*myDz)));
  if (factor < 1) stable=true;
  
  cerr << " vmax = " << vmax << ", dx = " << myDx << ", dz=" 
       << myDz << ", dt " << myDt << ", factor= " << factor << endl;

  return stable;
  
}
void RtmClass::adjtest(){

  // implement ADJOINT TEST HERE
  fprintf(stderr,"YOUR JOB");
}

void RtmClass::adjtestLaplac(){
  // YOUR JOB
  fprintf(stderr,"YOUR JOB");
  
}


void RtmClass::loop(bool adj, bool add, float* model, float* data){
  //cerr << " starting loop " << adj << endl;
  if (adj) memset(myModTmp,0,myNx*myNz*sizeof(float));
  else     memset(data,0,myNt*myNg*myNs*sizeof(float));

  time_t start;
  time_t end;
  double seconds;
  float** ptr;
  struct timeval t1, t2; // millisecond timing

  for (int is=0;is<myNs;is++){
    cerr << "shot = " << is+1 << " of " << myNs << endl;
    time(&start); // time in sec
    int totalTime= 0; // time in msec
    int totalTime2 = 0;    
    // initialize is-th source wavefield Ps
    memset(mySP0[0], 0, myNzpad*myNxpad*sizeof(float));
    memset(mySP1[0], 0, myNzpad*myNxpad*sizeof(float));
    memset(myGP0[0], 0, myNzpad*myNxpad*sizeof(float));
    memset(myGP1[0], 0, myNzpad*myNxpad*sizeof(float));
    if (adj){ // migration mm=Lt dd
      gettimeofday(&t1,NULL);
      for (int it=0;it<myNt;it++){ // in adj shots go forward, g reverse
	addSource(is,it,true);
	stepForward(mySP0,mySP1,false);
	//stepForward(mySP0,mySP1,false,c0,c11,c12,c21,c22);
	ptr=mySP0; mySP0=mySP1; mySP1=ptr;
	applySponge(mySP0);
	applySponge(mySP1);
	if (myFromBoundary) boundaryRW(mySP0,it, false);
	else                rw_snapshot(mySP0,it,false);
      }
      for (int it=myNt-1;it > -1;it--){
	// reverse time order, Img[]+=Ps[]*Pg[];
	// backpropagate receiver wavefield
	rcvrTimeSlice(myGP1,data,is,it,adj);
	if (myFromBoundary){
	  boundaryRW(mySP0, it, true);
	  //compareWavefields(it,sp0);
	  imaging(myModTmp,mySP0,myGP1,adj);
	  // reconstruct source wavefield Ps
	  ptr=mySP0; mySP0=mySP1; mySP1=ptr;
	  stepForward(mySP0,mySP1,false);
	  addSource(is,it,false);
	}
	else{
	  rw_snapshot(mySP0,it,true); // read source wavefield.
	  imaging(myModTmp,mySP0,myGP1,adj);
	}
	stepForward(myGP0,myGP1,false);
	ptr=myGP0; myGP0=myGP1; myGP1=ptr;
	applySponge(myGP0);
	applySponge(myGP1);
	
      }
      applyLaplac(myModTmp,model,adj);
      gettimeofday(&t2,NULL);	totalTime2 += getTime(t1,t2);
    }
    else{
      // in forward both wavefield go in same direction (no need for RWBndr)
      gettimeofday(&t1,NULL);
      applyLaplac(myModTmp,model,adj);
      for (int it=0; it<myNt;it++){
	// if ((myVerbose)&&(!(it%100))) cerr << it << endl;
	// forward time order Pg += Ps * Img;
	addSource(is,it,true);
	stepForward(mySP0,mySP1,false); // not sure why is false here
	ptr=mySP0;mySP0=mySP1;mySP1=ptr;
	applySponge(mySP0);
	applySponge(mySP1);
	imaging(myModTmp,mySP0,myGP1,adj);
	rcvrTimeSlice(myGP1,data,is,it,adj);
	stepForward(myGP0, myGP1, true); // true is forward
	ptr=myGP0; myGP0=myGP1; myGP1=ptr;
	applySponge(myGP1);
	applySponge(myGP0);
      }
      gettimeofday(&t2,NULL);	totalTime += getTime(t1,t2);      
    }

    cerr << "timeForward =" << totalTime << endl;
    cerr << "timeBackward =" << totalTime2 << endl;
    time(&end);
    seconds = difftime(end,start);
    sf_warning("shot %d finished: %f\n", is+1, seconds);

  }
}

void RtmClass::rcvrTimeSlice(float** g, float* data, int is, int it, bool adj){
  int gx, gz;
  if (adj){
    for (int ig=0;ig<myNg;ig++){
      gx=myGxz[ig]/myNz;
      gz=myGxz[ig]%myNz;
      g[gx+myNb][gz+myNbUp] += data[it+ig*myNt+is*myNt*myNg];
    }
  }
  else{
    for (int ig=0;ig<myNg;ig++){
      gx=myGxz[ig]/myNz;
      gz=myGxz[ig]%myNz;
      data[it+ig*myNt+is*myNt*myNg] += g[gx+myNb][gz+myNbUp];
    }
  }
    
}


void RtmClass::imaging(float* m, float** s, float **g, bool adj){
  if (adj){
    for (int i2=0;i2<myNx;i2++)
      for (int i1=0;i1<myNz;i1++)
	m[myNz*i2+i1] += 
	  s[i2+myNb][i1+myNbUp]*g[i2+myNb][i1+myNbUp];
  }
  else{
    for (int i2=0;i2<myNx;i2++)
      for (int i1=0;i1<myNz;i1++)
	g[i2+myNb][i1+myNbUp] +=
	  s[i2+myNb][i1+myNbUp]*m[i1+myNz*i2];

  }

}

int RtmClass::getTime(struct timeval& t1, struct timeval& t2){
  // return time in milliseconds;
  return ((t2.tv_sec - t1.tv_sec)* 1000 + (t2.tv_usec - t1.tv_usec)/1000);
}

float RtmClass::compareWavefields(int it, float **u){
// compare current wavefield with true wavefield saved in sp0array 
  float sum =0;
  for (int ix=myNb;ix<myNxpad-myNb;ix++){
    for (int iz=myNb;iz<myNzpad-myNb;iz++){
      if (myGP1[ix][iz]) sum += fabs((u[ix][iz]-mySp0array[it][ix][iz]));
    }
  }
  return sum;
}

int RtmClass::getNx(){ return myNx;}
int RtmClass::getNz(){ return myNz;}
int RtmClass::getNt(){ return myNt;}
int RtmClass::getNg(){ return myNg;}
int RtmClass::getNs(){ return myNs;}

void RtmClass::applyLaplac(float* r, float* p, int adj){
  int n2=myNx;
  int n1=myNz;
  int n12=myNx*myNz;
  bool disableLaplacian=true; // disable laplacian. 
  
  if (disableLaplacian){ 
    if (adj) memcpy(p,r,n12*sizeof(float));
    else     memcpy(r,p,n12*sizeof(float)); 
    return;
  }
  
  if  (adj) memset(p,0,n12*sizeof(float)); 
  else      memset(r,0,n12*sizeof(float));
  for  (int i2=0; i2 < n2; i2++) {
    for (int i1=0; i1 < n1; i1++) {
      int j = i1+i2*n1;
      if (i1 > 0) {
	if (adj) {
	  p[j-1] -= r[j];
	  p[j]   += r[j];
	} else {
	  r[j] += p[j] - p[j-1];
	}
      }
      if (i1 < n1-1) {
	if (adj) {
	  p[j+1] -= r[j];
	  p[j]   += r[j];
	} else {
	  r[j] += p[j] - p[j+1];
	}
      }
      
      if (i2 > 0) {
	if (adj) {
	  p[j-n1] -= r[j];
	  p[j]    += r[j];
	} else {
	  r[j] += p[j] - p[j-n1];
	}
      }
      if (i2 < n2-1) {
	if (adj) {
	  p[j+n1] -= r[j];
	  p[j]    += r[j];
	} else {
	  r[j] += p[j] - p[j+n1];
	}
      }
    }
  }

  return;
}
