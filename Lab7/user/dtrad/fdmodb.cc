#include <iostream>
#include <cstdlib>
#include <cstring> // need memcpy
#include <cmath>
#include "fdmodb.hh"
#include "MyAlloc.hpp"
#include <time.h>
#include <stdexcept>


#ifndef MARK
#define MARK fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif

// class for forward modeling with second, fourth and eight order staggered grid.
Fdmodb::Fdmodb(int nx, int nz, int nt, int nb, int nbup, int ng, int ns, float dx, float dz, float dt, float* vel, float fm, float amp, int method){
  cerr << " constructor ... " << endl;
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



  myEightOrder = false; // eight order has problems, high amp at center
  myFourthOrder = false;//true; // work best
  mySecondOrder = false;// faster but poorer

  if (method == 1) mySecondOrder = true;
  else if (method == 2) myFourthOrder = true;
  else if (method == 3) myEightOrder = true;
  else cerr << "method not implemented = " << method << endl;

  myP0 = 0;
  myP1 = 0;
  myDObs = 0;
  myTrans = 0;
  myVV    = 0;
  myV0    = 0;

  init(vel);
  if (!isStable(vel)){
    cerr << "WARNING********UNSTABLE*********" << endl;
    //throw runtime_error("UNSTABLE");
  }
  else cerr << "STABLE!" << endl;
  cerr << " done " << endl;

  ////for eight order
  myEC4xx = -1./560./(myDx*myDx);
  myEC3xx = 8./315./(myDx*myDx);
  myEC2xx = -1./5./(myDx*myDx);
  myEC1xx = 8./5./(myDx*myDx);
  myEC0xx = -205./72./(myDx*myDx);

  myEC4zz = -1./560./(myDz*myDz);
  myEC3zz = 8./315./(myDz*myDz);
  myEC2zz = -1./5./(myDz*myDz);
  myEC1zz = 8./5./(myDz*myDz);
  myEC0zz = -205./72./(myDz*myDz);

  
}
Fdmodb::~Fdmodb(){
  cerr << " destructor " << endl;
  if (myV0) del2d(myV0);
  if (myVV) del2d(myVV);
  if (myP0) del2d(myP0);
  if (myP1) del2d(myP1);
  if (myDObs) del1d(myDObs);
  if (myTrans) del1d(myTrans);
  if (myBndr)  del1d(myBndr);
  if (myWlt)   del1d(myWlt);
  if (mySxz)   del1d(mySxz);
  if (myGxz)   del1d(myGxz);

  cerr << " done " << endl;
}
void Fdmodb::init(float* vel){

  myNzpad=myNz+myNb+myNbUp;
  myNxpad=myNx+myNb+myNb;

  myC11 = 4./3./(myDz*myDz);
  myC21 = 4./3./(myDx*myDx);
  myC12 = -1./12./(myDz*myDz);
  myC22 = -1./12./(myDx*myDx);
  myC0  = -2.*(myC11+myC12+myC21+myC22);

  new2d(myV0,myNx,myNz);
  new2d(myVV,myNxpad,myNzpad);
  new2d(myP0,myNxpad,myNzpad);
  new2d(myP1,myNxpad,myNzpad);
  new1d(myDObs,myNg*myNt);
  new1d(myTrans,myNg*myNt);
  new1d(myBndr,myNb);
  new1d(myWlt,myNt);
  new1d(mySxz,myNs);
  new1d(myGxz,myNg);


  cerr << myNx << ", " << myNz << ", " << myNg << ", " << myNt << endl;
  memcpy(myV0[0],vel,myNx*myNz*sizeof(float));
  memset(myDObs,0,myNg*myNz*sizeof(float));
  memset(myTrans,0,myNg*myNz*sizeof(float));
  float tmp; 
  for (int ib=0;ib < myNb;ib++){
    tmp = 0.015*(myNb-ib);
    myBndr[ib]=expf(-tmp*tmp);
  }

  
  expand2d(myVV,myV0);

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

 
float* Fdmodb::getWavelet(){
  return myWlt;
}

void Fdmodb::expand2d(float** b, float** a){
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


void Fdmodb::window2d(float **a, float **b)
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

void Fdmodb::sg_init(int szbeg, int sxbeg, int jsz, int jsx, string type){
  if      (type == "shots") sg_init(mySxz,szbeg,sxbeg,jsz,jsx,myNs);
  else if (type == "rcvrs") sg_init(myGxz,szbeg,sxbeg,jsz,jsx,myNg);
}

void Fdmodb::sg_init(int* sxz, int szbeg, int sxbeg, int jsz, int jsx, int ns){
  int sx, sz;
  for (int is=0;is < ns;is++){
    sz=szbeg+is*jsz;
    sx=sxbeg+is*jsx;
    sxz[is]=myNz*sx+sz;
  }
}


void Fdmodb::modelling(int is){
  memset(myP1[0],0,myNzpad*myNxpad*sizeof(float));
  //initializing P1 of given size with 0 
  memset(myP0[0],0,myNzpad*myNxpad*sizeof(float));

  for (int it = 0; it< myNt;it++){
    
    if (!(it%1000)) cerr << "it = " << it << endl;   
    addSource(is,it);
    stepForward();
    circlePointers();
    if (is == 0) saveWave(it); // save first shot only
    applySponge();
    recordSeis(it);

  }

}


void Fdmodb::stepForward(){
  if      (mySecondOrder) stepForwardSO(myP0,myP1);
  else if (myFourthOrder) stepForwardFO(myP0,myP1);
  else if (myEightOrder)  stepForwardEO(myP1,myP0);
}

void Fdmodb::stepForwardSO(float** p0, float** p1)
/*< forward modeling step second order>*/
{
  // YOUR JOB


}

void Fdmodb::stepForwardFO(float** p0, float** p1)
/*< forward modeling step fourth order>*/
{

  // YOUR JOB
}

void Fdmodb::stepForwardEO(float** p0,float** p1)
//intuitive version of 8th order FDM
{
  // YOUR JOB
}

void Fdmodb::applySponge(){
  applySponge(myP0);
  applySponge(myP1);
}


void Fdmodb::applySponge(float**p0)
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
      p0[ix][ib ] *= myBndr[ib]; // upper sponge
    }
  }

}

void Fdmodb::addSource(int is, int it){
  addSource(&mySxz[is],myP1,1,&myWlt[it],true);
}


void Fdmodb::addSource(int *sxz, float **p, int ns, float *source, bool add)
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


void Fdmodb::recordSeis(int it){
    recordSeis(&myDObs[it*myNg], myGxz, myP0, myNg);
}


void Fdmodb::recordSeis(float *seis_it, int *gxz, float **p, int ng)
/*< record seismogram at time it into a vector length of ng >*/
{
  int ig, gx, gz;
  for(ig=0;ig<myNg; ig++){
    gx=gxz[ig]/myNz+myNb;
    gz=gxz[ig]%myNz+myNbUp;
    seis_it[ig]=p[gx][gz];
  }
}

void Fdmodb::transpose(float* trans){
  matrixTranspose(myDObs,trans,myNg,myNt);
}


void Fdmodb::matrixTranspose(float *matrix, float *trans, int n1, int n2)
/*< matrix transpose: matrix tansposed to be trans >*/
{
  int i1, i2;
  
  for(i2=0; i2<n2; i2++)
    for(i1=0; i1<n1; i1++)
      trans[i2+n2*i1]=matrix[i1+n1*i2];
}




void Fdmodb::circlePointers(){
  float** ptr = myP0;
  myP0 = myP1;
  myP1 = ptr;
  

}

void Fdmodb::setRecording(sf_file& wave, int ft, int jt){
  myWavefile = wave;
  myFTrecord = ft;
  myJTrecord = jt;
}

void Fdmodb::saveWave(int it){
  if (!((it-myFTrecord)%myJTrecord)){
    window2d(myV0,myP1);
    sf_floatwrite(myV0[0],myNx*myNz,myWavefile);
  }
  
}

bool Fdmodb::isStable(float* velp){
  bool stable = false;
  // YOUR JOB
  // return true if stable.


  return stable;
  
}

