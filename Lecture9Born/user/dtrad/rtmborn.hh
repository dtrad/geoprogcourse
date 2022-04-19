#ifndef RTMBORN_H
#define RTMBORN_H

#include <iostream>
#include <rsf.hh>
#include "rtmpar.hh"
using namespace std;
class RtmBorn{
public:
  RtmBorn(RtmPar* par, float* vel);
  ~RtmBorn();

  void sg_init(int szbeg, int sxbeg, int jsz, int jsx, string type);
  void transpose(float* trans);  
  float* getWavelet();
  float* getModel();
  void setRecording(sf_file& wavefile, int ft, int jt);
  void loop(bool adj, bool add, float* mod, float* dat);
  void initGeometry(int xbeg, int zbeg, int jx, int jz, int n, string type);
  void adjtest();
  void adjtestLaplac();
  void mute(int is, float* data);
  void muteoffset(int is, float* data);
  void setIter(int iter);
  int  getNx();
  int  getNz();
  int  getNt();
  int  getNg();
  int  getNs();
  float getDt();
  float getDx();
  float getDz();
  float getDg();
  double dot(float* a, float *b, size_t n);

private:
  void init(float* vels);  
  void sg_init(int* sxz, int szbeg, int sxbeg, int jsz, int jsx, int ns);
  void expand2d(float** b, float** a);
  void window2d(float **a, float **b);
  void boundaryRW(float** p, int it, bool read);
  void rw_snapshot(float** p, int it, bool read);
  void rcvrTimeSlice(float** g, float* data, int is, int it, bool adj);
  void stepForward(float** p0, float** p1, bool adj);
  void stepForwardEO(float** p0, float** p1, bool adj);
  void stepForwardFO(float** p0, float** p1, bool adj);
  void stepForwardSO(float** p0, float** p1, bool adj);
  void applySponge();
  void applySponge(float**p0);

  void addSource(int is, int it, bool add);
  void addSource(int *sxz, float **p, int ns, float *source, bool add);

  void matrixTranspose(float *matrix, float *trans, int n1, int n2);
  void saveWave(int it, float** p);
  bool isStable(float* velp);
  void setBndr(float*& bndr, int nb);
  void imaging(float* m, float** s, float **g, bool adj);
  int getTime(struct timeval& t1, struct timeval& t2);

  void applyLaplac(float* r, float* p, int adj);
  void applyLaplac(float** r, float** p, int adj);
  void applyLaplac(float* r, float* p, int adj, int n1, int n2);

  //Born
  void sourceTimeDerivative();
  void timeDerivative(float** p1, float** p2); // p_tt=p2+2*p1-p0;
  void copyWavefield(float** p, float** g); // g=p
  void imagingBorn(float* m, float** ds2, float**g, bool adj);

  // dimensions
  int myNxpad;
  int myNzpad;
  int myNx;
  int myNz;
  int myNb;
  int myNbUp;
  int myNt;
  int myNg;
  int myNs;
  
  // sampling
  float myDx;
  float myDz;
  float myDt;
  float myDg;

  // FD coeff
  float myC0;
  float myC11;
  float myC12;
  float myC21;
  float myC22;
  
  ////additions for EightOrder coeff
  float myEC4xx;
  float myEC3xx;
  float myEC2xx;
  float myEC1xx;
  float myEC0xx;

  float myEC4zz;
  float myEC3zz;
  float myEC2zz;
  float myEC1zz;
  float myEC0zz;
  

  // data, vel and model arrays
  float** myV0; // use as a temporary array nx,nz
  float** myVV;  // vel*vel*dt*dt
  float** myVel; // untouched vel

  float** mySP0;
  float** mySP1;
  float** myGP0;
  float** myGP1;
  float** mySPtt;   // added for Born
  float*** mySp0array;

  float*  myDObs;
  float*  myTrans;
  float*  myBndr;
  float*  myBndrUp;
  float*  myWlt;
  float*  myRWBndr;

  // source
  float myFm;
  float myAmp;

  // recording coordinates
  int*  mySxz;
  int*  myGxz;

  float myMaxOffset;
  
  // auxiliar
  float myPI;

  // temporary image for laplac
  float* myModTmp;

  sf_file myWavefile;
  int myFTrecord;
  int myJTrecord;
  bool myVerbose;

  int myNthreads;
  int myChunksize;
  int myFromBoundary;
  int myOption; // option for imaging condition (correlaiton or Born)
  int myApplyLaplac; // apply laplacian to the image 
  int myMute; // apply mute to first arrival.
  int myIter; // iteration number (0 is not LSMIG)
  int myNIter; // total number of iterations (if LSMIG)
  float myDWVel;
  
  // addition for eight order staggered grid
  bool myEightOrder;
  bool myFourthOrder;
  bool mySecondOrder;
  
  
};
#endif
