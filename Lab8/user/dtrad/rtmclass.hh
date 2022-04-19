#include <iostream>
#include <rsf.hh>
using namespace std;
class RtmClass{
public:
  RtmClass(int nx, int nz, int nt, int nb, int nbup, int ng, int ns, float dx, float dz, 
	 float dt, float* vel, float fm, float amp);
  ~RtmClass();

  void sg_init(int szbeg, int sxbeg, int jsz, int jsx, string type);
  void transpose(float* trans);  
  float* getWavelet();
  float* getModel();
  void setRecording(sf_file& wavefile, int ft, int jt);
  void loop(bool adj, bool add, float* mod, float* dat);
  void initGeometry(int xbeg, int zbeg, int jx, int jz, int n, string type);
  void adjtest();
  void adjtestLaplac();
  int  getNx();
  int  getNz();
  int  getNt();
  int  getNg();
  int  getNs();

private:
  void init(float* vels);  
  void sg_init(int* sxz, int szbeg, int sxbeg, int jsz, int jsx, int ns);
  void expand2d(float** b, float** a);
  void window2d(float **a, float **b);
  void boundaryRW(float** p, int it, bool read);
  void rw_snapshot(float** p, int it, bool read);
  void rcvrTimeSlice(float** g, float* data, int is, int it, bool adj);
  void stepForward(float** p0, float** p1, bool adj);
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
  float compareWavefields(int it, float **u);
  void applyLaplac(float* r, float* p, int adj);

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

  // FD coeff
  float myC0;
  float myC11;
  float myC12;
  float myC21;
  float myC22;

  // data, vel and model arrays
  float** myV0;
  float** myVV;

  float** mySP0;
  float** mySP1;
  float** myGP0;
  float** myGP1;
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

};
