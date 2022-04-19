#include <iostream>
#include <rsf.hh>


using namespace std;
class Fdmodb{
public:
  Fdmodb(int nx, int nz, int nt, int nb, int nbup, int ng, int ns, float dx, float dz, float dt, 
	 float* vel, float fm, float amp, int order);
  ~Fdmodb();

  void sg_init(int szbeg, int sxbeg, int jsz, int jsx, string type);
  void modelling(int is);
  void transpose(float* trans);
  void transposeAdd(float* trans); 
  float* getWavelet();
  void setRecording(sf_file& wavefile, int ft, int jt);
  

private:
  void init(float* vels);  
  void sg_init(int* sxz, int szbeg, int sxbeg, int jsz, int jsx, int ns);
  void expand2d(float** b, float** a);
  void window2d(float **a, float **b);
  
  void stepForward();

  void stepForwardFO(float** p0, float** p1);
  void stepForwardSO(float** p0, float** p1);
  void stepForwardEO(float** p0, float**p1);

  void applySponge();
  void applySponge(float**p0);


  void recordSeis(int it);
  void recordSeis(float *seis_it, int *gxz, float **p, int ng);

  void addSource(int is, int it);
  void addSource(int *sxz, float **p, int ns, float *source, bool add);

  void circlePointers();

 
  void matrixTranspose(float *matrix, float *trans, int n1, int n2);
  void saveWave(int it);
  bool isStable(float* velp);

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
  float** myP0;
  float** myP1;
  float*  myDObs;
  float*  myTrans;
  float*  myBndr;
  float*  myWlt;

  // source
  float myFm;
  float myAmp;

  // recording coordinates
  int*  mySxz;
  int*  myGxz;

  // auxiliar
  float myPI;

  // addition for eight order staggered grid
  bool myEightOrder;
  bool myFourthOrder;
  bool mySecondOrder;

  sf_file myWavefile;
  int myFTrecord;
  int myJTrecord;

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
};
