#ifndef RTMPAR_H
#define RTMPAR_H

#include <iostream>
#include <rsf.hh>
using namespace std;
typedef std::valarray<float> fmatrix;
class RtmPar{
public:
  RtmPar(iRSF& par);  
  ~RtmPar();
  void getModelSizes(iRSF& vel);
  void getDataSizes(iRSF& shots);
  void display();
  int nx;
  int nz;
  int ns;
  int ng;
  int nt;
  int nb;
  int nbup;
  float dx;
  float dz;
  float dt;
  float dxg;
  float amp;
  float fm;
  int sxbeg;
  int szbeg;
  int gxbeg;
  int gzbeg;
  int jsx;
  int jsz;
  int jgx;
  int jgz;
  int o1;
  int o2;
  int option;
  int applyLaplac;
  int mute;
  int niter;
  int order;
  float maxoffset;

};
#endif