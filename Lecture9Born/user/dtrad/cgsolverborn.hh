#ifndef CGSOLVERBORN_H
#define CGSOLVERBORN_H

#include <iostream>
#include <rsf.hh>
#include "MyAlloc.hpp"

using namespace std;
typedef std::valarray<float> fmatrix;

class RtmBorn;
class CGSolverBorn{
public:
  CGSolverBorn(RtmBorn* rtm);
  ~CGSolverBorn();

  float wpcgnr(fmatrix& data, fmatrix& model, fmatrix& residuals, fmatrix& predictions, int itercg);
  float cgsq(fmatrix& data, fmatrix& model, fmatrix& residuals, fmatrix& predictions, int itercg);
  float dot(fmatrix& data1, fmatrix& data2);
  float lsScale(fmatrix& data1, fmatrix& data2);
  
  // obsolete
  double dot(fmatrix& data1, fmatrix& data2, int nx, int nz);  
private:
  RtmBorn* myRtm;
  void setModelFile(int itercg);
  void setResidualFile(int itercg);
  void setResidualModelFile(int itercg);
  sf_file myModelFile;
  sf_file myResidualFile;
  sf_file myResidualModelFile;

};
#endif