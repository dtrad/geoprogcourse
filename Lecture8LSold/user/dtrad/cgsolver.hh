#include <iostream>
#include <rsf.hh>
#include "MyAlloc.hpp"

using namespace std;
typedef std::valarray<float> fmatrix;

class RtmClass;
class CGSolver{
public:
  CGSolver(RtmClass* rtm);
  ~CGSolver();

  float wpcgnr(fmatrix& data, fmatrix& model, fmatrix& residuals, fmatrix& predictions, int itercg);
  float dot(fmatrix& data1, fmatrix& data2);
  float lsScale(fmatrix& data1, fmatrix& data2);

  // obsolete
  double dot(fmatrix& data1, fmatrix& data2, int nx, int nz);  
private:
  RtmClass* myRtm;

};
