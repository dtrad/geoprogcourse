#ifndef MUTE_HPP
#define MUTE_HPP

using namespace std;

typedef struct {	/* Parameters for inversion  */
  float tmin; 
  float tmax;
  int ihmin;
  int ihmax;
  float slope;
  float threshold;
} mutemask_par;


class Mute{
public:
  Mute(mutemask_par par, int moveoutMin, int moveoutMax);
  ~Mute();
  void write_curve_mute();
  void init(int maxoffset, int rtmethod, float moveoutMin, 
	    float dq, int nq, float dt, int nt);
  float** getMask();
  float   getConversion();
private:
  float myTmin;
  float myTmax;
  float myMoveoutMin;
  float myMoveoutMax;
  float mySlope;
  float myThreshold;
  int  myIhmin;
  int  myIhmax;
  float myConversion;
  int  myNh;
  int  myNt;
  float myDh;
  float myDt;
  float** myMask;




};
#endif
