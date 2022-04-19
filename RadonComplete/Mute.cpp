#include <iostream>
#include <cstdio>
#include <cstdlib>
#include "Mute.hpp"
#include "su.h"

Mute::Mute(mutemask_par mutepar, int moveoutMin, int moveoutMax){
  myTmin = mutepar.tmin;
  myTmax = mutepar.tmax;
  
  myMoveoutMin = moveoutMin*1e-3; // convert to secs
  myMoveoutMax = moveoutMax*1e-3;
  mySlope = mutepar.slope;
  myThreshold = mutepar.threshold;
}

Mute::~Mute(){
  if (myMask) free2float(myMask);
}

void Mute::write_curve_mute(){
  
  FILE *fp;
  fp=efopen("curve1","w");
  
  fprintf(fp,"%f %i\n",myTmin, myIhmin);
  fprintf(fp,"%f %i\n",myTmin, myIhmax);
  fprintf(fp,"%f %i\n",myTmax, myIhmax);
  fprintf(fp,"%f %i\n",myTmax, myIhmin);
  fprintf(fp,"%f %i\n",myTmin, myIhmin);

  efclose(fp);
  
  return;
}

void Mute::init(int maxoffset, int rtmethod, float moveoutMin, 
		float dq, int nq, float dt, int nt){
  myNh = nq;
  myNt = nt;
  myDh = dq;
  myDt = dt;
  myConversion = maxoffset;
  if (rtmethod == 2) myConversion *= maxoffset;
  myMoveoutMin -=moveoutMin;
  myMoveoutMax -=moveoutMin;
  myIhmin = (myMoveoutMin/(myConversion*myDh));
  myIhmax = (myMoveoutMax/(myConversion*myDh));
  
  int itmin = (int) (myTmin/myDt);
  int itmax = (int) (myTmax/myDt);
  int shift;
  myMask = ealloc2float(myNt,myNh);
  for (int ih=0;ih<myNh;ih++){
    for (int it=0;it<myNt;it++){
       if ((it>itmin)&&(it<itmax)){
	shift=(int) ((it-itmin)*myDt*mySlope);
	if ((ih>myIhmin+shift)&&(ih<myIhmax+shift)) myMask[ih][it]=1;
       }
       else myMask[ih][it]=0;
    }
  }
  
  return;
  
}     

float** Mute::getMask(){ return myMask;}
float Mute::getConversion(){ return myConversion;}

