/* 
 * File:   radonclass.cpp
 * Author: dtrad
 * 
 * Created on October 7, 2017, 5:21 PM
 */

#include "radonclass.hpp"
#include "su.h"
void Atimesx(complex *y,complex **A,complex *x,int ny,int nx, int adj);

radonclass::radonclass(float dt, int nt, float moveoutmin, float moveoutmax, float maxoffset, float aperture, int rtmethod, float fmax) {
    myDt=dt;
    myNt=nt;
    myRtmethod=rtmethod;
    myFmax=fmax;
    // convert moveout to sec
    moveoutmin*=0.001;
    moveoutmax*=0.001;
    
    // initial estimates for Radon axis;
    // these parameters will be re-estimated as aperture changes but NQ will be fixed.
    myFactor=1;
    if (myRtmethod ==2){
        aperture=aperture*aperture;
        maxoffset=maxoffset*maxoffset;
    }
    myDq = myFactor*fabs(1./(fmax*aperture));
    myQmin=moveoutmin/maxoffset;
    myQmax=moveoutmax/maxoffset;
    myNq=(myQmax-myQmin)/myDq + 1;
    myQ=ealloc1float(myNq);
    for (int iq=0;iq<myNq;iq++) myQ[iq]=myQmin+myDq*iq;
    myOper=0; // allocated per group (depends on ntraces)
    myG=0;    // allocated per group (depends on ntraces)
    myModel = ealloc2float(myNt,myNq);
   
    
}
radonclass::~radonclass() {
    if (myQ) free1float(myQ);
    if (myOper) free2complex(myOper);
    if (myG)  free1float(myG);
    if (myModel) free2float(myModel);
}

void radonclass::initGroup(float* h, int nh){
    // for now assume dx is constant;
    myNh = nh;
    float dx=fabs(h[1]-h[0]);
    //adjustAxis(h,myNh,dx);
    if (myOper) free2complex(myOper);
    myOper = ealloc2complex(myNq,myNh);
    if (myG) free1float(myG);
    myG = ealloc1float(myNh);
    moveout(h,myG,myNh);
    
}


void radonclass::makeOperator(float w){
  // Operator 
  // Relates the cmp gather and the velocity gather in the f-x space.
    // input is a function of offset (g);
    float *g = myG;
    int nh = myNh;
    complex co;
    complex dco;
    complex phase, dphase;   
    for (int ih=0;ih<nh;ih++){
        phase.r=dphase.r=0;
        phase.i=(w*g[ih]*(myQ[0]-myDq));
        dphase.i=(w*g[ih]*myDq);
        co=exp(phase);
        dco=exp(dphase);
        for (int iq=0;iq<myNq;iq++){
            co*=dco;
            myOper[ih][iq]=conjg(co);                        
        }        
    }  
    return;
}

void radonclass::moveout(float* h, float*g, int nh){
  if (myRtmethod==1) // linear RT
    for (int ih=0;ih<nh;ih++) g[ih]=h[ih];
  else if (myRtmethod==2) // parabolic RT
    for (int ih=0;ih<nh;ih++) g[ih]=h[ih]*h[ih];
//  else if (myRtmethod==3) // pseudo hyperbolic (hyperbolic at one depth)
//    for (int ih=0;ih<nh;ih++) g[ih]=sqrt(h[ih]*h[ih]+depth*depth)-depth;   
    
}

void radonclass::applyOperator(complex* data, complex* model, int adj){
    Atimesx(data,myOper,model,myNh,myNq,adj);           
}

void radonclass::adjustAxis(float* h, int nh, float dx){
    // given offset axis it build the qaxis
    // define a factor for adjusting axis
    float xmin=h[0];
    float xmax=h[nh-1];
    
    if (myRtmethod==1) myDq=1./(myFmax*fabs(xmax-xmin));
    else if (myRtmethod==2)  myDq = 1./(myFmax*(xmax-xmin)*(xmax-xmin));
    
    myDq=myFactor*fabs(myDq);
    myQmax=myQmin+myDq*(myNq-1);
    for (int iq=0;iq<myNq;iq++) myQ[iq]=myQmin+myDq*iq;
    
    
}
complex** radonclass::getOperator(){return myOper;}
float*    radonclass::getMoveoutFunction(){return myG;}
float**   radonclass::getModel(){ return myModel;}
