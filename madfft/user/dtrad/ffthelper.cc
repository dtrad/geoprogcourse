/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ffthelper.cc
 * Author: dtrad
 * 
 * Created on February 20, 2019, 5:39 PM
 */
#include <cstring> // need memcpy and memset
#include "ffthelper.hh"
#include "MyAlloc.hpp"

FFThelper::FFThelper(int n1, int n2, float d1, vector<float>& bandpass) {
    myNt=n1;
    myNx=n2;
    myDt=d1;
    myFFT=new FFTclass(n1,d1,bandpass,false);
    myNFreq=myFFT->nf;
    myDataFreq=0;new2d(myDataFreq,myNFreq,myNx);
    memset(myDataFreq[0],0,myNFreq*n2*sizeof(complex));
    
}
void FFThelper::resize(int n2){
    if (n2==myNx) return;
    else myNx=n2;
    if (myDataFreq) del2d(myDataFreq);
    new2d(myDataFreq,myNFreq,myNx);
    return;
}

FFThelper::~FFThelper() {
    if (myDataFreq) del2d(myDataFreq);
    delete myFFT;
}



void FFThelper::apply(valarray<float>& data,int n2, int adj){
    if (adj){
        if (myNx<n2) resize(n2);        
        for (int i2=0;i2<myNx;i2++){
            myFFT->copyData(&data[i2*myNt],0,myNt);
            complex* tracec = myFFT->fftrc();
            for (int ifreq = 0;ifreq < myNFreq;ifreq++){
                myDataFreq[ifreq][i2] = tracec[ifreq];
            }
        }
    }
    else    
        for (int i2=0;i2<myNx;i2++){
            myFFT->copyData(myDataFreq,i2);
            float* trace = myFFT->ifftcr();                           
            memcpy(&data[i2*myNt],trace,myNt*sizeof(float));
        }                    
}
void FFThelper::filter(){
    for (int i2=0;i2<myNx;i2++){
        myFFT->copyData(myDataFreq,i2);
        myFFT->taper();
        myFFT->getDataFreq(myDataFreq,i2);
        
    }
}
