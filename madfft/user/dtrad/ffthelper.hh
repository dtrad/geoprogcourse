/*
 * Helper class for FFT of a group of traces using fftclass for each trace
 */

/* 
 * File:   ffthelper.hh
 * Author: dtrad
 *
 * Created on February 20, 2019, 5:39 PM
 */

#ifndef FFTHELPER_HH
#define FFTHELPER_HH
using namespace std;
#include <valarray>
#include "fftclass.hh"

class FFThelper {
public:
    FFThelper(int n1, int n2, float d1, vector<float>& bandpass);
    void apply(valarray<float>& in,int n2, int adj);
    void filter();
    void resize(int n2);
    ~FFThelper();
private:
    FFTclass* myFFT;
    int myNt;
    int myNx;
    int myNFreq;
    float myDt;
    complex** myDataFreq;
    valarray<float> myBandPass;
    
    
};

#endif /* FFTHELPER_HH */

