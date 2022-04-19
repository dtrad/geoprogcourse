#ifndef FFTclass_HPP
#define FFTclass_HPP
#include <iostream>
#include <stdio.h>
#include "complex.hh"
#include <fftw3.h>
#include <vector>


using namespace std;
#define complex MYCOMPLEX::complex
class FFTclass {
public:
    FFTclass(int nt, float dt, vector<float>& bandpass, bool pnumber = false);
    ~FFTclass();
    void display();
    complex* fftrc();
    float* ifftcr();
    bool copyData(float* data, int t0, int nt); // copy trace input to time workspace
    bool copyData(complex* data); // copy trace input to freq workspace
    bool copyData(complex** data, int ix); // copy freq order input
    bool getDataFreq(complex **data, int ix);
    bool sinc(float* data, int nt, int dt, float* datafine, int ntf, float dtf);
    void taper();
    void displayWeight();
    int nt;
    int nf;
    int nfft;
    float dt;
    float df;
    float onfft;
    complex* datafreq; // freq workspace (fftw)
    float* datareal; // time worspace  (fftw)
    float* wfreq;
    vector<float> myBandPass;
    float lowCut;
    int nmaxfreq; // highcut
    int nhighfreq; // highpass
    int nlowpfreq; // lowpass
    int nminfreq; // lowcut
private:
    float freqweight(int j, float df, float f1, float f2);
    fftwf_plan myForward;
    fftwf_plan myInverse;
    int npfaro(int nmin, int nmax); // TODO

};

#endif 
