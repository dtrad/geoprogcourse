#include "fftclass.hh"
#include <cstring> // Need for hp memcpy and memset
#include <stdexcept>

using namespace std;

FFTclass::FFTclass(int ntgiven, float dtgiven, vector<float>& bandpass, bool prime) {
    nt = ntgiven;
    dt = dtgiven;

    if (prime) nfft = npfaro(nt, 2 * nt); // transform length returns nt <= nfft <= 2*nt
    else {
        nfft = 2;
        while (nfft < nt) nfft *= 2;
    }



    nf = nfft / 2 + 1; // number of freq;
    if (nf > nt)
        throw runtime_error("FFTclass: nf should be less than nt ");

    onfft = 1. / nfft; // scale factor for ifft;
    df = 1. / (nfft * dt);
    myBandPass = bandpass;

    nmaxfreq = (int) (myBandPass[3] / df + 0.5); // highcut
    nhighfreq = (int) (myBandPass[2] / df + 0.5); // highpass
    nlowpfreq = (int) (myBandPass[1] / df + 0.5); // lowpass
    nminfreq = (int) (myBandPass[0] / df + 0.5); // lowcut
    datafreq = new complex[nf];
    datareal = new float[nfft];
    wfreq = new float[nf];

    myForward = fftwf_plan_dft_r2c_1d(nfft, datareal, (fftwf_complex*) datafreq, FFTW_MEASURE);
    myInverse = fftwf_plan_dft_c2r_1d(nfft, (fftwf_complex*) datafreq, datareal, FFTW_MEASURE);
    int it;
    for (it = 0; it < nlowpfreq; it++) wfreq[it] = 1 - freqweight(it, df, myBandPass[0], myBandPass[1]);
    for (it = nlowpfreq; it < nhighfreq; it++) wfreq[it] = 1.;
    for (it = nhighfreq; it < nmaxfreq; it++) wfreq[it] = freqweight(it, df, myBandPass[2], myBandPass[3]);
    for (it = nmaxfreq; it < nf; it++) wfreq[it] = 0;

}

FFTclass::~FFTclass() {
    fftwf_destroy_plan(myForward);
    fftwf_destroy_plan(myInverse);

    if (datareal) delete[] datareal;
    datareal = 0;
    if (datafreq) delete[] datafreq;
    datafreq = 0;
    if (wfreq) delete[] wfreq;
    wfreq = 0;
}

void FFTclass::display() {
    for (int i = 0; i < 4; i++) cout << "FFT bandpass=" << myBandPass[i] << ", ";
    cout << endl;
    cout << "FFTw -> ";
    fprintf(stdout, " nt = %d, nfft=%d, nf=%d, onfft=%f, dt=%f, df=%f, maxfreqcut=%f, nmaxfreq=%d\n",
            nt, nfft, nf, onfft, dt, df, myBandPass[3], nmaxfreq);

}

void FFTclass::displayWeight() {
    cout << "FFT weight[0]=" << wfreq[0] << " freq = 0 Hz " << endl;
    for (int it = 1; it < nf; it++)
        if (wfreq[it] != wfreq[it - 1])
            cout << "FFT weight[" << it << "]=" << wfreq[it] << " freq = " << it * df << " Hz " << endl;
    cout << "FFT weight[0]=" << wfreq[nf - 1] << " freq = " << nf * df << " Hz " << endl;
}

complex* FFTclass::fftrc() {
    fftwf_execute(myForward);
    return datafreq;
}

float* FFTclass::ifftcr() {
    fftwf_execute(myInverse);
    for (int it = 0; it < nt; it++) datareal[it] *= onfft;
    return datareal;
}

void FFTclass::taper() {
    for (int iw = 0; iw < nf; iw++) datafreq[iw] *= wfreq[iw];
}

bool FFTclass::copyData(float* data, int t0, int nt) {
    memset((void*) datareal, (int) '\0', nfft * sizeof (float));
    for (int it = 0; it < nt; it++) datareal[it] = data[it];
    return true;
}

// copy to workspace if input is in trace order

bool FFTclass::copyData(complex* data) {
    memcpy(datafreq, data, nfft * sizeof (complex));
    return true;
}

// copy to workspace if input is in frequency order

bool FFTclass::copyData(complex** data, int ix) {
    for (int ifreq = 0; ifreq < nf; ifreq++) {
        datafreq[ifreq] = data[ifreq][ix];
    }
    return true;
}

bool FFTclass::getDataFreq(complex **data, int ix){
    for (int ifreq = 0; ifreq < nf; ifreq++) {
        data[ifreq][ix]=datafreq[ifreq];
    }
    return true;
}
bool FFTclass::sinc(float* data, int nt, int dt, float* datafine, int ntf, float dtf) {
    // not implemented yet
    cout << "function to calculate sinc interpolation in Fourier domain " << endl;

    return true;

}

float FFTclass::freqweight(int j, float df, float f1, float f2) {
    // return weight for each frequency
    // int j           freq index
    // float df        freq increment
    // float f1        taper off freq
    // float f2        freq beyond which all components are zero
    float w;
    float f = j*df;
    if (f <= f1) return (1.);
    if (f >= f2) return (0.);
    w = (f2 - f) / (f2 - f1);
    return (w);
}

int FFTclass::npfaro(int nmin, int nmax) {
    throw runtime_error("not implemented yet. Need to port from Seismic Unix");

}

