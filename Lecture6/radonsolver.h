/* 
 * File:   radonsolver.h
 * Author: dtrad
 *
 * Created on October 8, 2017, 8:20 AM
 */

#ifndef RADONSOLVER_H
#define	RADONSOLVER_H
#include "radonclass.hpp"

void radonsolver(radonclass* rad, float** data, char* solver, float sigmad);
void radoninv(radonclass* rad, float** data);
void fft_parameters(int nt, float dt, int *pnfft, int *pnf, float *pdf);
void fftgo_xt2fx(int sign,float **d,complex  **m, int nh, int nt, float dt, int nfft, int nf);
int fftback_fx2xt(int sign,float **d,complex  **m, int nh, int nt, float dt, int nfft, int nf);
float freqweight(int j, float df, float f1, float f2);
float radon_toeplitz(complex *d, complex **L, complex *m, float sigmad, int nh, int nq);

#endif	/* RADONSOLVER_H */

