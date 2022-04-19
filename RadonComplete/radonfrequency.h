#ifndef RADONFREQUENCY_H
#define RADONFREQUENCY_H


#ifndef TRACE
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif


#include "Mute.hpp"
#include "su.h"

void radonsolver0_subtract(float **data, float *h, int nh,float *t, int nt, float dt, float **model, float *q, int nq, float *vel, int itercg, int iter_end, float step, float eps2, float eps1, float quantil, int norm, float factor, float smute, float nmofactor, int rtmethod, float depth, float fmax, char *solver, float **M, int mute);


void radon_param_init(float rtmethod, float fmax, float factor,
                      float aperture, float dx, 
                      float moveoutmin, float moveoutmax, float maxoffset,
                      float& dq, int& nq, float& qmin, float& qmax);

void interpovv (int nt, int ncdp, float *cdp, float **ovv, float cdpt, float *ovvt);

void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv);

/************************************************************************/

void mutemask2(float **M, int nh, int nt, float dt, mutemask_par mutepar);



/************************************************************************/

#endif









