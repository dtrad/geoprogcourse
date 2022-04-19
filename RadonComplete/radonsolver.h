#ifndef RADONSOLVER_H
#define RADONSOLVER_H

#ifndef TRACE
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif



#include"su.h"
#include"dan.h"
#include "inversion_par.h"


void radonsolver(float *pos, int nh, float **data, float *t, int nt, float dt,
                 float **model, float *q, int nq, float dq, float eps1, float eps2, 
                 float eps, float fmax, float *Wd, int itercg, int iter_end, 
                 int norm,float step, int testadj, int rtmethod, float depth, 
                  char *solver);

void radon_matrix(complex **l, float *h,float *q,int nh,int nq,float w);

void radon_moveout(float *h, float *g, int nh, int rtmethod, float depth);

void hrrti(float **d,float *pos,float dt,float **m,float *q, float fmax, 
	   int nt, int  nh, int nq, int rtmethod, float depth);

int taper(float **data, int nt, int nh, int ntaper,int flag);

int taper(float *data, int nh, int ntaper);

void nmo(float *d,float *m,float *t,float h,float *vel,int adj,int nt,float dt,
	 float smute);

void weights_inv(complex *m, int nx, int norm, float sigmam, float *Wm, int iter);

void deviations(complex *m, int nx, complex *d, int ny, int norm, float quantil1, 
		float quantil2, float *sigmam, float *sigmad);

void dataweigths(float *pos, int nh, float *Wd, int add);

void radon_matrix_cgfft(complex **L,complex *RC, int nh, int nq, int nf2, float *Cd);

float wtcgls(complex *b,complex **L, complex *x,float *Wm,
	     float *Wd,int nh,int nq, float tol, float step, int itercg);

float radon_cgfft(complex *d2, complex **L, complex *RC, complex *m2, int nh, int nq, int nf2, float *Wm, float *Wd, float eps, int itercg, float step);

float radon_toeplitz(complex *d, complex **L, complex *m, float sigmad, int nh, int nq);

float radon_cgtoep(complex *d, complex **L, complex *m, float sigmad, int nh, int nq);

void weights_window_inv(complex **m, int buffer, int nq, int freq, int norm, float sigmam, float *Wm, int iter);


#endif









