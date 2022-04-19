#ifndef GEOPROG_H
#define GEOPROG_H

#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)

#define TRACE2 { \
  fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__); \
  fflush(stderr); \
}
#define TRACEF(a) { \
  fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__); \
  fprintf a; \
  fflush(stderr); \
}

#include "su.h"
#include "cwp.h"
#include "segy.h"
//#include "Complex.h"
//#include "Dcomplex.h"

segy cleansegy(segy tr);
complex cdot(complex *x,complex *y,int n);
float rcdot(int n, complex *a, complex *b);
float dot(int n, float *a, float *b);
void conv(int lx, int ifx, float *x, int ly, int ify, float *y,
	  int lz, int ifz, float *z);
#endif











