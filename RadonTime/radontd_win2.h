#include "inversion_par.h"

#ifndef RADONTD_WIN_H
#define RADONTD_WIN_H

#ifndef TRACE
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif


/* FUNCTION PROTOTYPES */

void radontd_sparse(float *t, float *q, float *h, float **m,float **d, int nt, int nh, int nq, float dt, float *vel, float dperv, float pervmin, float t0, inv_par inv, int centralq, int filtout, int nw, float fpeak, int typewav, int LI, float parmute, int mute, int itm, int plot);

void radontd_sparse_win(float *t, float *q, float *h, float **m,float **d, int nt, int nh, int nq, float dt, float **vgrid, inv_par inv, int dataprec, int nw, float fpeak, int typewav, int LI, float parmute, int mute, int itm, int plot);

void build_index_slowness(float *t, float *h, float *q, float **vel, int nt, int nh, int nq,unsigned int **index);

void irreg_slowness_axis(int nq, int nt, float pervmin, float dperv, float *t, float *vel,float *q, float **vgrid, int centralq);

float wpcgnr(void (*oper) (float *,float *,unsigned int **, int, int ,int ,int,int), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wd, float *M, unsigned int **index, inv_par inv);

float wpcgnr(void (*oper) (float *,float *,unsigned int **, int, int ,int ,int,int, float *wavelet, int nw), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wd, float *M, unsigned int **index, inv_par inv, float *wavelet, int nw);


void plotgather(float **d, int nh, int nt, float dt);

int read_ascii_file(const char *name,float *x);

void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv);

void save_gather(float **d, int nh, int nt, float dt, char* name);

void save2dfile(float **d, int nh, int nt, float dt, const char *s);
void save2dfilepower(float **d, int nh, int nt, float dt, float power, const char *s);
void interpovv (int nt, int ncdp, float *cdp, float **ovv, float cdpt, float *ovvt);

//void radonhyp(float *m,float *t, float *h, float *q, float *d, float **vel,int adj,int nt, int nh, int nq);
void radonhyp(float *m, float *d, unsigned int **index, int adj, int nt, int nh, int nq, int nsparse, float *wavelet, int nw);
void radonhyp(float *m, float *d, unsigned int **index, int adj, int nt, int nh, int nq, int nsparse);
void radonhyp_sinc(float *m,float *t, float *h, float *q, float *d, float **vel,int adj,int nt, int nh, int nq);

void subtract_multiples(void  (*oper) (float *, float *, unsigned int **, int , int , int, 
				   int, int, float *, int), 
		    float **data, float **model, unsigned int **index, int adj, int nt, 
		    int nh, int nq, int nsparse, float *wavelet, int nw, float parmute, 
		    float *q, float *h, int itm, int plot);


void subtract_multiples(void  (*oper) (float *, float *, unsigned int **, int , int , int, 
				   int, int), 
		    float **data, float **model, unsigned int **index, int adj, int nt, 
		    int nh, int nq, int nsparse, float parmute, float *q, float *h, int itm,
		    int plot);

/********************* Protoypes for sinc **************************/

void radontd_sparse(float *t, float *q, float *h, float **m,float **d, int nt, int nh, int nq, float dt, float *vel, float dperv, float pervmin, float t0, inv_par inv, int centralq, int filtout, int nw, float fpeak, int typewav, int LI);

void radonhyp(float *m, float *d, float **index, int adj, int nt, int nh, int nq, int ns);

void build_index_slowness(float *t, float *h, float *q, float **slow2, int nt, int nh, int nq,
			  float **index);

float wpcgnr(void (*oper) (float *, float *,float **, int, int ,int ,int,int), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wd, float *M, float **index, inv_par inv);

float testadjop(void (*oper) (float *,float *,float **,int ,int ,int, int, int),float **index,int nt, int nh, int nq, int nsparse);

float wpcgnr_m0(void (*oper) (float *,float *,unsigned int **, int, int ,int ,int,int), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wd, float *M, unsigned int **index, inv_par inv);

/****************************************************************************************/

/********************* Protoypes for LI **************************/


void build_index_slowness_li(float *t, float *h, float *q, float  **vel, int nt, int nh, int nq, unsigned int **index);


void radonhyp_li(float *m, float *d, unsigned int **index, int adj, int nt, int nh, int nq, int ns);

void radonhyp_li(float *m, float *d, unsigned int **index, int adj, int nt, int nh, int nq, int ns, float *wavelet, int nw);

/****************************************************************************************/


void contruc_2(int conj,int add, int nx, float *xx, int nb, float *bb,float *yy);

int get_wavelet(float *wavelet,const char *name,int nw, int type, float dt, float fpeak);

/***********************Protoptypes fo radontd_sparse_interp ****************************/

void radontd_sparse(float *t, float *q, float *h, float **m,float **d, int nt, int nh, int nq, float dt, float *vel, float dperv, float pervmin, float t0,  inv_par inv, int centralq, int filtout, int nw, float fpeak, int typewav, int LI);

float wpcgnr(void (*oper) (float *,float *,float **, int, int ,int ,int,int, float *wavelet, int nw), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wd, float *M, float **index, inv_par inv, float *wavelet, int nw);


float testadjop(void (*oper) (float *,float *,unsigned int **,int, int ,int, int, int),unsigned int **index,int nt, int nh, int nq, int nsparse);

float testadjop(void (*oper) (float *,float *,unsigned int **,int ,int ,int, int, int, float *wavelet, int nw),unsigned int **index,int nt, int nh, int nq, int nsparse, float *wavelet,int nw);

float testadjop_conv(void (*oper) (int, int, int, float *,int, float *,float *),int nb, float *bb, int nx, int ny);

void radonhyp(float *m, float *d, float **index, int adj, int nt, int nh, int nq, int nsparse);

void radonhyp(float *m, float *d, float **index, int adj, int nt, int nh, int nq, int nsparse, float *wavelet, int nw);

/****************************************************************************************/
void smute_gather(float **d, int nt, int nh, float *t, float *h, float *vel, float smute);

void weights_td(float *m, int nx, int norm, float sigmam, float *Wm, int iter);
void filterOutlier(int dataprec, float* h, int nh, float* t, int nt, float dt, float** d, float** Wd); 
void deviations(float *m, int nx, float *d, int ny, int norm, float quantil1, float quantil2, float *sigmam, float *sigmad);
int taper(float **data, int nt, int nh, int ntaper,int flag);
float dot(int n, float *a, float *b);
void save_gather(float **d, int nh, int nt, float dt, char* name);
void save_gather(float **d, int nh, float *h, int nt, float dt, char* name);

#endif







