#ifndef INVERSION_PAR_H
#define INVERSION_PAR_H

#define MAXITER 100
/* TYPEDEFS */
typedef struct {	/* Parameters for inversion  */

  int itercg;	/* Number maximum of iterations in CG*/

  int iter_end;	/* Number maximum of External Iterations */

  float eps1;	/* noise covariance    */

  float eps2 ;	/* model covariance */

  float eps;	/* small number for tolerance */

  float step;	/* step length factor  */

  int norm;     /* norm to use in model weights; */

  int restart;  /* always set to 1 for now */

  int taperflag; /* If set applies taper to 5 outer traces */

  int mute;      /* if set applies mute */

  float parmute; /* if mute is set it controls the length of muting */

  float J[MAXITER];
  
} inv_par;


typedef struct {	/* Parameters for inversion  */

  int taperflag; /* If set applies taper to 5 outer traces */

  int mute;      /* if set applies mute */

  float parmute; /* if mute is set it controls the length of muting */

  int testadj;   /* if set test adjoint */

} options_par;


#endif
























