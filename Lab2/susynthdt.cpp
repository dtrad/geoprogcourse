/* Copyright (c) University of British Columbia, 1999.*/
/* SUSYNTH2:  $Date: December 1999  */

#include "su.h"
#include "segy.h"

#define TRACE { \
  fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__); \
  fflush(stderr); \
}

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUSYNTHDT    > output                                               ",
  " Required parameters:                                                ",
  "                                                                     ",
  "                                                                     ",
  "                                                                     ",
  "                                                                     ",
  " Optional parameters                                                 ",
  "                                                                     ",
  "  Example 1                                                          ",
  NULL};

/* Credits:
 *	
/**************** end self doc ***********************************/

int main(int argc, char **argv) {
    int nx;
    float dt;
    segy tr;

    // Initialize 

    initargs(argc, argv);
    requestdoc(0);
    fprintf(stderr, "argc=%d\n", argc);

    if (!getparfloat("dt", &dt)) dt = 0.004;
    if (!getparint("nx", &nx)) nx = 30;
    
    for (int ix=0;ix<nx;ix++){


      tr.dt = (unsigned short) (dt * 1e6);
      tr.ntr = (int) nx;
      puttr(&tr);
    }
    


    return EXIT_SUCCESS;
}














