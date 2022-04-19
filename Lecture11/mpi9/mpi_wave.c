/* ------------------------------------------------------------------------
 * MPI Concurrent Wave Equation - C Version
 * Point-to-Point Communications Example
 *
 * FILE: mpi_wave.c
 * OTHER FILES: draw_wave.c
 * DESCRIPTION:
 *   This program implements the concurrent wave equation described 
 *   in Chapter 5 of Fox et al., 1988, Solving Problems on Concurrent
 *   Processors, vol 1.  
 *
 *   A vibrating string is decomposed into points.  Each processor is 
 *   responsible for updating the amplitude of a number of points over
 *   time.
 *	
 *   At each iteration, each processor exchanges boundary points with
 *   nearest neighbors.  This version uses low level sends and receives
 *   to exchange boundary points.
 *
 * AUTHOR: Blaise Barney - adapted from original by Roslyn Leibensperger
 *   Converted to MPI: George L. Gusciora (1/25/95)
 * LAST REVISED: 10/14/98 Blaise Barney
 * ------------------------------------------------------------------------ */
#include "mpi.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

#define MASTER 0
#define MAXPROC 8
#define MAXPOINTS 1000
#define MAXSTEPS 10000
#define PI 3.14159265

int E_RtoL = 10;
int E_LtoR = 20;
int E_OUT1 = 30;
int E_OUT2 = 40;

void init_master(void);
void init_workers(void);
void init_line(void);
void update (int left, int right);
void output_master(void);
void output_workers(void);
extern void draw_wave(double *);

int	taskid,               /* task ID */
	numtasks,             /* number of processes */
	nsteps,               /* number of time steps */
	tpoints,              /* total points along string */
	npoints,              /* number of points handled by this processor */
	first,                /* index of 1st point handled by this processor */
	rcode;                /* generic return code */
double	etime,                /* elapsed time in seconds */
	values[MAXPOINTS+2],  /* values at time t */
	oldval[MAXPOINTS+2],  /* values at time (t-dt) */
	newval[MAXPOINTS+2];  /* values at time (t+dt) */
/*  ------------------------------------------------------------------------
 *	Master obtains input values from user
 *  ------------------------------------------------------------------------ */
void init_master(void)
{
   char tchar[8];
   int buffer[2];

   /* set number of points, number of iterations */
   tpoints = 0;
   nsteps = 0;
   while ((tpoints < numtasks) || (tpoints > MAXPOINTS))
   {
      printf("Enter number of points along vibrating string\n");
      scanf("%s", tchar);
      tpoints = atoi(tchar);
      if ((tpoints < numtasks) || (tpoints > MAXPOINTS))
         printf("enter value between %d and %d\n", numtasks, MAXPOINTS);
   }
   while ((nsteps < 1) || (nsteps > MAXSTEPS))
   {
      printf("Enter number of time steps\n");
      scanf("%s", tchar);
      nsteps = atoi(tchar);
      if ((nsteps < 1) || (nsteps > MAXSTEPS))
         printf("enter value between 1 and %d\n", MAXSTEPS);
   }

   printf("%d: points = %d, steps = %d\n", taskid, tpoints, nsteps);
   /* then broadcast time advance parameter, total points, time steps */
   buffer[0] = tpoints;
   buffer[1] = nsteps;
   MPI_Bcast(buffer, 2, MPI_INT, MASTER, MPI_COMM_WORLD);
}

/*  -------------------------------------------------------------------------
 *	Workers receive input values from master
 *  -------------------------------------------------------------------------*/
void init_workers(void)
{
   int buffer[2];

   /* receive total points, time steps */
   MPI_Bcast(buffer, 2, MPI_INT, MASTER, MPI_COMM_WORLD);
   tpoints = buffer[0];
   nsteps = buffer[1];
}

/*  ------------------------------------------------------------------------
 *     All processes initialize points on line
 *  --------------------------------------------------------------------- */
void init_line(void)
{
   int nmin, nleft, npts, i, j, k;
   double x, fac;

   /* calculate initial values based on sine curve */
   nmin = tpoints/numtasks;
   nleft = tpoints%numtasks;
   fac = 2.0 * PI;
   for (i = 0, k = 0; i < numtasks; i++)
   {
      npts = (i < nleft) ? nmin + 1 : nmin;
      if (taskid == i)
      {
         first = k + 1;
         npoints = npts;
         printf ("%d: first = %d, npoints = %d\n", taskid, first, npts);
         for (j = 1; j <= npts; j++, k++)
         {
            x = (double)k/(double)(tpoints - 1);
            values[j] = sin (fac * x);
         } 
      }
      else k += npts;
   }
   for (i = 1; i <= npoints; i++) oldval[i] = values[i];
}

/*  -------------------------------------------------------------------------
 *	All processes update their points a specified number of times 
 *  -------------------------------------------------------------------------*/
void update(int left, int right)
{
   int i, j, id_rtol, id_ltor, nbytes;
   clock_t start;
   double dtime, c, dx, tau, sqtau;
   MPI_Request request;
   MPI_Status status;

   dtime = 0.3;
   c = 1.0;
   dx = 1.0;
   tau = (c * dtime / dx);
   sqtau = tau * tau;

   start = clock();
   /* update values for each point along string */
   for (i = 1; i <= nsteps; i++)
   {
      /* exchange data with "left-hand" neighbor */
      if (first != 1) 
      {
         MPI_Isend(&values[1], 1, MPI_DOUBLE, left, E_RtoL, MPI_COMM_WORLD,
                   &request);
         MPI_Recv(&values[0], 1, MPI_DOUBLE, left, E_LtoR, MPI_COMM_WORLD, &status);
         MPI_Wait(&request, &status);
      }
      /* exchange data with "right-hand" neighbor */
      if (first + npoints -1 != tpoints)
      {
         MPI_Isend(&values[npoints], 1, MPI_DOUBLE, right, E_LtoR,
                    MPI_COMM_WORLD, &request);
         MPI_Recv(&values[npoints+1], 1, MPI_DOUBLE, right, E_RtoL,
                   MPI_COMM_WORLD, &status);
         MPI_Wait(&request, &status);
      }
      /* update points along line */
      for (j = 1; j <= npoints; j++)
      {
         /* global endpoints */
         if ((first + j - 1 == 1) || (first + j - 1 == tpoints))
            newval[j] = 0.0;
         else
            /* use wave equation to update points */
            newval[j] = (2.0 * values[j]) - oldval[j]
                      + (sqtau * (values[j-1] - (2.0 * values[j]) + values[j+1]));
      }
      for (j = 1; j <= npoints; j++)
      {
         oldval[j] = values[j];
         values[j] = newval[j];
      }
   }
   etime = (double) (clock() - start);
}

/*  ------------------------------------------------------------------------
 *	Master receives results from workers and prints
 *  ------------------------------------------------------------------------ */
void output_master(void) 
{
   int i, source, nbytes, start, npts, tpts, buffer[2];
   double results[MAXPOINTS], tarray[MAXPROC];
   double tmin = 10000.0;
   double tmax = 0.0;
   double tsum = 0.0;
   MPI_Status status;
 
   /* get values DONTCARE and ALLGRP first */
   /* store worker's results in results array */
   for (i = 1; i < numtasks; i++)
   {
      /* receive number of points and first point */
      MPI_Recv(buffer, 2, MPI_INT, MPI_ANY_SOURCE, E_OUT1,
               MPI_COMM_WORLD, &status);
      start = buffer[0];
      npts = buffer[1];
      /* receive results */
      MPI_Recv(&results[start-1], npts, MPI_DOUBLE, MPI_ANY_SOURCE, E_OUT2,
                MPI_COMM_WORLD, &status);
   }

   /* store master's results in results array */
   for (i = first; i < first + npoints; i++)
      results[i-1] = values[i];

   tpts = (tpoints < 10) ? tpoints: 10; 
   printf("first %d points (for validation):\n", tpts);
   for (i = 0; i < tpts; i++) printf("%4.2f  ", results[i]);
   printf("\n");

   /* gather elapsed time from workers */
   MPI_Gather(&etime, 1, MPI_DOUBLE,
               tarray, 1, MPI_DOUBLE,
               MASTER, MPI_COMM_WORLD);
   for (i = 0; i < numtasks; i++)
   {
      tarray[i] = tarray[i]/1.0e6;
      tsum += tarray[i];
      if (tarray[i] < tmin) tmin = tarray[i];
      if (tarray[i] > tmax) tmax = tarray[i];
   }
   printf("\n");
   printf("Minimum, maximum, and average time to update wave = \n");
   printf("%5.3f, %5.3f, %5.3f seconds\n\n", tmin, tmax, tsum/numtasks);

   /* display results with draw_wave routine */
   draw_wave(&results[0]);
}

/*  -------------------------------------------------------------------------
 *	Workers send the updated values to the master
 *  -------------------------------------------------------------------------*/
 
void output_workers(void)
{
   int buffer[2], nbytes, wait;
   double tarray[MAXPROC];
   MPI_Request request1, request2;
   MPI_Status status;

   /* send first point and number of points handled to master */
   buffer[0] = first;
   buffer[1] = npoints;
   MPI_Isend(&buffer, 2, MPI_INT, MASTER, E_OUT1,
              MPI_COMM_WORLD, &request1);
   /* send results to master */
   MPI_Isend(&values[1], npoints, MPI_DOUBLE, MASTER, E_OUT2,
              MPI_COMM_WORLD, &request2);
   MPI_Wait(&request1, &status);
   MPI_Wait(&request2, &status);

   /* send elapsed time to master*/
   MPI_Gather(&etime, 1, MPI_DOUBLE,
              tarray, 1, MPI_DOUBLE,
              MASTER, MPI_COMM_WORLD);
}

/*  ------------------------------------------------------------------------
 *	Main program
 *  ------------------------------------------------------------------------ */

int main(argc,argv)
int argc;
char *argv[];
{
   int left, right, rc;

   /* enroll process */
   rc = MPI_Init(&argc,&argv);
   rc|= MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
   rc|= MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
   if (rc != MPI_SUCCESS)
      printf ("error initializing MPI and obtaining task ID information\n");
   else
      printf ("task ID = %d\n", taskid);

   /* determine left and right neighbors */
   if (taskid == numtasks-1)
      right = 0;
   else
      right = taskid + 1;

   if (taskid == 0)
      left = numtasks - 1;
   else
      left = taskid - 1;

   /* get program parameters and initialize wave values */
   if (taskid == MASTER)
      init_master();
   else
      init_workers();

   init_line();

   /* update values along the line for nstep time steps */
   update(left, right);

   /* master collects results from workers and prints */
   if (taskid == MASTER)
      output_master();
   else
      output_workers();

   MPI_Finalize();

   return 0;
}

