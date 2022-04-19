// rtmborn. 
// Daniel Trad, 2016. 
// Based on Penliang Yang program user/pyang/prtm2d.c

#include <iostream>
#include <cstdlib>
#include <cstring> // need memcpy
#include <cmath>
#include "rtmborn.hh"
#include "MyAlloc.hpp"
#include <time.h>
#include <stdexcept>
#include <sys/time.h>
#include <valarray>

#ifdef _OPENMP
#include <omp.h>
#endif

// #ifdef _OPENMP
// #undef _OPENMP
// #endif

#ifndef MARK
#define MARK fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif


// class for forward modeling with second, fourth and eight order staggered grid.

RtmBorn::RtmBorn(RtmPar* par, float* vel) {
    cerr << " RTM class constructor " << endl;
    myDz = par->dz;
    myDx = par->dx;
    myDt = par->dt;

    myNx = par->nx;
    myNz = par->nz;
    myNt = par->nt;

    myNg = par->ng;
    myNs = par->ns;

    myFm = par->fm;
    myAmp = par->amp;

    myPI = acos(-1.);
    myNb = par->nb;
    myNbUp = par->nbup; // nb at surface. =0 to produce multiples.
    myNIter = par->niter;

    mySP0 = 0;
    mySP1 = 0;
    myGP0 = 0;
    myGP1 = 0;
    mySPtt = 0;

    myDObs = 0;
    myTrans = 0;
    myVV = 0;
    myV0 = 0;
    myVel = 0;

    myBndr = 0;
    myBndrUp = 0;
    myWlt = 0;
    mySxz = 0;
    myGxz = 0;
    myModTmp = 0;
    myVerbose = true;
    myNthreads = 1;
    myChunksize = 1;
    myFromBoundary = 1;
    myIter = 0; // iteration number (used when changing operator with iter)

    myEightOrder = false;
    myFourthOrder = false;
    mySecondOrder = false;
    int order = par->order;
    if (order == 2) mySecondOrder = true;
    else if (order == 4) myFourthOrder = true;
    else if (order == 8) {
        myEightOrder = true;
        myFromBoundary = 0;
    } else {
        cerr << "method not implemented = " << order << endl;
        throw;
    }

    cerr << "Using " << order << "th order " << endl;

    // several options for forward modeling
    // 0 correl
    // 1 correl with Laplacian 
    // 2 Born with laplacian for p_tt calculation
    // 3 Born with proper second derivative
    myOption = par->option;
    myApplyLaplac = par->applyLaplac;
    cerr << "Using option = " << myOption << endl;
    if (myApplyLaplac) cerr << "apply Laplacian to image " << endl;
    myMute = par->mute;
    if (myMute) cerr << "Apply mute to the first arrival" << endl;
    myMaxOffset = par->maxoffset;
    if (myMaxOffset) cerr << "Limit to maximum offset of " << myMaxOffset << endl;

    //estimate direct wave velocity
    float tempsum = 0;
    for (int ix = 0; ix < myNx; ix++) tempsum += vel[ix * myNz + 3]; // replace later 3 by source depth
    myDWVel = tempsum / myNx;

#ifdef _OPENMP
    myNthreads = omp_get_max_threads();
    cerr << "Using " << myNthreads << " threads " << endl;
    myChunksize = (myNx + myNb + myNb) / myNthreads / 8;
    cerr << "chunksize" << myChunksize << endl;
#endif


    init(vel);
    if (!isStable(vel)) {
        cerr << "WARNING********UNSTABLE*********" << endl;
        //throw runtime_error("UNSTABLE");
    } else cerr << "STABLE!" << endl;
    cerr << " ********* RTM constructor done ******** " << endl;


}

RtmBorn::~RtmBorn() {
    cerr << " calling RTM destructor " << endl;

    if (myV0) del2d(myV0);
    if (myVV) del2d(myVV);
    if (myVel) del2d(myVel);
    if (mySP0) del2d(mySP0);
    if (mySP1) del2d(mySP1);

    if (myGP0) del2d(myGP0);
    if (myGP1) del2d(myGP1);

    if (myDObs) del1d(myDObs);
    if (mySPtt) del2d(mySPtt); // source second time derivative

    if (myTrans) del1d(myTrans);
    if (myBndr) del1d(myBndr);
    if (myBndrUp) del1d(myBndrUp);
    if (myWlt) del1d(myWlt);
    if (mySxz) del1d(mySxz);
    if (myGxz) del1d(myGxz);


    if (myModTmp) del1d(myModTmp);

    if (myRWBndr) del1d(myRWBndr);
    if (!myFromBoundary) if (mySp0array) del3d(mySp0array);

    cerr << " RTM destructor done " << endl;
}

void RtmBorn::init(float* vel) {

    myNzpad = myNz + myNb + myNbUp;
    myNxpad = myNx + myNb + myNb;


    myC11 = 4. / 3. / (myDz * myDz);
    myC21 = 4. / 3. / (myDx * myDx);
    myC12 = -1. / 12. / (myDz * myDz);
    myC22 = -1. / 12. / (myDx * myDx);
    myC0 = -2. * (myC11 + myC12 + myC21 + myC22);

    ////for eight order
    myEC4xx = -1. / 560. / (myDx * myDx);
    myEC3xx = 8. / 315. / (myDx * myDx);
    myEC2xx = -1. / 5. / (myDx * myDx);
    myEC1xx = 8. / 5. / (myDx * myDx);
    myEC0xx = -205. / 72. / (myDx * myDx);

    myEC4zz = -1. / 560. / (myDz * myDz);
    myEC3zz = 8. / 315. / (myDz * myDz);
    myEC2zz = -1. / 5. / (myDz * myDz);
    myEC1zz = 8. / 5. / (myDz * myDz);
    myEC0zz = -205. / 72. / (myDz * myDz);


    new2d(myV0, myNx, myNz);
    new2d(myVel, myNx, myNz);
    new2d(myVV, myNxpad, myNzpad);
    new2d(mySP0, myNxpad, myNzpad);
    new2d(mySP1, myNxpad, myNzpad);
    new2d(myGP0, myNxpad, myNzpad);
    new2d(myGP1, myNxpad, myNzpad);

    // added for Born 
    new2d(mySPtt, myNxpad, myNzpad);

    if (!myFromBoundary) new3d(mySp0array, myNt, myNxpad, myNzpad);

    new1d(myDObs, myNg * myNt);
    new1d(myTrans, myNg * myNt);
    new1d(myModTmp, myNx * myNz);
    new1d(myWlt, myNt);
    new1d(mySxz, myNs);
    new1d(myGxz, myNg);
    new1d(myRWBndr, myNt * 4 * (myNx + myNz));

    memcpy(myV0[0], vel, myNx * myNz * sizeof (float));
    memcpy(myVel[0], vel, myNx * myNz * sizeof (float));
    memset(myDObs, 0, myNg * myNt * sizeof (float));
    memset(myTrans, 0, myNg * myNt * sizeof (float));

    setBndr(myBndr, myNb);
    setBndr(myBndrUp, myNbUp);

    expand2d(myVV, myV0);
    float tmp;
    for (int ix = 0; ix < myNxpad; ix++) {
        for (int iz = 0; iz < myNzpad; iz++) {
            tmp = myVV[ix][iz] * myDt;
            myVV[ix][iz] = tmp*tmp;
        }
    }

    for (int it = 0; it < myNt; it++) {
        tmp = myPI * myFm * (it * myDt - 1.0 / myFm);
        tmp = tmp*tmp;
        myWlt[it] = myAmp * (1. - 2. * tmp) * expf(-tmp);
    }
}

void RtmBorn::setBndr(float*& bndr, int nb) {
    // initialize sponge ABC coefficients 
    new1d(bndr, nb);
    float t;
    for (int ib = 0; ib < nb; ib++) {
        t = 0.015 * (nb - 1 - ib);
        bndr[ib] = expf(-t * t);
    }
}

float* RtmBorn::getWavelet() {
    return myWlt;
}

void RtmBorn::expand2d(float** b, float** a) {
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for (int ix = 0; ix < myNx; ix++) {
        for (int iz = 0; iz < myNz; iz++) {
            b[myNb + ix][myNbUp + iz] = a[ix][iz];
        }
    }


    for (int ix = 0; ix < myNxpad; ix++) {
        for (int iz = 0; iz < myNbUp; iz++) { // top
            b[ix][ iz] = b[ix][myNbUp ];
        }

        for (int iz = 0; iz < myNb; iz++) {
            b[ix][myNzpad - iz - 1] = b[ix][myNzpad - myNb - 1]; /* bottom*/
        }
    }
    for (int ix = 0; ix < myNb; ix++) {
        for (int iz = 0; iz < myNzpad; iz++) {
            b[ix ][iz] = b[myNb ][iz]; /* left */
            b[myNxpad - ix - 1 ][iz] = b[myNxpad - myNb - 1 ][iz]; /* right */
        }
    }

}

void RtmBorn::window2d(float **a, float **b)
/*< window 'b' to 'a': source(b)-->destination(a) >*/ {
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for (int ix = 0; ix < myNx; ix++) {
        for (int iz = 0; iz < myNz; iz++) {
            a[ix][iz] = b[myNb + ix][myNbUp + iz];
        }
    }
}

void RtmBorn::rw_snapshot(float** p, int it, bool read) {
    // read/write snapshot completely 
    // if read=true read, else write

    if (!read) {
        for (int ix = 0; ix < myNxpad; ix++) {
            for (int iz = 0; iz < myNzpad; iz++) {
                mySp0array[it][ix][iz] = p[ix][iz];
            }
        }
    } else {
        for (int ix = 0; ix < myNxpad; ix++) {
            for (int iz = 0; iz < myNzpad; iz++) {
                p[ix][iz] = mySp0array[it][ix][iz];
            }
        }
    }

}

void RtmBorn::boundaryRW(float** p, int it, bool read) {

    float* spo = &myRWBndr[it * 4 * (myNx + myNz)];
    // read/write using effective boundary saving strategy
    // if read = true, read the boundary out, else save/write boundary
    if (read) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int ix = 0; ix < myNx; ix++) {
            for (int iz = 0; iz < 2; iz++) {
                p[ix + myNb][iz - 2 + myNbUp] = spo[iz + 4 * ix];
                p[ix + myNb][iz + myNz + myNbUp] = spo[iz + 2 + 4 * ix];
            }
        }
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int iz = 0; iz < myNz; iz++) {
            for (int ix = 0; ix < 2; ix++) {
                p[ix - 2 + myNb][iz + myNbUp] = spo[4 * myNx + iz + myNz * ix];
                p[ix + myNx + myNb][iz + myNbUp] = spo[4 * myNx + iz + myNz * (ix + 2)];
            }
        }
    } else {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int ix = 0; ix < myNx; ix++) {
            for (int iz = 0; iz < 2; iz++) {
                spo[iz + 4 * ix] = p[ix + myNb][iz - 2 + myNbUp];
                spo[iz + 2 + 4 * ix] = p[ix + myNb][iz + myNz + myNbUp];
            }
        }
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int iz = 0; iz < myNz; iz++) {
            for (int ix = 0; ix < 2; ix++) {
                spo[4 * myNx + iz + myNz * ix] = p[ix - 2 + myNb][iz + myNbUp];
                spo[4 * myNx + iz + myNz * (ix + 2)] = p[ix + myNx + myNb][iz + myNbUp];
            }
        }
    }

}

void RtmBorn::initGeometry(int xbeg, int zbeg, int jx, int jz, int n, string type) {
    // set geometry ////////////////////////////////////
    if (!(xbeg >= 0 && zbeg >= 0 && xbeg + (n - 1) * jx < myNx && zbeg + (n - 1) * jz < myNz)) {
        cerr << type << " exceeds the computing zone" << endl;
        exit(1);
    }
    sg_init(zbeg, xbeg, jz, jx, type);
    if (type == "rcvrs") myDg = jx * myDx;
}

void RtmBorn::sg_init(int szbeg, int sxbeg, int jsz, int jsx, string type) {
    if (type == "shots") sg_init(mySxz, szbeg, sxbeg, jsz, jsx, myNs);
    else if (type == "rcvrs") sg_init(myGxz, szbeg, sxbeg, jsz, jsx, myNg);
}

void RtmBorn::sg_init(int* sxz, int szbeg, int sxbeg, int jsz, int jsx, int ns) {
    int sx, sz;
    for (int is = 0; is < ns; is++) {
        sz = szbeg + is*jsz;
        sx = sxbeg + is*jsx;
        sxz[is] = myNz * sx + sz;
    }
}

void RtmBorn::stepForward(float** p0, float** p1, bool adj) {
    if (mySecondOrder) stepForwardSO(p0, p1, adj);
    else if (myFourthOrder) stepForwardFO(p0, p1, adj);
    else if (myEightOrder) stepForwardEO(p0, p1, adj);
}

void RtmBorn::stepForwardEO(float** p0, float** p1, bool adj) {
    //intuitive version of 8th order FDM
    if (adj) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int ix = 4; ix < myNxpad - 4; ix++)
            for (int iz = 4; iz < myNzpad - 4; iz++) {
                float tmpx = (
                        myEC4xx * (myVV[ix-4][iz]*p1[ix-4][iz] + myVV[ix+4][iz]*p1[ix+4][iz]) +
                        myEC3xx * (myVV[ix-3][iz]*p1[ix-3][iz] + myVV[ix+3][iz]*p1[ix+3][iz]) +
                        myEC2xx * (myVV[ix-2][iz]*p1[ix-2][iz] + myVV[ix+2][iz]*p1[ix+2][iz]) +
                        myEC1xx * (myVV[ix-1][iz]*p1[ix-1][iz] + myVV[ix+1][iz]*p1[ix+1][iz]) +
                        myEC0xx * myVV[ix][iz]*p1[ix][iz]
                        );

                float tmpz = (
                        myEC4zz * (myVV[ix][iz-4]*p1[ix][iz-4] + myVV[ix][iz+4]*p1[ix][iz+4]) +
                        myEC3zz * (myVV[ix][iz-3]*p1[ix][iz-3] + myVV[ix][iz+3]*p1[ix][iz+3]) +
                        myEC2zz * (myVV[ix][iz-2]*p1[ix][iz-2] + myVV[ix][iz+2]*p1[ix][iz+2]) +
                        myEC1zz * (myVV[ix][iz-1]*p1[ix][iz-1] + myVV[ix][iz+1]*p1[ix][iz+1]) +
                        myEC0zz * myVV[ix][iz]*p1[ix][iz ]
                        );

                p0[ix][iz] = 2 * p1[ix][iz] - p0[ix][iz] + (tmpx + tmpz);
            }
    } else {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int ix = 4; ix < myNxpad - 4; ix++)
            for (int iz = 4; iz < myNzpad - 4; iz++) {
                float tmpx = (
                        myEC4xx * (p1[ix - 4][iz] + p1[ix + 4][iz]) +
                        myEC3xx * (p1[ix - 3][iz] + p1[ix + 3][iz]) +
                        myEC2xx * (p1[ix - 2][iz] + p1[ix + 2][iz]) +
                        myEC1xx * (p1[ix - 1][iz] + p1[ix + 1][iz]) +
                        myEC0xx * p1[ix ][iz]
                        );

                float tmpz = (
                        myEC4zz * (p1[ix][iz - 4] + p1[ix][iz + 4]) +
                        myEC3zz * (p1[ix][iz - 3] + p1[ix][iz + 3]) +
                        myEC2zz * (p1[ix][iz - 2] + p1[ix][iz + 2]) +
                        myEC1zz * (p1[ix][iz - 1] + p1[ix][iz + 1]) +
                        myEC0zz * p1[ix][iz ]
                        );

                p0[ix][iz] = 2 * p1[ix][iz] - p0[ix][iz] + myVV[ix][iz]*(tmpx + tmpz);
            }
    }

}

void RtmBorn::stepForwardFO(float** p0, float** p1, bool adj)
/*< original forward modeling step fourth order>*/ {

    if (adj) {
#ifdef _OPENMP
#pragma omp parallel for 
#endif
        for (int ix = 2; ix < myNxpad - 2; ix++)
            for (int iz = 2; iz < myNzpad - 2; iz++) {
                p0[ix][iz] = 2 * p1[ix][iz] - p0[ix][iz]
                        + myC0 * myVV[ix][iz] * p1[ix][iz]
                        + myC11 * (myVV[ix][iz + 1] * p1[ix][iz + 1] + myVV[ix][iz - 1] * p1[ix][iz - 1])
                        + myC12 * (myVV[ix][iz + 2] * p1[ix][iz + 2] + myVV[ix][iz - 2] * p1[ix][iz - 2])
                        + myC21 * (myVV[ix + 1][iz] * p1[ix + 1][iz] + myVV[ix - 1][iz] * p1[ix - 1][iz])
                        + myC22 * (myVV[ix + 2][iz] * p1[ix + 2][iz] + myVV[ix - 2][iz] * p1[ix - 2][iz]);
            }
    } else {
#ifdef _OPENMP
#pragma omp parallel for 
#endif
        for (int ix = 2; ix < myNxpad - 2; ix++)
            for (int iz = 2; iz < myNzpad - 2; iz++) {
                float tmp = myC0 * p1[ix][iz] +
                        myC11 * (p1[ix][iz - 1] + p1[ix][iz + 1]) +
                        myC12 * (p1[ix][iz - 2] + p1[ix][iz + 2]) +
                        myC21 * (p1[ix - 1][iz] + p1[ix + 1][iz]) +
                        myC22 * (p1[ix - 2][iz] + p1[ix + 2][iz]);
                p0[ix][iz] = 2 * p1[ix][iz] - p0[ix][iz] + myVV[ix][iz] * tmp;
            }
    }
}

void RtmBorn::stepForwardSO(float** p0, float** p1, bool adj) {
    /*< forward modeling step second order>*/
    float xscale = 1. / (myDx * myDx);
    float zscale = 1. / (myDz * myDz);

    if (adj) { // not implemented correctly yet
#ifdef _OPENMP
#pragma omp parallel for 
#endif
        for (int ix = 1; ix < myNxpad - 1; ix++) {
            for (int iz = 1; iz < myNzpad - 1; iz++) {
                float tmpx = myVV[ix+1][iz]*p1[ix + 1][iz] + myVV[ix - 1][iz]*p1[ix - 1][iz] - 2 * myVV[ix][iz]*p1[ix][iz];
                float tmpz = myVV[ix][iz + 1]*p1[ix][iz + 1] + myVV[ix][iz - 1]*p1[ix][iz - 1] - 2 * myVV[ix][iz]*p1[ix][iz];
                p0[ix][iz] = 2 * p1[ix][iz] - p0[ix][iz]
                        + (xscale * tmpx + zscale * tmpz);
            }
        }

    } else {
#ifdef _OPENMP
#pragma omp parallel for 
#endif
        for (int ix = 1; ix < myNxpad - 1; ix++)
            for (int iz = 1; iz < myNzpad - 1; iz++) {
                float tmpx = p1[ix + 1][iz] + p1[ix - 1][iz] - 2 * p1[ix][iz];
                float tmpz = p1[ix][iz + 1] + p1[ix][iz - 1] - 2 * p1[ix][iz];
                p0[ix][iz] = 2 * p1[ix][iz] - p0[ix][iz]
                        + myVV[ix][iz]*(xscale * tmpx + zscale * tmpz);
            }
    }
}

void RtmBorn::applySponge(float**p0)
/*< apply absorbing boundary condition >*/ {
    int ix, iz, ib, ibx, ibz;
    float w;

#ifdef _OPENMP
#pragma omp parallel for  private(ib,iz,ix,ibz,ibx,w)  
#endif
    for (ib = 0; ib < myNb; ib++) {
        w = myBndr[ib];

        ibz = myNzpad - ib - 1;
        for (ix = 0; ix < myNxpad; ix++) {
            p0[ix][ibz] *= w; /* bottom sponge */
        }

        ibx = myNxpad - ib - 1;
        for (iz = 0; iz < myNzpad; iz++) {
            p0[ib ][iz] *= w; /*   left sponge */
            p0[ibx][iz] *= w; /*  right sponge */
        }
    }

    for (ib = 0; ib < myNbUp; ib++) {
        for (ix = 0; ix < myNxpad; ix++) {
            p0[ix][ib ] *= myBndrUp[ib]; // upper sponge
        }
    }

}

void RtmBorn::addSource(int is, int it, bool add) {
    addSource(&mySxz[is], mySP1, 1, &myWlt[it], add);
}

void RtmBorn::addSource(int *sxz, float **p, int ns, float *source, bool add)
/*< add source term >*/ {
    int is, sx, sz;
    if (add) {
        for (is = 0; is < ns; is++) {
            sx = sxz[is] / myNz + myNb;
            sz = sxz[is] % myNz + myNbUp;
            p[sx][sz] += source[is];
        }
    } else {
        for (is = 0; is < ns; is++) {
            sx = sxz[is] / myNz + myNb;
            sz = sxz[is] % myNz + myNbUp;
            p[sx][sz] -= source[is];
        }
    }
}

void RtmBorn::transpose(float* trans) {
    matrixTranspose(myDObs, trans, myNg, myNt);
}

void RtmBorn::matrixTranspose(float *matrix, float *trans, int n1, int n2)
/*< matrix transpose: matrix tansposed to be trans >*/ {
    int i1, i2;

    for (i2 = 0; i2 < n2; i2++)
        for (i1 = 0; i1 < n1; i1++)
            trans[i2 + n2 * i1] = matrix[i1 + n1 * i2];
}

void RtmBorn::setRecording(sf_file& wave, int ft, int jt) {
    myWavefile = wave;
    myFTrecord = ft;
    myJTrecord = jt;
}

void RtmBorn::saveWave(int it, float** p) {
    if (!((it - myFTrecord) % myJTrecord)) {
        window2d(myV0, p);
        sf_floatwrite(myV0[0], myNx*myNz, myWavefile);
    }

}

bool RtmBorn::isStable(float* velp) {
    float vmax = 0;
    size_t n = myNx*myNz;
    for (size_t i = 0; i < n; i++) vmax = max(vmax, velp[i]);

    bool stable = false;
    float factor = (vmax * myDt * sqrt(1. / (myDx * myDx) + 1. / (myDz * myDz)));
    if (factor < 1) stable = true;

    cerr << " vmax = " << vmax << ", dx = " << myDx << ", dz="
            << myDz << ", dt " << myDt << ", factor= " << factor << endl;

    return stable;

}

void RtmBorn::adjtest() {
    unsigned int dataSize = myNg * myNt*myNs;
    unsigned int modelSize = myNx*myNz;
    float* data1 = new float[dataSize];
    float* data2 = new float[dataSize];
    float* model1 = new float[modelSize];
    float* model2 = new float[modelSize];
    srand48(time(NULL));
    for (unsigned int i = 0; i < dataSize; i++) data1[i] = (drand48() - 0.5);
    for (unsigned int i = 0; i < modelSize; i++) model2[i] = (drand48() - 0.5);
    cerr << "data1 size (random) " << dot(data1,data1,dataSize) << endl;
    cerr << "model2 size (random) " << dot(model2,model2,modelSize) << endl;
    
    memset(model1, 0, modelSize * sizeof (float));
    memset(data2, 0, dataSize * sizeof (float));
    cerr << "calculating adjoint " << endl;
    loop(true, false, model1, data1);
    cerr << "calculating forward " << endl;
    loop(false, false, model2, data2);

    double dotdata2=dot(data2,data2,dataSize);
    double dotmodel1=dot(model1,model1,modelSize);
    double dotdata=dot(data1,data2,dataSize);
    double dotmodel=dot(model1,model2,modelSize);
    cerr << "data2 size (operator) " << dotdata2 << endl;
    cerr << "model1 size (operator) " << dotmodel1 << endl;
    cerr << "data1*data2 size (random*operator) " << dotdata << endl;
    cerr << "model1*model2 size (operator*randon) " << dotmodel << endl;

    fprintf(stderr, "Adjoint test RTM operator: \n dot=%g, dotdata=%g dotmodel=%g\n",
            dotdata / dotmodel, dotdata, dotmodel);

    del1d(data1);
    del1d(data2);
    del1d(model1);
    del1d(model2);
}

void RtmBorn::adjtestLaplac() {
    unsigned int dataSize = myNx*myNz;
    unsigned int modelSize = myNx*myNz;
    float* data1 = new float[dataSize];
    float* data2 = new float[dataSize];
    float* model1 = new float[modelSize];
    float* model2 = new float[modelSize];
    srand48(time(NULL));
    for (unsigned int i = 0; i < dataSize; i++) data1[i] = (drand48() - 0.5);
    for (unsigned int i = 0; i < modelSize; i++) model2[i] = (drand48() - 0.5);
    memset(model1, 0, modelSize * sizeof (float));
    memset(data2, 0, dataSize * sizeof (float));
    applyLaplac(data1, model1, true);
    applyLaplac(data2, model2, false);

    double sum = 0;
    for (unsigned int i = 0; i < dataSize; i++) sum += data1[i] * data2[i];
    double dotdata = sum;

    sum = 0;
    for (unsigned int i = 0; i < modelSize; i++) sum += model1[i] * model2[i];
    double dotmodel = sum;

    fprintf(stderr, "Adjoint test Laplacian \n dot=%g, dotdata=%g dotmodel=%g\n",
            dotdata / dotmodel, dotdata, dotmodel);

    del1d(data1);
    del1d(data2);
    del1d(model1);
    del1d(model2);



}

void RtmBorn::loop(bool adj, bool add, float* model, float* data) {
    //cerr << " starting loop " << adj << endl;
    //if (adj) memset(model,0,myNx*myNz*sizeof(float));  
    if (adj) memset(myModTmp, 0, myNx * myNz * sizeof (float));
    else memset(data, 0, myNt * myNg * myNs * sizeof (float));

    time_t start;
    time_t end;
    double seconds;
    float** ptr;
    struct timeval t1, t2; // millisecond timing

    for (int is = 0; is < myNs; is++) {
        cerr << "shot = " << is + 1 << " of " << myNs << endl;
        time(&start); // time in sec
        int totalTime = 0; // time in msec
        int totalTime2 = 0;
        // initialize is-th source wavefield Ps
        memset(mySP0[0], 0, myNzpad * myNxpad * sizeof (float));
        memset(mySP1[0], 0, myNzpad * myNxpad * sizeof (float));
        memset(myGP0[0], 0, myNzpad * myNxpad * sizeof (float));
        memset(myGP1[0], 0, myNzpad * myNxpad * sizeof (float));
        if (adj) { // migration mm=Lt dd
            gettimeofday(&t1, NULL);
            for (int it = 0; it < myNt; it++) { // in adj shots go forward, g reverse
                addSource(is, it, true);
                stepForward(mySP0, mySP1, false);
                ptr = mySP0;
                mySP0 = mySP1;
                mySP1 = ptr;
                applySponge(mySP0);
                applySponge(mySP1);
                if (myFromBoundary) boundaryRW(mySP0, it, false);
                else rw_snapshot(mySP0, it, false);
            }
            if (myMute) mute(is, data); // test
            for (int it = myNt - 1; it > -1; it--) {
                // reverse time order, Img[]+=Ps[]*Pg[];
                // backpropagate receiver wavefield
                rcvrTimeSlice(myGP1, data, is, it, adj);
                if (myFromBoundary) boundaryRW(mySP0, it, true);
                else rw_snapshot(mySP0, it, true); // read source wavefield.
                // calculate second time derivative 
                sourceTimeDerivative();
                // use second time derivative of source instead of source
                if (myOption < 2) imaging(myModTmp, mySP0, myGP1, adj);
                else imaging(myModTmp, mySPtt, myGP1, adj);
                // reconstruct source wavefield Ps
                ptr = mySP0;
                mySP0 = mySP1;
                mySP1 = ptr;
                stepForward(mySP0, mySP1, false);
                addSource(is, it, false);
                stepForward(myGP0, myGP1, false);
                ptr = myGP0;
                myGP0 = myGP1;
                myGP1 = ptr;
                applySponge(myGP0);
                applySponge(myGP1);
            }
            applyLaplac(myModTmp, model, adj);
            gettimeofday(&t2, NULL);
            totalTime2 += getTime(t1, t2);
            if (0) {
                std::valarray<float> mm(myNz * myNx);
                for (int ix = 0; ix < myNx; ix++) for (int iz = 0; iz < myNz; iz++) mm[ix * myNz + iz] = model[ix * myNz + iz];
                cerr << "adjoint in rtmborn is = " << (mm * mm).sum() << endl;
            }
        } else {
            // in forward both wavefield go in same direction (no need for RWBndr)
            gettimeofday(&t1, NULL);
            applyLaplac(myModTmp, model, adj);
            for (int it = 0; it < myNt; it++) {
                // forward time order Pg += Ps * Img;
                addSource(is, it, true);
                if (myOption == 3) copyWavefield(mySP0, mySPtt); // SPtt=SP0 before overwriting
                stepForward(mySP0, mySP1, false);
                ptr = mySP0;
                mySP0 = mySP1;
                mySP1 = ptr;
                applySponge(mySP0);
                applySponge(mySP1);
                sourceTimeDerivative();
                if (myOption < 2) imaging(myModTmp, mySP0, myGP1, adj);
                else if (myOption == 2) imaging(myModTmp, mySPtt, myGP1, adj);
                else imagingBorn(myModTmp, mySPtt, myGP1, adj);
                rcvrTimeSlice(myGP1, data, is, it, adj);
                stepForward(myGP0, myGP1, true);
                ptr = myGP0;
                myGP0 = myGP1;
                myGP1 = ptr;
                applySponge(myGP1);
                applySponge(myGP0);
            }
            if (myMute) mute(is, data); // test
            gettimeofday(&t2, NULL);
            totalTime += getTime(t1, t2);
        }

        cerr << "timeForward =" << totalTime << endl;
        cerr << "timeBackward =" << totalTime2 << endl;
        time(&end);
        seconds = difftime(end, start);
        sf_warning("shot %d finished: %f\n", is + 1, seconds);

    }
}

void RtmBorn::rcvrTimeSlice(float** g, float* data, int is, int it, bool adj) {
    int gx, gz;
    // addition to limit maximum offset;
    if (myMaxOffset) {
        int sx = mySxz[is] / myNz;
        if (adj) {
            //cerr << "is=" << is << ", sx=" << sx << "sz = " << mySxz[is] % myNz << endl;
            for (int ig = 0; ig < myNg; ig++) {
                gx = myGxz[ig] / myNz;
                gz = myGxz[ig] % myNz;
                //cerr << "ig=" << ig << " gx= " << gx << " gz=" << gz << endl;
                if ((fabs(gx - sx) * myDx) < myMaxOffset)
                    g[gx + myNb][gz + myNbUp] += data[it + ig * myNt + is * myNt * myNg];
            }
        } else {
            for (int ig = 0; ig < myNg; ig++) {
                gx = myGxz[ig] / myNz;
                gz = myGxz[ig] % myNz;
                if ((fabs(gx - sx) * myDx) < myMaxOffset) data[it + ig * myNt + is * myNt * myNg] += g[gx + myNb][gz + myNbUp];
            }
        }

    } else {
        if (adj) {
            for (int ig = 0; ig < myNg; ig++) {
                gx = myGxz[ig] / myNz;
                gz = myGxz[ig] % myNz;
                g[gx + myNb][gz + myNbUp] += data[it + ig * myNt + is * myNt * myNg];
            }
        } else {
            for (int ig = 0; ig < myNg; ig++) {
                gx = myGxz[ig] / myNz;
                gz = myGxz[ig] % myNz;
                data[it + ig * myNt + is * myNt * myNg] += g[gx + myNb][gz + myNbUp];
            }
        }
    }
}

void RtmBorn::imaging(float* m, float** s, float **g, bool adj) {
    if (adj) {
        for (int i2 = 0; i2 < myNx; i2++)
            for (int i1 = 0; i1 < myNz; i1++)
                m[myNz * i2 + i1] +=
                    s[i2 + myNb][i1 + myNbUp] * g[i2 + myNb][i1 + myNbUp];
    } else {
        for (int i2 = 0; i2 < myNx; i2++)
            for (int i1 = 0; i1 < myNz; i1++)
                g[i2 + myNb][i1 + myNbUp] +=
                    s[i2 + myNb][i1 + myNbUp] * m[i1 + myNz * i2];
    }
}

void RtmBorn::imagingBorn(float* m, float** s, float **g, bool adj) {
    if (adj) {
        for (int i2 = 0; i2 < myNx; i2++)
            for (int i1 = 0; i1 < myNz; i1++)
                m[myNz * i2 + i1] +=
                    s[i2 + myNb][i1 + myNbUp] * g[i2 + myNb][i1 + myNbUp];
    } else {
        for (int i2 = 0; i2 < myNx; i2++)
            for (int i1 = 0; i1 < myNz; i1++)
                g[i2 + myNb][i1 + myNbUp] +=
                    (s[i2 + myNb][i1 + myNbUp] * m[i1 + myNz * i2] / myVel[i2][i1] / myVel[i2][i1]);
    }
}

int RtmBorn::getTime(struct timeval& t1, struct timeval& t2) {
    // return time in milliseconds;
    return ((t2.tv_sec - t1.tv_sec)* 1000 + (t2.tv_usec - t1.tv_usec) / 1000);
}

void RtmBorn::copyWavefield(float** p, float** g) {
    int nxlim = myNx + myNb;
    int nzlim = myNz + myNbUp;
    for (int ix = myNb; ix < nxlim; ix++) {
        for (int iz = myNbUp; iz < nzlim; iz++) {
            g[ix][iz] = p[ix][iz];
        }
    }
    return;
}

int RtmBorn::getNx() {
    return myNx;
}

int RtmBorn::getNz() {
    return myNz;
}

int RtmBorn::getNt() {
    return myNt;
}

int RtmBorn::getNg() {
    return myNg;
}

int RtmBorn::getNs() {
    return myNs;
}

float RtmBorn::getDx() {
    return myDx;
}

float RtmBorn::getDz() {
    return myDz;
}

float RtmBorn::getDt() {
    return myDt;
}

float RtmBorn::getDg() {
    return myDg;
}


// mute direct wave

void RtmBorn::mute(int is, float* data) {
    float factor = 1; // amplification factor for mute;
    // first time mute is only direct wave.
    //if (myIter) factor=(myNIter/(myIter+1.0));
    cerr << "factor = " << factor << ", iter= " << myIter << endl;
    int gx, gz, sx, sz;
    for (int ig = 0; ig < myNg; ig++) {
        gx = myGxz[ig] / myNz;
        gz = myGxz[ig] % myNz;
        sx = mySxz[is] / myNz;
        sz = mySxz[is] % myNz;
        float offset = abs(sx - gx) * myDx;
        float depth = abs(sz - gz) * myDz;
        float time = (offset / myDWVel + depth / myDWVel);
        int itime = time / myDt;
        int itfinal = itime + 200 * factor;
        if (itfinal > myNt) itfinal = myNt;
        // if (!(ig%10)) 
        //   cerr << "ig=" << ig << "offset = " << offset << ", depth" << depth 
        // 	   << "time=" << time << "itime=" << itime << endl;
        for (int it = itime; it < itfinal; it++)
            data[it + ig * myNt + is * myNt * myNg] = 0;
    }

}

void RtmBorn::muteoffset(int is, float* data) {
    int sx = mySxz[is] / myNz;
    int gx;
    for (int ig = 0; ig < myNg; ig++) {
        gx = myGxz[ig] / myNz;
        float offset = abs(sx - gx) * myDx;
        if (offset >= myMaxOffset)
            for (int it = 0; it < myNt; it++)
                data[it + ig * myNt + is * myNt * myNg] = 0;

    }
}
// functions to calculate source derivative.

void RtmBorn::sourceTimeDerivative() {
    // calculate source time derivative and store in mySPtt
    if (myOption == 2) applyLaplac(mySPtt, mySP0, false);
    else timeDerivative(mySP0, mySP1);
    return;

}

void RtmBorn::timeDerivative(float** p0, float** p1) {
    float** p_old = mySPtt; // older wavefield stored already;
    int nxlimit = myNx + myNb;
    int nzlimit = myNz + myNbUp;
    for (int i2 = myNb; i2 < nxlimit; i2++)
        for (int i1 = myNbUp; i1 < nzlimit; i1++)
            mySPtt[i2][i1] = (p1[i2][i1] - 2. * p0[i2][i1] + p_old[i2][i1]);
    return;

}

void RtmBorn::applyLaplac(float* r, float* p, int adj) {
    // laplacian wrapper for images (size nx,nz)
    // for options that do not use the laplacian need to copy and exit
    if (!myApplyLaplac) { // bypass laplacian in the image
        if (adj) memcpy(p, r, myNx * myNz * sizeof (float));
        else memcpy(r, p, myNx * myNz * sizeof (float));
        cerr << "skip laplacian" << endl;
        return;
    } else applyLaplac(r, p, adj, myNz, myNx);
    cerr << "laplacian applied to image with adj" << adj << endl;
    return;
}

void RtmBorn::applyLaplac(float** r, float** p, int adj) {
    // laplacian for wavefields (size,nxpad,nzpad);
    applyLaplac(r[0], p[0], false, myNzpad, myNxpad);
    for (int i = 0; i < myNxpad; i++)
        for (int j = 0; j < myNzpad; j++)
            r[i][j] = -r[i][j];

}

void RtmBorn::applyLaplac(float* r, float* p, int adj, int n1, int n2) {
    int n12 = n1*n2;
    if (adj) memset(p, 0, n12 * sizeof (float));
    else memset(r, 0, n12 * sizeof (float));
    for (int i2 = 0; i2 < n2; i2++) {
        for (int i1 = 0; i1 < n1; i1++) {
            int j = i1 + i2*n1;
            if (i1 > 0) {
                if (adj) {
                    p[j - 1] -= r[j];
                    p[j] += r[j];
                } else {
                    r[j] += p[j] - p[j - 1];
                }
            }
            if (i1 < n1 - 1) {
                if (adj) {
                    p[j + 1] -= r[j];
                    p[j] += r[j];
                } else {
                    r[j] += p[j] - p[j + 1];
                }
            }

            if (i2 > 0) {
                if (adj) {
                    p[j - n1] -= r[j];
                    p[j] += r[j];
                } else {
                    r[j] += p[j] - p[j - n1];
                }
            }
            if (i2 < n2 - 1) {
                if (adj) {
                    p[j + n1] -= r[j];
                    p[j] += r[j];
                } else {
                    r[j] += p[j] - p[j + n1];
                }
            }
        }
    }

    return;
}

void RtmBorn::setIter(int iter) {
    myIter = iter;
}

double RtmBorn::dot(float* a, float* b, size_t n){
  double sum=0;
  for (size_t i=0;i< n;i++) sum += a[i]*b[i];
  return sum;
}
