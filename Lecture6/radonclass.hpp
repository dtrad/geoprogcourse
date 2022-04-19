/* 
 * File:   radonclass.hpp
 * Author: dtrad
 *
 * Created on October 7, 2017, 5:21 PM
 */

#ifndef RADONCLASS_HPP
#define	RADONCLASS_HPP
#include "Complex.h"

class radonclass {
public:
    radonclass(float dt, int nt, float moveoutmin, float moveoutmax, float maxoffset, float aperture, int rtmethod, float fmax);
    void initGroup(float* h, int nh);
    void makeOperator(float w);
    void applyOperator(complex* data, complex* model, int adj);
    void adjustAxis(float* h, int nh, float dx);
    void moveout(float* h, float*g, int nh);
    float**    getModel();
    complex** getOperator();
    float*    getMoveoutFunction();
    ~radonclass();

    float myFmax;    
    float myDt;   
    int myNt;
    int myNh;
    int myNq;    
    
private:
    float *myG;
    float* myQ;
    float myQmin;
    float myQmax;
    float myFactor;
    float myDq;
    complex** myOper;
    float**   myModel;
    int myRtmethod;

};

#endif	/* RADONCLASS_HPP */

