/* Interface to numerical solvers 
 *   This function solves the system of equations 
     (FH F + I sigmad ) m = FH d 
          
   Daniel Trad, UBC, 2000. 
*/

#include "radonsolver.h"
#include "radonclass.hpp"
#include "su.h"

void radonsolver(radonclass* radon, float** data, char* solver, float sigmad){
    // calculate frequency slices and transpose from x,t to f,x
    float dt = radon->myDt;
    float fmax= radon->myFmax;
    int nt = radon->myNt;
    int nh = radon->myNh;
    int nq = radon->myNq;
    int nfft = 0; // total number of frequencies.
    int nf = 0;   // effective number of frequencies (usually half nfft)
    float df = 0; // frequency interval.
    complex czero=(0,0);
    float ** model=radon->getModel();
    complex ** d2=ealloc2complex(nh,nt);
    complex ** m2=ealloc2complex(nq,nt);

    float* dh=ealloc1float(nh);
    float* Wm=ealloc1float(nq);
    float* g=radon->getMoveoutFunction();
    float* Wd=ealloc1float(nh);

    // set Wd and Wm to 1 for now.
    for (int ih=0;ih<nh;ih++) Wd[ih]=1;
    for (int iq=0;iq<nq;iq++) Wm[iq]=1;
   
    // calculate temporal frequency slices, transpose for efficiency.
    fft_parameters(nt,dt,&nfft,&nf,&df);
    fftgo_xt2fx(-1,data,d2,nh,nt,dt,nfft,nf);
    int maxfreq=(int) (fmax/df);
    if (!fmax) maxfreq=nf;

    fprintf(stderr,"maxfreq=%d, dt=%f, df=%f, nfft=%d nf=%d \n",maxfreq,dt,df,nfft,nf);

    // array to save cost function value at each frequency    
    for (int freq=1;freq<maxfreq;freq++){
        float w=2*PI*freq*df;
        float wa=freqweight(freq,df,fmax-10,fmax);
        radon->makeOperator(w);
        fprintf(stderr,"freq=%f\n",freq*df);
        if (STREQ(solver,"adj_")) radon->applyOperator(d2[freq],m2[freq],true);
        else if (STREQ(solver,"toep")) radon_toeplitz(d2[freq],radon->getOperator(),m2[freq],sigmad,radon->myNh,radon->myNq);
        else if (STREQ(solver,"cgls")) fprintf(stderr,"Your job ...");      
        if (STREQ(solver,"adj_")) for (int iq=0;iq<nq;iq++) m2[freq][iq]/=nh;
        if ((wa<1)&&(wa>0)) for (int iq=0;iq<nq;iq++)  m2[freq][iq]*=wa;
    }
   
    for (int iq=0;iq<nq;iq++) m2[0][iq]=czero;  //dc can not be recovered  
    for (int freq=maxfreq;freq<nf;freq++)  
        for (int iq=0;iq<nq;iq++) m2[freq][iq]=czero;
    fftback_fx2xt(1,model,m2,nq,nt,dt,nfft,nf);  
    free1float(Wd);
    free1float(Wm);
    free1float(dh);
    free2complex(m2);
    free2complex(d2);
    
}

void radoninv(radonclass* radon, float** data){
    // calculate frequency slices and transpose from x,t to f,x
    float dt = radon->myDt;
    float fmax= radon->myFmax;
    int nt = radon->myNt;
    int nh = radon->myNh;
    int nq = radon->myNq;
    int nfft = 0; // total number of frequencies.
    int nf = 0;   // effective number of frequencies (usually half nfft)
    float df = 0; // frequency interval.
    complex czero=(0,0);
    float** model=radon->getModel();
    complex ** d2=ealloc2complex(nh,nt);
    complex ** m2=ealloc2complex(nq,nt);

    float* dh=ealloc1float(nh);
    float* Wm=ealloc1float(nq);
    float* g=radon->getMoveoutFunction();
    float* Wd=ealloc1float(nh);

   
    // calculate temporal frequency slices, transpose for efficiency.
    fft_parameters(nt,dt,&nfft,&nf,&df);
    fftgo_xt2fx(-1,model,m2,nq,nt,dt,nfft,nf);
    int maxfreq=(int) (fmax/df);
    if (!fmax) maxfreq=nf;

    fprintf(stderr,"maxfreq=%d, dt=%f, df=%f, nfft=%d nf=%d \n",maxfreq,dt,df,nfft,nf);
    
    
    for (int freq=1;freq<maxfreq;freq++){
        float w=2*PI*freq*df;
        float wa=freqweight(freq,df,fmax-10,fmax);
        radon->makeOperator(w);
        fprintf(stderr,"freq=%f\n",freq*df);
        radon->applyOperator(d2[freq],m2[freq],false);
        if ((wa<1)&&(wa>0)) for (int ih=0;ih<nh;ih++)  d2[freq][ih]*=wa;
    }
   
    for (int ih=0;ih<nh;ih++) d2[0][ih]=czero;  //dc can not be recovered  
    for (int freq=maxfreq;freq<nf;freq++)  
        for (int ih=0;ih<nh;ih++) d2[freq][ih]=czero;
    fftback_fx2xt(1,data,d2,nh,nt,dt,nfft,nf);  
    
    free1float(Wd);
    free1float(Wm);
    free1float(dh);
    free2complex(m2);
    free2complex(d2);
    
}
float freqweight(int j, float df, float f1, float f2)
/*******************************************************************
return weight for each frequency
******************************************************************
Function parameters:

int j		freq index
float df	freq increment
float f1	taper off freq
float f2	freq beyond which all components are zero
*******************************************************************
Author: John Anderson (visitor to CSM from Mobil) Spring 1993
*******************************************************************/
{
	float w;
	float f=j*df;
        //fprintf(stderr,"f2=%f,f1=%f\n,f=%f\n",f1,f2,f);
	if(f<=f1) return (1.);
	if(f>=f2) return (0.);

	w = (f2-f)/(f2-f1);
	return (w);
}



