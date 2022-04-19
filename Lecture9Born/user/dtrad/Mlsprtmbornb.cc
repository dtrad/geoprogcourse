// 2D prestack least-squares RTM .
//   Copyright (C) 2010 University of Texas at Austin
//   Modified from sflsprtm.  Daniel Trad - 2017-University of Calgary
//  
//   This program is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 2 of the License, or
//   (at your option) any later version.
//  
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//  
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software
//   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

//   Based on sflsprtm program from Penliang Yang (user/pyang/sflsprtm)

//   Daniel Trad, University of Calgary, 2016: 
//   1) Implemented in C++
//   2) Added classes rtmclass, cgsolver
//   3) Added adjoint test, laplacian operator, optional boundary condition, wavefield saving
//   4) Added optional input and output to predict data.
//   5) Rearranged operator to pass adjoint test.
//   6) Added option to save whole wavefield for comparison
//   7) Added options for using second derivative of source wavefield.
//   8) Added mute for direct wave
//   9) Added output for model and residual iterations
//   10) Added image domain inversion

#include <time.h>

#include <rsf.hh>
#include <vector>
#include "rtmborn.hh"
#include "cgsolverborn.hh"


#ifndef MARK
#define MARK fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__);cerr << flush;
#endif

using namespace std;
void getmypar(iRSF& par);
bool getImageIn(std::valarray<float>& mod, int nx, int nz);
void setImageOut(oRSF& img, RtmPar* par);
void setDataOut(oRSF& data, RtmPar* par);

int main(int argc, char** argv)
{   
    sf_init(argc,argv); // Initialize RSF

    iRSF par(0);
    RtmPar* mypar=new RtmPar(par);
    int   ns;par.get("ns",ns,1);
    
    // parameters from command line
    int predict;par.get("predict",predict,0);
    int testadj; par.get("testadj",testadj,0);
    int niter = mypar->niter;
    int imgdomain;par.get("imgdomain",imgdomain,0); // solve in image domain insted of data
    // read velocity
    iRSF vel("vel");    
    mypar->getModelSizes(vel); 
    int nz=mypar->nz;
    int nx=mypar->nx;
    int m = nx*nz;
    std::valarray<float> vels(m);
    vel >> vels;
    
    // read input
    iRSF shots("in");   
    mypar->getDataSizes(shots);
    // test to overwrite dominant frequency from data
    float fmnew=0;
    par.get("fm",fmnew,0);
    if (fmnew > mypar->fm) mypar->fm=fmnew;
    mypar->display();
    int nt=mypar->nt;
    int ng=mypar->ng;
    ns=mypar->ns;
    
    // define here array to be used for model space
    std::valarray<float> mod(nx*nz); // lsrtm
    std::valarray<float> modb(nx*nz); // rtm
    mod = 0; 
    modb = 0;
    int givenImage;par.get("givenimage",givenImage,0);
    bool haveImage=false; // flag to check if image is valid
    // optional input to be use for prediction.
    if (givenImage) haveImage=getImageIn(mod,nx,nz);

    // initialize output 
    oRSF imag("out");
    setImageOut(imag,mypar);
    imag.put("title","LSRTM");

    // initialize secondary migration output
    oRSF imagb("imgrtm");
    setImageOut(imagb,mypar);

    // optional secondary output with predictions
    oRSF datp("datp");
    setDataOut(datp,mypar);


    std::valarray<float> dat(ns*ng*nt);  // input data
    std::valarray<float> dat2(ns*ng*nt); // predictions
    std::valarray<float> residuals(ns*ng*nt);

    shots >> dat;
    
    RtmBorn* rtm = new RtmBorn(mypar,&(vels[0])); 
    rtm->initGeometry(mypar->sxbeg,mypar->szbeg,mypar->jsx,mypar->jsz,ns,"shots");
    rtm->initGeometry(mypar->gxbeg,mypar->gzbeg,mypar->jgx,mypar->jgz,ng,"rcvrs");

    CGSolverBorn* cg= new CGSolverBorn(rtm);
    cg->dot(dat,dat,ns*ng,nt);
    if (mypar->mute) for (int is=0;is<ns;is++) rtm->mute(is,&dat[0]);
    if (mypar->maxoffset) for (int is=0;is<ns;is++) rtm->muteoffset(is,&dat[0]);
    // save a wavefield
    cerr << "creating wave before and after reflection" << endl;
    sf_file wave = sf_output("wave");
    int ft = 0;   // first time to record.
    int jt = 100; // interval between snapshots
    sf_putint(wave,"n1",nz);
    sf_putint(wave,"n2",nx);
    sf_putint(wave,"n3",ns*(nt-ft)/jt);
    sf_putstring(wave,"title","wavefield");
    sf_putstring(wave,"label1","Depth");
    sf_putstring(wave,"unit1","m");
    sf_putstring(wave,"label2","Distance");
    rtm->setRecording(wave,ft,jt);
    
    if (testadj){ 
      rtm->adjtestLaplac();
      rtm->adjtest();
      imag << mod;
      datp << dat2;
      delete rtm;
      exit(0);
    }
    
    cerr << "before CG model has size = " << cg->dot(mod,mod,nx,nz) << endl;
    if (!haveImage) { // if not given calculate
      cerr << "**** Start migration " << endl;
      cerr << " data norm before adjoint is " << (dat*dat).sum() << endl;
      cerr << " model norm before adjoint is " << (modb*modb).sum() << endl;
      rtm->loop(true,false,&modb[0],&dat[0]);
      cerr << " model norm from adjoint is " << (modb*modb).sum() << endl;
      imagb << modb; // RTM (adjoint or simple migration)
      //oRSF* p = &(imagb);delete p;
      if (niter){
        if (!imgdomain){
            cerr << "data domain LSRTM" << endl;
            cg->wpcgnr(dat,mod,residuals,dat2,niter);
        } 
        if (imgdomain){ 
            mod = modb; // start with migration
            cerr << "img domain LSRTM" << endl;
            cg->cgsq(dat,mod,residuals,dat2,niter);
        }
      }
      else{
          cerr << " no iterations, just output adjoint " << endl;
          mod=modb;
      }
    }
    else{ 
      cerr << "use given image, do not calculate " << endl;
    }

    cerr << "After CG model has size = " << cg->dot(mod,mod,nx,nz);

    // write contents of mod to output file.
    imag << mod;   // LSRTM


    if (predict){ // if given output for predictions.
      dat2 = 0;
      rtm->loop(false,false,&mod[0],&dat2[0]);
      //dat2 = dat -dat2; // uncomment this line to output residuals
      datp << dat2;
    }


    sf_fileclose(wave);
    cerr << "migration ends " << endl;
    delete cg;
    cerr << "deleted cg" << endl;
    delete rtm;
    cerr << "deleted rtm" << endl;
    
    exit(0);
}

void getmypar(iRSF& par){
  int ns;par.get("ns",ns,1);
  cerr << "****ns " << ns << endl;
  
}

bool getImageIn(std::valarray<float>& mod, int nx, int nz){
  iRSF imagin("imagin");
  float o1in=0;
  float o2in=0;
  int nxin=0;
  int nyin=1;
  int nzin=0;
  float dzin=1;
  float dxin=1;
  imagin.get("n1",nzin);
  imagin.get("n2",nxin);
  imagin.get("n3",nyin);
  imagin.get("d1",dzin);
  imagin.get("d2",dxin);
  imagin.get("o1",o1in);
  imagin.get("o2",o2in);
  if ((nx==nxin)&&(nz==nzin)){
    cerr << "have valid image, do not re-calculate" << endl;
    imagin >> mod; // use this image to predict;
  }
  return true;
}

void setImageOut(oRSF& img, RtmPar* par){
    img.put("n1",par->nz);
    img.put("n2",par->nx);
    img.put("n3",1);
    img.put("d1",par->dz);
    img.put("d2",par->dx);
    img.put("o1",par->o1);
    img.put("o2",par->o2);
    img.put("label1","Depth");
    img.put("label2","Distance");
    img.put("unit1","m");
    img.put("title","RTM");
}

void setDataOut(oRSF& data, RtmPar* par){
    data.put("n1",par->nt);
    data.put("n2",par->ng);
    data.put("n3",par->ns);
    data.put("d1",par->dt);
    data.put("d2",par->dxg);
    data.put("amp",par->amp);
    data.put("fm",par->fm);
    data.put("ng",par->ng);
    data.put("szbeg",par->szbeg);
    data.put("sxbeg",par->sxbeg);
    data.put("gzbeg",par->gzbeg);
    data.put("gxbeg",par->gxbeg);
    data.put("jsx",par->jsx);
    data.put("jsz",par->jsz);
    data.put("jgx",par->jgx);
    data.put("jgz",par->jgz);
    data.put("label1","Time");
    data.put("unit1","sec");
    data.put("title","Prediction");
}
