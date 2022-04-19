// 2D prestack least-squares RTM .
//   Copyright (C) 2010 University of Texas at Austin
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
//
//
//   Modified from sflsprtm to pass adjoint test and converted to C++.  
//   Daniel Trad - 2016-University of Calgary. 
//

#include <time.h>
#include <rsf.hh>
#include <vector>
#include "rtmclass.hh"
#include "cgsolver.hh"

#ifndef MARK
#define MARK fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__);cerr << flush;
#endif

using namespace std;

int main(int argc, char** argv)
{   
    sf_init(argc,argv); // Initialize RSF

    iRSF par(0);

    // parameters read from data
    float dt;
    int nt;
    int ng;
    int jsx;
    int jgx;
    int jsz;
    int jgz;
    int sxbeg;
    int szbeg;
    int gxbeg;
    int gzbeg;

    // parameters that may exist in data or in command line
    float amp;par.get("amp",amp,1000.);
    float fm;par.get("fm",fm,10);
    int   ns;par.get("ns",ns,1);

    
    // parameters from command line
    int predict;par.get("predict",predict,0);
    int nb; par.get("nb",nb,30);
    int nbup; par.get("nbup",nbup,30);
    int testadj; par.get("testadj",testadj,0);
    int niter; par.get("niter",niter,1);

    // read velocity
    iRSF vel("vel");
    int nz,nx;
    vel.get("n1",nz);
    vel.get("n2",nx);
    int m = nx*nz;
    std::valarray<float> vels(m);
    vel >> vels;
    float dz, dx;
    vel.get("d1",dz);
    vel.get("d2",dx);
    
    // read input
    iRSF shots("in");
    float dxg;
    shots.get("n1",nt);
    shots.get("n2",ng);
    shots.get("n3",ns);
    shots.get("d1",dt);
    shots.get("d2",dxg);
    shots.get("amp",amp);
    shots.get("fm",fm);
    shots.get("ng",ng);
    shots.get("szbeg",szbeg);
    shots.get("sxbeg",sxbeg);
    shots.get("gzbeg",gzbeg);
    shots.get("gxbeg",gxbeg);
    shots.get("jsx",jsx);
    shots.get("jsz",jsz);
    shots.get("jgx",jgx);
    shots.get("jgz",jgz);
    
    // define here array to be used for model space
    std::valarray<float> mod(nx*nz); // lsrtm
    std::valarray<float> modb(nx*nz); // rtm
    mod = 0; 
    modb = 0;
    int givenImage;par.get("givenimage",givenImage,0);
    bool haveImage=false; // flag to check if image is valid
    // optional input to be use for prediction.
    if (givenImage){
      // if parameter givenImage is true, 
      // read the image from imagin
      // and copy the contents to mod;
      // set the flag haveImage to true
      // and use for predictions
      cerr << "YOUR JOB" << endl;

    }

    // initialize output 
    oRSF imag("out");
    float o1 = 0;
    float o2 = 0;
    imag.put("n1",nz);
    imag.put("n2",nx);
    imag.put("n3",1);
    imag.put("d1",dz);
    imag.put("d2",dx);
    imag.put("o1",o1);
    imag.put("o2",o2);
    imag.put("label1","Depth");
    imag.put("label2","Distance");

    // initialize secondary migration output
    oRSF imagb("imgrtm");
    imagb.put("n1",nz);
    imagb.put("n2",nx);
    imagb.put("n3",1);
    imagb.put("d1",dz);
    imagb.put("d2",dx);
    imagb.put("o1",o1);
    imagb.put("o2",o2);
    imagb.put("label1","Depth");
    imagb.put("label2","Distance");

    // optional secondary output with predictions
    oRSF datp("datp");
    datp.put("n1",nt);
    datp.put("n2",ng);
    datp.put("n3",ns);
    datp.put("d1",dt);
    datp.put("d2",dxg);
    datp.put("amp",amp);
    datp.put("fm",fm);
    datp.put("ng",ng);
    datp.put("szbeg",szbeg);
    datp.put("sxbeg",sxbeg);
    datp.put("gzbeg",gzbeg);
    datp.put("gxbeg",gxbeg);
    datp.put("jsx",jsx);
    datp.put("jsz",jsz);
    datp.put("jgx",jgx);
    datp.put("jgz",jgz);

    std::valarray<float> dat(ns*ng*nt);  // input data
    std::valarray<float> dat2(ns*ng*nt); // predictions
    std::valarray<float> residuals(ns*ng*nt);

    shots >> dat;

    RtmClass* rtm = new RtmClass(nx,nz,nt,nb,nbup,ng,ns,dx,dz,dt,&(vels[0]),fm,amp); 
    rtm->initGeometry(sxbeg,szbeg,jsx,jsz,ns,"shots");
    rtm->initGeometry(gxbeg,gzbeg,jgx,jgz,ng,"rcvrs");

    CGSolver* cg= new CGSolver(rtm);
    cg->dot(dat,dat,ns*ng,nt);
    
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
      rtm->loop(true,false,&modb[0],&dat[0]);
      if (niter) cg->wpcgnr(dat,mod,residuals,dat2,niter);
      else mod=modb;
    }
    else{ 
      cerr << "use given image, do not calculate " << endl;
    }

    cerr << "After CG model has size = " << cg->dot(mod,mod,nx,nz);

    // write contents of mod to output file.
    imag << mod;   // LSRTM
    imagb << modb; // RTM (adjoint or simple migration)

    if (predict){ // if given output for predictions.
      dat2 = 0;
      rtm->loop(false,false,&mod[0],&dat2[0]);
      //dat2 = dat -dat2; // uncomment this line to output residuals
      datp << dat2;
    }

    cerr << "migration ends " << endl;
    delete cg;
    delete rtm;
    
    exit(0);
}
