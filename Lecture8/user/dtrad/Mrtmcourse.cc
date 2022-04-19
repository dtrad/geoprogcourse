// 2D prestack RTM .
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
#include <time.h>
#include <rsf.hh>
#include <vector>
#include "rtmclass.hh"

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
    std::valarray<float> mod(nx*nz); // rtm
    mod = 0; 
    int givenImage;par.get("givenimage",givenImage,0);
    bool haveImage=false; // flag to check if image is valid
    // optional input to be use for prediction.
    if (givenImage){
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
	haveImage=true;
	imagin >> mod; // use this image to predict;
      }
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
    shots >> dat;

    RtmClass* rtm = new RtmClass(nx,nz,nt,nb,nbup,ng,ns,dx,dz,dt,&(vels[0]),fm,amp); 
    rtm->initGeometry(sxbeg,szbeg,jsx,jsz,ns,"shots");
    rtm->initGeometry(gxbeg,gzbeg,jgx,jgz,ng,"rcvrs");

    // save a wavefield
    cerr << "creating wave before and after reflection" << endl;
    sf_file wave = sf_output("wave");
    int ft = 0;   // first time to record.
    int jt = 100; // interval between snapshots
    sf_putint(wave,"n1",nz);
    sf_putint(wave,"n2",nx);
    sf_putint(wave,"n3",ns*(nt-ft)/jt);
    sf_putstring(wave,"title=","wavefield");
    rtm->setRecording(wave,ft,jt);
    if (testadj){ 
      rtm->adjtest();
      imag << mod;
      datp << dat2;
      delete rtm;
      exit(0);
    }
    if (!haveImage) { // if not given calculate
      cerr << "**** Start migration " << endl;
      rtm->loop(true,false,&mod[0],&dat[0]);
    }
    else{ 
      cerr << "use given image, do not calculate " << endl;
    }

    // write contents of mod to output file.
    imag << mod;   // RTM
    if (predict){ // if given output for predictions.
      cerr << "Calculate prediction from image" << endl;
      dat2 = 0;
      rtm->loop(false,false,&mod[0],&dat2[0]);
      datp << dat2;
    }
    sf_fileclose(wave);
    cerr << "migration ends " << endl;
    delete rtm;
    exit(0);
}
