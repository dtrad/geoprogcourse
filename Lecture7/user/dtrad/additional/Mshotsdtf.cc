// Data modeling by Finite Differences
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
#include "fdmodb.hh"
// based on Y.Peng programs.
// implemented with class fdmodb which uses acoustic FD with different order (2,4 and 8).
// Daniel Trad, March 10, 2016

using namespace std;

int main(int argc, char** argv)
{   
    sf_init(argc,argv); // Initialize RSF

    iRSF par(0);

    float dt;
    par.get("dt",dt,0.001); // time step

    int nt;
    par.get("nt",nt);

    float amp;
    par.get("amp",amp,1000.);
    
    float fm;
    par.get("fm",fm,10);
    
    int nb;
    par.get("nb",nb,30);

    int nbup;
    par.get("nbup",nbup,0);

    int ns;
    par.get("ns",ns,1);

    int ng;
    par.get("ng",ng,100);

    int jsx;
    par.get("jsx",jsx,1);
    
    int jgx;
    par.get("jgx",jgx,1);

    int jsz;
    par.get("jsz",jsz,0);
    
    int jgz;
    par.get("jgz",jgz,0);

    int sxbeg;
    par.get("sxbeg",sxbeg,0);
    
    int szbeg;
    par.get("szbeg",szbeg,0);
    
    int gxbeg;
    par.get("gxbeg",gxbeg,0);

    int gzbeg;
    par.get("gzbeg",gzbeg,0);

    bool csdgather;
    par.get("csdgather",csdgather,false);

    int method;
    par.get("method",method,0); // eight order, more precise
    


    // read input
    iRSF vel("in");
    int nz,nx;
    vel.get("n1",nz);
    vel.get("n2",nx);
    int m = nx*nz;
    std::valarray<float> vels(m);
    vel >> vels;
    float dz, dx;
    vel.get("d1",dz);
    vel.get("d2",dx);
    
    // read density (assume same size)
    std::valarray<float> rho(m);
    if (1){
      iRSF rhofile("rho");
      int nzrho=nz;
      int nxrho=nx;
    
      rhofile.get("n1",nzrho);
      rhofile.get("n2",nxrho);

      if ((nzrho != nz)||(nxrho != nx)){ // how to check for no input?
	cerr << " rho has wrong size " << endl;
	exit(0);
      }
      rhofile >> rho;
    }
    
      
      
    

    // set output
    oRSF shots;
    shots.put("n1",nt);
    shots.put("n2",ng);
    shots.put("n3",ns);
    shots.put("d1",dt);
    shots.put("d2",jgx*dx);
    shots.put("d3",jsx*dx);
    shots.put("o1",0); // need to calculate based on shot depth?
    shots.put("o2",0); 
    shots.put("o3",0); 

    shots.put("o2",gxbeg*dx);
    shots.put("o3",sxbeg*dx);

    shots.put("label1","Time");
    shots.put("label2","Lateral");
    shots.put("label3","Shot");
    shots.put("unit1","sec");
    shots.put("unit2","m");
    shots.put("amp",amp);
    shots.put("fm",fm);
    shots.put("ng",ng);
    shots.put("szbeg",szbeg);
    shots.put("sxbeg",sxbeg);
    shots.put("gzbeg",gzbeg);
    shots.put("gxbeg",gxbeg);
    shots.put("jsx",jsx);
    shots.put("jsz",jsz);
    shots.put("jgx",jgx);
    shots.put("jgz",jgz);
    shots.put("csdgather",csdgather?1:0);
    shots.put("title","shots");


    Fdmodb* fdmodb = new Fdmodb(nx,nz,nt,nb,nbup,ng,ns,dx,dz,dt,&(vels[0]),&(rho[0]),fm,amp, method); 


    // example outputing an auxiliary file
    oRSF source("source.rsf");
    source.put("n1",nt);
    source.put("n2",1);
    source.put("d1",dt);
    source.put("label1","time");
    source.put("unit1","s");
    source.put("title","source");
    
    
    std::valarray<float> wavelet(nt);
    float* wav = fdmodb->getWavelet();
    for (int it=0;it<nt;it++) wavelet[it]=wav[it];
    source << wavelet;
   
    cerr << " source saved " << endl;

    // oRSF wavefield("wavp.rsf");
    // wavefield.put("n1",nz);
    // wavefield.put("n2",nx);
    

    // hard coded parameters, change later
    int ft = 0; // first time to record
    int jt = 100; // interval between snapshot
    sf_file wave = sf_output("wave.rsf");
    sf_putint(wave,"n1",nz);
    sf_putint(wave,"n2",nx);
    sf_putint(wave,"n3",ns*(nt-ft)/jt);
    fdmodb->setRecording(wave,ft,jt);
    
    //fdmodb->saveWave(wavefield);

    // should this be encapsulated in a class or method?
    // set source geometry ////////////////////////////////////
    if (!(sxbeg>=0 && szbeg>=0 && sxbeg+(ns-1)*jsx<nx && szbeg+(ns-1)*jsz<nz)){
      cerr << "last shot at = " << sxbeg + (ns-1) * jsx << " last cell is " << nx << endl;
      cerr << "max number of shots with this geometry is " << (nx - sxbeg)/jsx + 1 << endl; 
      sf_error("sources exceeds the computing zone!\n"); 
      exit(1);
    }
    fdmodb->sg_init(szbeg, sxbeg, jsz, jsx, "shots");

    // set receiver geometry
    if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx && gzbeg+(ng-1)*jgz<nz)){	
      sf_error("geophones exceeds the computing zone!\n"); 
      exit(1);
    }
    float distx=sxbeg-gxbeg;
    float distz=szbeg-gzbeg;
    // set offsets 
    if (csdgather &&
	!((sxbeg+(ns-1)*jsx)+(ng-1)*jgx-distx <nx  &&
	  (szbeg+(ns-1)*jsz)+(ng-1)*jgz-distz <nz)){
      sf_error("geophones exceeds the computing zone!\n"); 
      exit(1);
    }
    fdmodb->sg_init(gzbeg, gxbeg, jgz, jgx, "rcvrs");

    //////////////////////////////////////////////////////////
    time_t start;
    time_t end;
    double seconds;
    std::valarray<float> data(ng*nt);
    //for (int ix=0;ix< (ng*nt);ix++) data[ix]=1;
    for (int is=0;is<ns;is++){
      time(&start);
      if (csdgather){
	gxbeg=sxbeg+is*jsx-distx;
	fdmodb->sg_init(gzbeg,gxbeg,jgz,jgx,"rcvrs");
      }
      cerr << "shot " << is << endl;
      fdmodb->modelling(is);
      fdmodb->transpose(&data[0]);
      shots << data;
      
      time(&end);
      seconds = difftime(end,start);
      sf_warning("shot %d finished: %f\n", is+1, seconds);
    }
    delete fdmodb;
    exit(0);
}
