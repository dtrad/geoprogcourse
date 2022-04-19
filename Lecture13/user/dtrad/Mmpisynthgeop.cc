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
#include <valarray>
#include <rsf.hh>
#include "fdmod.hh"
// based on Y.Peng programs.
// implemented with class fdmod which uses acoustic FD with different order (2,4 and 8).
// Daniel Trad, November, 2018, GOPH699 Geophysical Programming

using namespace std;

#ifdef USE_MPI
#include <mpi.h>
#endif

int main(int argc, char** argv)
{ 
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int nslaves = world_size -1;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  cerr << "rank =" << rank << endl;
#endif

  float dt;
  int nt;  
  float amp;
  float fm;
  int nb;
  int nbup;
  int ns;
  int ng;
  int jsx;
  int jgx;
  int jsz;
  int jgz;
  int sxbeg;
  int szbeg;
  int gxbeg;
  int gzbeg;
  bool csdgather=false;
  int order;
  int nz,nx;
  
  float dz, dx;
  int master=0;
  std::valarray<float> vels;

  if (rank==0){
    sf_init(argc,argv);
    iRSF par(0);
    par.get("dt",dt,0.001); // time step
    par.get("nt",nt);
    cerr << "nt=" << nt << endl;
    par.get("amp",amp,1.);
    par.get("fm",fm,20);
    par.get("nb",nb,30);
    par.get("nbup",nbup,30);
    par.get("ns",ns,1);
    par.get("ng",ng,100);
    par.get("jsx",jsx,1);
    par.get("jgx",jgx,1);
    par.get("jsz",jsz,0);
    par.get("jgz",jgz,0);
    par.get("sxbeg",sxbeg,1);
    par.get("szbeg",szbeg,1);
    par.get("gxbeg",gxbeg,1);
    par.get("gzbeg",gzbeg,1);
    par.get("csdgather",csdgather,false);
    par.get("order",order,4); // eight order, more precise
    // read input
    iRSF vel("vel");
    vel.get("n1",nz);
    vel.get("n2",nx);
    cerr << "master rank=" << rank << "nx=" << nx << "nz=" << nz << endl;
    vels.resize(nx*nz);
    vel >> vels;
    vel.get("d1",dz);
    vel.get("d2",dx);
  }
  { // if any rank 
    MPI_Bcast(&nt,1,MPI_INT,master,MPI_COMM_WORLD);
    MPI_Bcast(&nb,1,MPI_INT,master,MPI_COMM_WORLD);
    MPI_Bcast(&nbup,1,MPI_INT,master,MPI_COMM_WORLD);
    MPI_Bcast(&ns,1,MPI_INT,master,MPI_COMM_WORLD);
    MPI_Bcast(&ng,1,MPI_INT,master,MPI_COMM_WORLD);
    MPI_Bcast(&jsx,1,MPI_INT,master,MPI_COMM_WORLD);
    MPI_Bcast(&jgx,1,MPI_INT,master,MPI_COMM_WORLD);
    MPI_Bcast(&jsz,1,MPI_INT,master,MPI_COMM_WORLD);
    MPI_Bcast(&jgz,1,MPI_INT,master,MPI_COMM_WORLD);
    MPI_Bcast(&sxbeg,1,MPI_INT,master,MPI_COMM_WORLD);
    MPI_Bcast(&szbeg,1,MPI_INT,master,MPI_COMM_WORLD);
    MPI_Bcast(&gxbeg,1,MPI_INT,master,MPI_COMM_WORLD);
    MPI_Bcast(&gzbeg,1,MPI_INT,master,MPI_COMM_WORLD);
    MPI_Bcast(&order,1,MPI_INT,master,MPI_COMM_WORLD);
    MPI_Bcast(&nx,1,MPI_INT,master,MPI_COMM_WORLD);
    MPI_Bcast(&nz,1,MPI_INT,master,MPI_COMM_WORLD);
    MPI_Bcast(&dt,1,MPI_FLOAT,master,MPI_COMM_WORLD);
    MPI_Bcast(&dx,1,MPI_FLOAT,master,MPI_COMM_WORLD);
    MPI_Bcast(&dz,1,MPI_FLOAT,master,MPI_COMM_WORLD);
    MPI_Bcast(&fm,1,MPI_FLOAT,master,MPI_COMM_WORLD);
    MPI_Bcast(&amp,1,MPI_FLOAT,master,MPI_COMM_WORLD);
    if (rank) vels.resize(nx*nz);
    MPI_Bcast(&vels[0],nx*nz,MPI_FLOAT,master,MPI_COMM_WORLD);
  }
  if (rank == 0){

    // set output
    oRSF shots("shots");
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


    oRSF source("source.rsf");
    source.put("n1",nt);
    source.put("n2",1);
    source.put("d1",dt);
    source.put("label1","time");
    source.put("unit1","s");
    source.put("title","source");

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
    cerr << "rank=" << rank << "nx=" << nx << endl;
    cerr << "rank=" << rank << "nz=" << nz << endl;
    cerr << "rank=" << rank << "nt=" << nt << endl;
    cerr << "rank=" << rank << "nb=" << nb << endl;
    cerr << "rank=" << rank << "nbup=" << nbup << endl;
    cerr << "rank=" << rank << "ng=" << ng << endl;
    cerr << "rank=" << rank << "ns=" << ns << endl;
    cerr << "rank=" << rank << "dx=" << dx << endl;
    cerr << "rank=" << rank << "dz=" << dz << endl;
    cerr << "rank=" << rank << "dt=" << dt << endl;
    cerr << "rank=" << rank << "fm=" << fm << endl;
    cerr << "rank=" << rank << "amp=" << amp << endl;
    cerr << "rank=" << rank << "order=" << order << endl;
    Fdmodb* fdmodb = new Fdmodb(nx,nz,nt,nb,nbup,ng,ns,dx,dz,dt,&(vels[0]),fm,amp,order); 
    std::valarray<float> wavelet(nt);
    float* wav = fdmodb->getWavelet();
    for (int it=0;it<nt;it++) wavelet[it]=wav[it];
    source << wavelet;
    cerr << " source saved " << endl;

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
    delete fdmodb;
    //std::valarray<float> data(ng*nt);

    std::valarray<float> sum(ns*ng*nt);sum=0;
    std::valarray<float> data(ns*ng*nt);
    //reduce data
    //float* data=0;
    MPI_Reduce(&data[0],&sum[0],ns*ng*nt,MPI_FLOAT,MPI_SUM,master,MPI_COMM_WORLD);
    shots << sum;
    
  }
  else{
    cerr << "rank=" << rank << "nx=" << nx << endl;
    cerr << "rank=" << rank << "nz=" << nz << endl;
    cerr << "rank=" << rank << "nt=" << nt << endl;
    cerr << "rank=" << rank << "nb=" << nb << endl;
    cerr << "rank=" << rank << "nbup=" << nbup << endl;
    cerr << "rank=" << rank << "ng=" << ng << endl;
    cerr << "rank=" << rank << "ns=" << ns << endl;
    cerr << "rank=" << rank << "dx=" << dx << endl;
    cerr << "rank=" << rank << "dz=" << dz << endl;
    cerr << "rank=" << rank << "dt=" << dt << endl;
    cerr << "rank=" << rank << "fm=" << fm << endl;
    cerr << "rank=" << rank << "amp=" << amp << endl;
    cerr << "rank=" << rank << "order=" << order << endl;

    Fdmodb* fdmodb = new Fdmodb(nx,nz,nt,nb,nbup,ng,ns,dx,dz,dt,&(vels[0]),fm,amp,order); 
    fdmodb->sg_init(szbeg, sxbeg, jsz, jsx, "shots");
    fdmodb->sg_init(gzbeg, gxbeg, jgz, jgx, "rcvrs");
    float distx=sxbeg-gxbeg;
    //////////////////////////////////////////////////////////
    time_t start;
    time_t end;
    double seconds;
    std::valarray<float> data(ns*ng*nt);
    data=0;
    for (int is=0;is<ns;is++){
      bool cond=((is%nslaves+1)==rank);
      cerr << "rank=" << rank << "shot=" <<  is << " nslaves=" 
	   << nslaves << " condition=" << cond <<endl;
      if ((cond)&&(1)){
	
	time(&start);
	if (csdgather){
	  gxbeg=sxbeg+is*jsx-distx;
	  fdmodb->sg_init(gzbeg,gxbeg,jgz,jgx,"rcvrs");
	}
	cerr << "rank=" << rank << "->shot " << is << endl;
	fdmodb->modelling(is);
	fdmodb->transpose(&data[is*ng*nt]);
	time(&end);
	seconds = difftime(end,start);
	sf_warning("shot %d finished: %f\n", is+1, seconds);
      }
    }
    delete fdmodb;
    // reduce data
    float* sum=0;
    //std::valarray<float> sum(ns*ng*nt);
    MPI_Reduce(&data[0],sum,ns*ng*nt,MPI_FLOAT,MPI_SUM,master,MPI_COMM_WORLD);
  }


#ifdef USE_MPI
  MPI_Finalize();
#endif

  exit(0);
}
