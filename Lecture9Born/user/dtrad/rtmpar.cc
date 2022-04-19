#include "rtmpar.hh"
RtmPar::RtmPar(iRSF& par){
  par.get("ns",ns,1);
  par.get("option",option,2);  // default take second time derivative by laplacian
  par.get("laplac",applyLaplac,1); // default apply Laplacian to image
  par.get("mute",mute,1); // apply mute to the first arrival
  par.get("nb",nb,30);
  par.get("nbup",nbup,30);
  par.get("niter",niter,1);
  par.get("maxoffset",maxoffset,0);
  par.get("order",order,4);
  
}
void RtmPar::getModelSizes(iRSF& vel){
  vel.get("n1",nz);
  vel.get("n2",nx);
  vel.get("d1",dz);
  vel.get("d2",dx);
  vel.get("o1",o1);
  vel.get("o2",o2);
}

void RtmPar::getDataSizes(iRSF& shots){
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
}

void RtmPar::display(){
  cerr << "input parameters " << endl;
  cerr << "nx=" << nx << endl;
  cerr << "nz=" << nz << endl;
  cerr << "ng=" << ng << endl;
  cerr << "nt=" << nt << endl;
  cerr << "ns=" << ns << endl;
  cerr << "nb=" << nb << endl;
  cerr << "nbup=" << nbup << endl;
  cerr << "dt=" << dt << endl;
  cerr << "dx=" << dx << endl;
  cerr << "dz=" << dz << endl;
  cerr << "dxg=" << dxg << endl;
  cerr << "amp=" << amp << endl;
  cerr << "fm=" << fm << endl;
  cerr << "o1=" << o1 << endl;
  cerr << "o2=" << o2 << endl;
  cerr << "applyLaplac=" << applyLaplac << endl;
  cerr << "mute=" << mute << endl;
  cerr << "option=" << option << endl;
  cerr << "order="  << order << endl;

}
