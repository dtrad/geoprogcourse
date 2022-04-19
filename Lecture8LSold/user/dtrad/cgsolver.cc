#include "cgsolver.hh"
#include "rtmclass.hh"

#ifndef MARK
#define MARK print("%s @ %u\n",__FILE__,__LINE__)
#endif


CGSolver::CGSolver(RtmClass* rtm){
  myRtm = rtm;
}

CGSolver::~CGSolver(){
  
}

float CGSolver::wpcgnr(fmatrix& data,fmatrix& model,fmatrix& residuals,fmatrix& predictions, int itercg){

  // arguments are data, output model, final residuals, final predictions.
  // predictions can be calculated from residuals directly but need th memory space for calculations

  float J=0;
  int nx=myRtm->getNx();
  int nz=myRtm->getNz();
  int nt=myRtm->getNt();
  int ng=myRtm->getNg();
  int ns=myRtm->getNs();
  int ny=ng*ns;
  float normb=0;
  float eps=1e-7;
  float alpha;
  float alphaden;
  float alphanum;
  float alphanumold;
  float beta;

  // store residual and update norms
  std::valarray<float> rho(itercg+1);
  std::valarray<float> eta(itercg+1);

  // temporary model spaces;
  std::valarray<float> g(nx*nz);
  std::valarray<float> s(nx*nz);
  std::valarray<float> z(nx*nz);

  // first do migration
  model=0;
  residuals = data;
  float normdata=dot(residuals,residuals,ny,nt);
  cerr << "normdata = " << normdata << endl;
  myRtm->loop(true,false,&g[0],&residuals[0]); 
  z=g; // // replace later z=Wm*g;
  s=z; 
  cerr << "normmodel = " << dot(s,s) << endl;
  alphanum=dot(z,g,nx,nz);
  normb=dot(z,z);
  rho[0]=1; // rho is normalized to normb rho[k]=rho[k]/dot(z,z)
  float rhold=eps*2; // force first iteration to always occur
  int k=0;
  cerr << "residuals in internal iteration " << k << " = " << rho[k] << " eta = " << eta[k] << endl;
  while ((rhold > eps)&&(k< itercg)){
    k++;
    myRtm->loop(false,false,&s[0],&predictions[0]);
    alphaden=dot(predictions,predictions,ny,nt);
    alpha=alphanum/alphaden;
    if (1){
      cerr << "alpha = " << alpha << endl;
      cerr << "alphanum (energy of migration)= " << alphanum << endl; 
      cerr << "alphanum (energy of prediction)= " << alphaden << endl; 
      //if no exact adjoint replace step size with some linear search
      float scaleGlobal = lsScale(residuals,predictions); // scale that would LS match amplitudes
      cerr << "scaleGlobal = " << scaleGlobal << endl;
      cerr << "relation alpha/scaleGlobal = " << alpha/scaleGlobal << endl;
      //alpha=scaleGlobal; // change to lsscale
    }
    cerr << "iter " << k << " modelDot before = " << dot(model,model,nx,nz) << endl;
    cerr << "iter " << k << " model s  energy = " << dot(s,s,nx,nz) << endl;
    model=model+alpha*s;
    cerr << "iter " << k << " modelDot after = " << dot(model,model,nx,nz) << endl;
    residuals=residuals-alpha*predictions;
    myRtm->loop(true,false,&g[0],&residuals[0]);
    z=g;
    float tmpdot=dot(z,z,nx,nz);
    rho[k]=tmpdot/normb;
    eta[k]=dot(residuals,residuals,ny,nt)/normdata;
    cerr << "residuals in internal iteration " << k << " model update = " 
	 << rho[k] << " residuals = " << eta[k] << endl;
    alphanumold=alphanum;
    alphanum=dot(z,g,nx,nz);
    beta=alphanum/alphanumold;
    fprintf(stderr,"beta=%e, alphanum=%e, alphanumold=%e \n",beta,alphanum,alphanumold);
    s=z+beta*s;
  }
  int maxk=k+1;
  for (int k=0;k<maxk;k++){
    cerr << " residuals in iteration " << k << " model update = " 
	 << rho[k] << " residuals = " << eta[k] << endl; 
  }
  return J;
}

float CGSolver::dot(fmatrix& data1, fmatrix& data2){
  return (data1*data2).sum();
}

float CGSolver::lsScale(fmatrix& data1, fmatrix& data2){
  return (dot(data1,data2)/dot(data2,data2));
}

double CGSolver::dot(fmatrix& data1, fmatrix& data2, int nx, int nz){
  int index;
  double sum=0;
  for (int ix=0;ix<nx;ix++){
    for (int iz=0;iz<nz;iz++){
      index=ix*nz+iz;
      sum += (data1[index]*data2[index]);
    }
  }

  // test;
  double sum2= (data1*data2).sum();
  cerr << "sum=" << sum << ", sum2 =" << sum2 << endl;
  return sum;
      

}
