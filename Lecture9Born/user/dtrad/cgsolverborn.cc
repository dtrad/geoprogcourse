#include "cgsolverborn.hh"
#include "rtmborn.hh"

#ifndef MARK
#define MARK print("%s @ %u\n",__FILE__,__LINE__)
#endif


CGSolverBorn::CGSolverBorn(RtmBorn* rtm){
  myRtm = rtm;
  myModelFile=0;         // used for both 
  myResidualFile=0;      // used for DATA DOMAIN residual
  myResidualModelFile=0; // used for IMAGE DOMAIN residual
  
}

CGSolverBorn::~CGSolverBorn(){
  cerr << "closing files " << endl;
  if (myModelFile){    sf_fileclose(myModelFile);cerr << " close modelfile" << endl;}
  if (myResidualFile){ sf_fileclose(myResidualFile);cerr << "close residualfile" << endl;}
  if (myResidualModelFile){ sf_fileclose(myResidualModelFile);cerr << "close residual model file " << endl;}
}

float CGSolverBorn::wpcgnr(fmatrix& data,fmatrix& model,fmatrix& residuals,fmatrix& predictions, int itercg){

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
  std::valarray<double> alphav(itercg+1);

  // temporary model spaces;
  std::valarray<float> g(nx*nz);
  std::valarray<float> s(nx*nz);
  std::valarray<float> z(nx*nz);

  setModelFile(itercg);
  setResidualFile(itercg+1);

  // first do migration
  model=0;
  residuals = data;
  float normdata=dot(residuals,residuals,ny,nt);
  cerr << "normdata = " << normdata << endl;
  myRtm->loop(true,false,&g[0],&residuals[0]); 
  // this first model will be calculated after scaling by alpha

  z=g; // // replace later z=Wm*g;
  s=z; 
  cerr << "normmodel = " << dot(s,s) << endl;
  alphanum=dot(z,g,nx,nz);
  normb=dot(z,z);
  rho[0]=1; // rho is normalized to normb rho[k]=dot(model,model)/normb;
  eta[0]=1; // eta is normalized to normdata (eta[k]=dot(residuals,residuals)/normdata;
  alphav[0]=1;
  float rhold=eps*2; // force first iteration to always occur
  int k=0;
  cerr << "residuals in internal iteration " << k << " = " << rho[k] << " eta = " << eta[k] << endl;
  while ((rhold > eps)&&(k< itercg)){
    k++;
    myRtm->setIter(k);
    myRtm->loop(false,false,&s[0],&predictions[0]);
    alphaden=dot(predictions,predictions,ny,nt);
    alpha=alphanum/alphaden;
    alphav[k]=alpha;
    if (1){
      cerr << "alpha = " << alpha << endl;
      cerr << "alphanum (energy of migration)= " << alphanum << endl; 
      cerr << "alphaden (energy of prediction)= " << alphaden << endl; 
      //if no exact adjoint replace step size with some linear search
      float scaleGlobal = lsScale(residuals,predictions); // scale that would LS match amplitudes
      cerr << "scaleGlobal = " << scaleGlobal << endl;
      cerr << "relation alpha/scaleGlobal = " << alpha/scaleGlobal << endl;
      //alpha=scaleGlobal; // change to lsscale
    }
    cerr << "iter " << k << " modelDot before = " << dot(model,model,nx,nz) << endl;
    cerr << "iter " << k << " model s  energy = " << dot(s,s,nx,nz) << endl;
    model=model+alpha*s;
    sf_floatwrite(&model[0],nx*nz,myModelFile);
    cerr << "iter " << k << " modelDot after = " << dot(model,model,nx,nz) << endl;
    if (k==1) sf_floatwrite(&residuals[0],nt*ng*ns,myResidualFile);
    residuals=residuals-alpha*predictions;
    sf_floatwrite(&residuals[0],nt*ng*ns,myResidualFile);
    myRtm->loop(true,false,&g[0],&residuals[0]);
    z=g;
    float tmpdot=dot(z,z,nx,nz);
    rho[k]=tmpdot/normb;
    eta[k]=dot(residuals,residuals,ny,nt)/normdata;
    cerr << "residuals in internal iteration " << k << " model update = " 
	 << rho[k] << " residuals = " << eta[k] << endl;
    cerr << "energy of residuals = " << dot(residuals,residuals,ny,nt) << endl;
    cerr << "energy of normdata  = " << normdata << endl;
    cerr << "energy of updates = " << tmpdot << endl;
    cerr << "energy of migration = " << normb << endl;
    alphanumold=alphanum;
    alphanum=dot(z,g,nx,nz);
    beta=alphanum/alphanumold;
    fprintf(stderr,"beta=%e, alphanum=%e, alphanumold=%e \n",beta,alphanum,alphanumold);
    s=z+beta*s;
  }
  int maxk=k+1;
  for (int k=0;k<maxk;k++){
    cerr << " normalized residuals in iteration " << k << " model update = " 
	 << rho[k] << " residuals = " << eta[k] << endl; 
  }
  for (int k=0;k<maxk;k++){
    cerr << " residuals in iteration " << k << " model update = " 
	 << rho[k]*normb << " residuals = " << eta[k]*normdata << endl; 
  }
  for (int k=0;k<maxk;k++){
    cerr << " alpha in iteration " << k << " = " << alphav[k] << endl; 
  }

  return J;
}
float CGSolverBorn::cgsq(fmatrix& data,fmatrix& model,fmatrix& residuals,fmatrix& predictions, int itercg){

  // argument is migrated data (stored in output model), 
  // data, residuals and predictions are not really used for now.
  
  float J=0;
  int nx=myRtm->getNx();
  int nz=myRtm->getNz();
  float alpha; // need float for template valarray
  double alphaden;
  double eps=1e-5;
  float  eps2=1e-10;
  // store residual and update norms
  std::valarray<float> rho(itercg+1);
  // temporary model spaces;
  std::valarray<float> Ap(nx*nz);
  std::valarray<float> p(nx*nz); // model space prediction
  std::valarray<float> r(nx*nz); //model space residual
  
  setModelFile(itercg);
  setResidualModelFile(itercg+1);
  // save first model version (migration)
  sf_floatwrite(&model[0],nx*nz,myModelFile);
  // migration already done and store in model.
  myRtm->loop(false,false,&model[0],&predictions[0]);
  predictions*=eps2;
  myRtm->loop(true,false,&Ap[0],&predictions[0]);
  r = model-Ap;
  p = r;
  sf_floatwrite(&r[0],nx*nz,myResidualModelFile);
  double rsold0=dot(r,r,nx,nz);
  double rsold=rsold0;
  cerr << "normdata (image space)= " << rsold << endl;
  rho[0]=1; // rho is normalized to normb rho[k]=rho[k]/dot(z,z)
  int k=0;
  cerr << "residuals in internal iteration " << k << " = " << rho[k]  << endl;
  while ((rsold > eps)&&(k< itercg)){
    k++;
    // operator is (modeling + migration)
    myRtm->loop(false,false,&p[0],&predictions[0]); // modeling
    predictions*=eps2;
    myRtm->loop(true,false,&Ap[0],&predictions[0]); // migration
    
    alphaden=dot(p,Ap,nx,nz);
    alpha=rsold/alphaden;
    cerr << "alphaden (p,Ap) " << alphaden << endl;
    cerr << "alpha " << alpha << endl;
    model=model+alpha*p;
    r=r-alpha*Ap;
    sf_floatwrite(&r[0],nx*nz,myResidualModelFile);
    float rsnew=dot(r,r,nx,nz);
    sf_floatwrite(&model[0],nx*nz,myModelFile);
    cerr << "iter " << k << " modelDot after = " << dot(model,model,nx,nz) << endl;
    if (sqrt(rsnew) < eps){ break;cerr << "******BREAK**********" << endl;}
    rho[k]=rsnew/rsold0;
    float beta=(float) (rsnew/rsold);
    p = r + beta * p;
    rsold = rsnew;
  }
  int maxk=k+1;
  for (int k=0;k<maxk;k++){
    cerr << " residuals in iteration " << k << " model update = " << rho[k] << endl; 
  }

  return J;
}

float CGSolverBorn::dot(fmatrix& data1, fmatrix& data2){
  return (data1*data2).sum();
}

float CGSolverBorn::lsScale(fmatrix& data1, fmatrix& data2){
  return (dot(data1,data2)/dot(data2,data2));
}

double CGSolverBorn::dot(fmatrix& data1, fmatrix& data2, int nx, int nz){
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

void CGSolverBorn::setModelFile(int itercg){
  myModelFile = sf_output("modeliter.rsf");
  sf_putint(myModelFile,"n1",myRtm->getNz());
  sf_putint(myModelFile,"n2",myRtm->getNx());
  sf_putint(myModelFile,"n3",itercg);
  sf_putfloat(myModelFile,"d1",myRtm->getDz());
  sf_putfloat(myModelFile,"d2",myRtm->getDx());
  sf_putfloat(myModelFile,"d3",1);

  sf_putstring(myModelFile,"title","LSRTM");
  sf_putstring(myModelFile,"label3","iterations");
  sf_putstring(myModelFile,"label2","Distance");
  sf_putstring(myModelFile,"label1","Depth");
  sf_putstring(myModelFile,"unit1","m");
}


void CGSolverBorn::setResidualFile(int itercg){
  myResidualFile = sf_output("resiter.rsf");
  sf_putint(myResidualFile,"n1",myRtm->getNt());
  sf_putint(myResidualFile,"n2",myRtm->getNg()*myRtm->getNs());
  sf_putint(myResidualFile,"n3",itercg);
  sf_putfloat(myResidualFile,"d1",myRtm->getDt());
  sf_putfloat(myResidualFile,"d2",myRtm->getDg());
  sf_putfloat(myResidualFile,"d3",1);

  sf_putstring(myResidualFile,"title","Residuals");
  sf_putstring(myResidualFile,"label3","iterations");
  sf_putstring(myResidualFile,"label2","Distance");
  sf_putstring(myResidualFile,"label1","Time");
  sf_putstring(myResidualFile,"unit1","s");
}

void CGSolverBorn::setResidualModelFile(int itercg){
  myResidualModelFile = sf_output("resmodeliter.rsf");
  sf_putint(myResidualModelFile,"n1",myRtm->getNz());
  sf_putint(myResidualModelFile,"n2",myRtm->getNx());
  sf_putint(myResidualModelFile,"n3",itercg);
  sf_putfloat(myResidualModelFile,"d1",myRtm->getDz());
  sf_putfloat(myResidualModelFile,"d2",myRtm->getDx());
  sf_putfloat(myResidualModelFile,"d3",1);

  sf_putstring(myResidualModelFile,"title","Residuals Image domain");
  sf_putstring(myResidualModelFile,"label3","iterations");
  sf_putstring(myResidualModelFile,"label2","Distance");
  sf_putstring(myResidualModelFile,"label1","Depth");
  sf_putstring(myResidualModelFile,"unit1","m");
}
