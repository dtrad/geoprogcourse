#include "su.h"
#include "segy.h"
#include "header.h"
#include "geoprog.h"

void save_gather(float **d, int nh, int nt, float dt, char* name);
void smoothing(float **d,int nt,int nx,int nl,int nr,int flag);
void fftgo0(int sign,float **d,complex  **m, int nh, int nt, float dt, int *nf0);
int fftback0(int sign,float **d,complex  **m, int nh, int nt, float dt, int nf0);
complex cdot(complex *x,complex *y,int n);
/*********************** self documentation **********************/
char *sdoc[] = {
  " suadapdiff file1 file2 > output   					",
  "                                      				",
  " subtract file2 from file1 using and scale factor, i.e.  		",
  " file1 - c file2 = output.            				",
  " c is an scale factor calculated to minimize the square error	",
  " Optional parameter                                                  ",  
  " scale=0   if  scale is given then it is not calculated              ",
  " normalize=0 =1 does not subtract but just normalize file2 by        ",
  "                using the least square calculated scale              ",
  " 	   	   (useful to correct distorted amplitudes)             ",
  "					                		",   
  " global=0 (default) one scale for all traces                         ",
  "       =1           one scale for each trace                         ",
  "       =2           one scale for each time slice                    ",
  "       =3           one scale for each sample with smoothing         ",
  "       =4           one scale for each frequency slice               ",
  " verbose=0 (default) =1 print extra info                             ",
  " smooth=0            =1 smooth scales (only global=3)                ",
  " smoothl=11          smooth length for left and right                ",
  " plot=0              =1 plotting for inputs and output               ",
  NULL};
/* Credits:
 *	Daniel Trad. UBC 1999. UoC, 2017
 * Trace header fields accessed: ns, dt, offset, ntr if not null
**************** end self doc ***********************************/

segy tr,tr2; 


int main(int argc, char **argv)
{
  int j,ih;
  register int it;
  int nt;
  int ntr;
  int nh;
  int nh2;
  
  float dt;
  float scale;
  int plot;
  int normalize;
  int global;
  int verbose=0;
  int smooth=0;
  // Initialize 
  struct smooth{
    // smoothing
    int nl;  //  npoints left hand side
    int nr;  //  npoints left hand side
    int flag;  // 1 rectangular, 2 triangular
  } parsmooth;
  int smoothl;
  
  initargs(argc, argv);
  requestdoc(1);
   
  if (!getparint("plot",&plot)) plot=0;
  if (!getparint("normalize",&normalize)) normalize=0;
  if (!getparfloat("scale",&scale)) scale=0;
  if (!getparint("global",&global)) global=4;
  if (!getparint("verbose",&verbose)) verbose=1;
  if (!getparint("smooth",&smooth)) smooth=0;
  if (!getparint("smoothl",&smoothl)) smoothl=5;

  parsmooth.nl=smoothl;
  parsmooth.nr=smoothl;
  parsmooth.flag=2;
  if (verbose) fprintf(stderr,"nl=%d, nr=%d, flag=%d\n",parsmooth.nl,parsmooth.nr,parsmooth.flag);
  
  /* Open two files given as arguments for reading */
  FILE* fp1 = efopen(argv[1], "r");  // file pointer to first input
  FILE* fp2 = efopen(argv[2], "r");  // file pointer to second input
  FILE* tracefp = etmpfile();        // file pointer to temp file for data
  FILE* headerfp = etmpfile();       // file pointer to temp file for headers
  
  if (!fgettr(fp1,&tr)) err("can't read first trace");
  
  dt = 1e-6*(float) tr.dt;
  nt = (int) tr.ns;  
  if (!(ntr=(unsigned long int) tr.ntr)){
      ntr=1;
      while (fgettr(fp1,&tr)) ntr++;
      fprintf(stderr,"found ntr = %d traces in file 1\n",ntr);
      rewind(fp1);
      if (!fgettr(fp1,&tr)) err("can't read first trace after scanning");
  }

      
  float *h=ealloc1float(ntr);
  // store data and header from first input in disk
  j=0;
  do{   // Loop over traces    
    efwrite(&tr,HDRBYTES,1,headerfp);
    efwrite(tr.data,FSIZE, nt, tracefp);  	   
    h[j]=tr.offset;
    j++;
  }while(fgettr(fp1, &tr));
  nh=j;

  fprintf(stderr,"nh=%d,nt=%d\n",nh,nt);

  float **d=ealloc2float(nt,nh);
  float **d2=ealloc2float(nt,nh);
  float **d3=ealloc2float(nt,nh); // diff, used for plots
  if ((!d)||(!d2)||(!d3)) err("can't allocate memory");

  if (!fgettr(fp2,&tr)) err("can't read first trace");
  // read second input data to memory
  ih=0;
  do{   // Loop over traces    
    memcpy((void *) d2[ih],(const void *) tr.data,nt*sizeof(float));
    ih++;
  }while(fgettr(fp2, &tr));
  nh2=ih;

  if (plot){ save_gather(d2,nh,nt,dt,"d2");system("suxwigb  < d2 title=\"d2\" perc=99 & ");}
  

  fprintf(stderr,"nh2=%d\n",nh2);

  if (nh!=nh2) err("nh is different from nh2\n");

  // read header and data from first input
  erewind(tracefp);
  erewind(headerfp);
  for (ih=0;ih<nh;ih++){
    efread(tr.data,FSIZE, nt, tracefp);   
    memcpy((void *) d[ih],(const void *) tr.data,nt*sizeof(float));
  }

  if (plot){ save_gather(d,nh,nt,dt,"d");system("suxwigb  < d title=\"d\" perc=99 & ");}
  
  //depending on user choice perform different dataflow
  if (global == 0){

    fprintf(stderr,"******** calculate one scale for all traces ********* \n");

    float dpd=dot(nh*nt,d[0],d2[0]);
    float dpdp=dot(nh*nt,d2[0],d2[0]);
    if (dpdp==0) err("second input is empty, cannot calculate ");
    /* if scale=0 (default) then calculate the least square scale */
    if (scale==0) scale=dpd/dpdp;
    fprintf(stderr,"scale===>%f\n",scale);

    erewind(tracefp);
    erewind(headerfp);

    for (ih=0;ih<nh;ih++){
      efread(&tr,HDRBYTES,1,headerfp);
      efread(tr.data,FSIZE, nt, tracefp);
      if (normalize==0) // default
	for (it=0;it<nt;it++) tr.data[it]-=scale*d2[ih][it];
      else 
	for (it=0;it<nt;it++) tr.data[it]=scale*d2[ih][it];
    
      puttr(&tr);
    }
  }
  else if (global == 1){
    fprintf(stderr,"******** calculate one scale for each trace ********* \n");
    float dpd=0;
    float dpdp=0;
    float* scalev=ealloc1float(nh);
    
    for (ih=0;ih<nh;ih++){
      
      dpd=dot(nt,d[ih],d2[ih]);
      dpdp=dot(nt,d2[ih],d2[ih]);
      
      /* if scale=0 (default) then calculate the least square scale */
      if (dpdp) scalev[ih]=dpd/dpdp;
      else      scalev[ih]=1;
      if (verbose) fprintf(stderr,"scalev[%d]===>%f\n",ih,scalev[ih]);
    }
    
    erewind(tracefp);
    erewind(headerfp);

    for (ih=0;ih<nh;ih++){
      efread(&tr,HDRBYTES,1,headerfp);
      efread(tr.data,FSIZE, nt, tracefp);
      if (normalize==0) // default
	for (it=0;it<nt;it++) tr.data[it]-=scalev[ih]*d2[ih][it];
      else 
	for (it=0;it<nt;it++) tr.data[it]=scalev[ih]*d2[ih][it];  
      puttr(&tr);
    }
    free1float(scalev);
  }
  else if (global == 2){
    fprintf(stderr,"******** calculate one scale for time slice ********* \n");
    float dpd=0;
    float dpdp=0;
    float* scalev=ealloc1float(nt);
    float *temp1 = ealloc1float(nh);
    float *temp2 = ealloc1float(nh);
    for (it=0;it<nt;it++){
      for (ih=0;ih<nh;ih++){
	temp1[ih]=d[ih][it];
	temp2[ih]=d2[ih][it];
      }      
      dpd=dot(nh,temp1,temp2);
      dpdp=dot(nh,temp2,temp2);
      
      /* if scale=0 (default) then calculate the least square scale */
      if (dpdp) scalev[it]=dpd/dpdp;
      else      scalev[it]=1;
      fprintf(stderr,"scalev[%d]===>%f\n",it,scalev[it]);
    }
    
    erewind(tracefp);
    erewind(headerfp);

    for (ih=0;ih<nh;ih++){
      efread(&tr,HDRBYTES,1,headerfp);
      efread(tr.data,FSIZE, nt, tracefp);
      if (normalize==0) // default
	for (it=0;it<nt;it++) tr.data[it]-=scalev[it]*d2[ih][it];
      else 
	for (it=0;it<nt;it++) tr.data[it]=scalev[it]*d2[ih][it];
      
      puttr(&tr);
    }
    
    free1float(temp1);
    free1float(temp2);
    free1float(scalev);
  }
  else if (global == 3){
    fprintf(stderr,"******** calculate one scale for each sample ********* \n");
    float s1,s2;
    float** scalev=ealloc2float(nt,nh);
    for (it=0;it<nt;it++){
      for (ih=0;ih<nh;ih++){
	s1=d[ih][it];
	s2=d2[ih][it];
	if (s2) scalev[ih][it]= s1*s2/(s2*s2);
        else    scalev[ih][it]=1;
      }
    }

    if (smooth) smoothing(scalev,nt,nh,parsmooth.nl,parsmooth.nr,parsmooth.flag);
    erewind(tracefp);
    erewind(headerfp);
    
    for (ih=0;ih<nh;ih++){
      efread(&tr,HDRBYTES,1,headerfp);
      efread(tr.data,FSIZE, nt, tracefp);
      if (normalize==0) // default
	for (it=0;it<nt;it++){ tr.data[it]-=scalev[ih][it]*d2[ih][it];d3[ih][it]=tr.data[it];}
      else 
	for (it=0;it<nt;it++){ tr.data[it]=scalev[ih][it]*d2[ih][it];d3[ih][it]=tr.data[it];}
      
      puttr(&tr);
    }

    if (plot){ 
      save_gather(d3,nh,nt,dt,"d3");system("suxwigb  < d3 title=\"d3\" perc=99 & ");
      save_gather(scalev,nh,nt,dt,"scalev");system("suximage  < scalev title=\"scalev\" legend=1 & ");
    }
    free2float(scalev);
    
  }
  else if (global == 4){
    if (!tr.dt) err("dt header field must be set");
    if (!tr.ns) err("ns header field must be set");
    if (!tr.ntr) tr.ntr = ntr;
    
    fprintf(stderr,"******** calculate one scale for each frequency slice ********* \n");
    float dpd, dpdp, df;
    
    int nf0,ifreq,maxfreq, nfreq;
    float scale = 0;
    complex** m1 = ealloc2complex(nh,nt); // transpose from (ix,it) to (iw,ix)
    complex** m2 = ealloc2complex(nh,nt);
    float fmax=0;

    fftgo0(-1,d ,m1,nh,nt,dt,&nf0);
    fftgo0(-1,d2,m2,nh,nt,dt,&nf0);
    
    nfreq=nf0/2;
    df=1/(nf0*dt);
    if (fmax==0) maxfreq=nfreq;
    else         maxfreq=(int) (fmax/df);

    for (ifreq=0;ifreq<maxfreq;ifreq++){
      dpd =(cdot(m1[ifreq],m2[ifreq],nh)).r;
      dpdp=(cdot(m2[ifreq],m2[ifreq],nh)).r;
      
      /* if scale=0 (default) then calculate the least square scale */
      if (dpdp) scale=dpd/dpdp;
      else      scale=1;
      if (verbose) fprintf(stderr,"scale[%d]===>%f\n",ifreq,scale);
      for (ih=0;ih<nh;ih++) m2[ifreq][ih]*=scale;
    }

    fftback0(1,d ,m1,nh,nt,dt,nf0);
    fftback0(1,d2,m2,nh,nt,dt,nf0);

    erewind(tracefp);
    erewind(headerfp);
    
    for (ih=0;ih<nh;ih++){
      efread(&tr,HDRBYTES,1,headerfp);
      efread(tr.data,FSIZE, nt, tracefp);
      if (normalize==0){ // default
	for (it=0;it<nt;it++){ 
	  tr.data[it]-=d2[ih][it];
	  d3[ih][it]=tr.data[it];
	}
      }
      else{
	for (it=0;it<nt;it++){ 
	  tr.data[it]=d2[ih][it];
	  d3[ih][it]=tr.data[it];
	}
      }
      
      puttr(&tr);
    }

    if (plot){ 
      save_gather(d3,nh,nt,dt,"d3");
      system("suxwigb  < d3 title=\"d3\" perc=99 & ");
    }

    free2complex(m1);
    free2complex(m2);    
  }  
  free2float(d);
  free2float(d2);
  free2float(d3);
  
  efclose(fp1);
  efclose(fp2);
  efclose(tracefp);
  efclose(headerfp);
  free1float(h);

  return EXIT_SUCCESS;
}




















