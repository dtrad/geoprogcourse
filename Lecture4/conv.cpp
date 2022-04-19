void conv(int lx, int ifx, float *x, int ly, int ify, float *y,
	  int lz, int ifz, float *z){
  int ilx=ifx+lx-1,ily=ify+ly-1,ilz=ifz+lz-1,i,j,jlow,jhigh;
  float sum;
  
  x -= ifx;  y -= ify;  z -= ifz;
  for (i=ifz; i<=ilz; ++i) {
    jlow = i-ily;  if (jlow<ifx) jlow = ifx;
    jhigh = i-ify;  if (jhigh>ilx) jhigh = ilx;
    for (j=jlow,sum=0.0; j<=jhigh; ++j) sum += x[j]*y[i-j];
    z[i] = sum;
  }
}
