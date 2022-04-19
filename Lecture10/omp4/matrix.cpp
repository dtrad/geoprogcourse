#include <stdio.h> 
#include <omp.h>
#include <iostream>
#include "MyAlloc.hpp"

using namespace std;

void printA(int **A, int n, int m){
  cout << "A = " << endl;
  for (int i=0;i<n;i++){
    for (int j=0;j<m;j++){
      cout << A[i][j] << " ";
    }
    cout << endl;
  }
}

void printV(int *v, int n){
  cout << "vector = " << endl;
  for (int i=0;i<n;i++){
    cout << v[i] << " ";
  }
  cout << endl;
}

bool compareResults(int* a, int *b, int n){
  bool equal=true;
  for (int i=0;i<n;i++)
    if (a[i]!=b[i]){ 
      equal=false;
      //cout << " error " << ((float)(a[i]-b[i]))/a[i] << endl;;
    }
  return equal;
}

void trackThread(int i){
    cout << "thread " << omp_get_thread_num() << " i=" << i << endl;
    cout.flush();
}

int main()
{
  double elapsedTimeSec=-1.;
  double t1, t2;
  int n=10000;
  int m=10000;
  int **A=0;
  new2d(A,n,m);
  int NUMT= omp_get_max_threads();
  cout << "Available " << NUMT << " threads " << endl;
  NUMT-=0; // change 0 for 1 or 2 if your system has enough threads
  cout << "using " << NUMT << " threads " << endl;
  int *x=0;
  new1d(x,m);
  int *b=0;
  new1d(b,n);
  int *c=0;
  new1d(c,n);
 
  srand((unsigned)time(NULL));
  
  cout << "initialize A with random numbers" <<  endl;
  // calculate multiplication in one thread
  t1=omp_get_wtime();
  // Is it rand() function thread safe?
  //#pragma omp parallel for num_threads(NUMT)
  for (int i=0;i<n;i++){
    for (int j=0;j<m;j++){
      A[i][j]= rand()%9 + 1;
    }
  }
  t2=omp_get_wtime();
  cout << "initialization =" << t2 - t1 << " secs" << endl;
  for (int i=0;i<m;i++) x[i]=i;
  //printA(A,n,m);
  //printV(x,m);
  
  cout << " calculate in serial mode" <<  endl;
  // calculate multiplication in one thread
  t1=omp_get_wtime();
  for (int i=0;i<n;i++){
    b[i]=0;
    for (int j=0;j<m;j++)
      b[i]+=(A[i][j]*x[j]);
  }
  t2=omp_get_wtime();
  cout << "serial multiplication =" << t2 - t1 << " secs" << endl;
  // parallel multiplication


  // Now back to threads. This time we initialize threads in the loop itself
  cout << " calculate in parallel mode " << endl;
  t1=omp_get_wtime();
#pragma omp parallel for num_threads(NUMT)
  for (int i=0;i<n;i++){
    //trackThread(i);
    c[i]=0;
    for (int j=0;j<m;j++)
      c[i]+=(A[i][j]*x[j]);
  }
  t2=omp_get_wtime();
  cout << "parallel multiplication =" << t2 - t1 << " secs" << endl;


  if (compareResults(b,c,n))
    cout << "parallel and serial are equal" << endl;
  else
    cout << "ERROR " << endl;



  // A test to see to race conditions  
  cout << " calculate in parallel mode a second time with race condition. " << endl;
  t1=omp_get_wtime();
  for (int i=0;i<n;i++) c[i]=0;
#pragma omp parallel for num_threads(NUMT)
  for (int j=0;j<m;j++){
    //trackThread(j);
    for (int i=0;i<n;i++){
      //#pragma omp critical
      c[i]+=(A[i][j]*x[j]);
    }
  }
  t2=omp_get_wtime();
  cout << "parallel with race condition  =" << t2 - t1 << " secs" << endl;

  if (compareResults(b,c,n))
    cout << "Second test parallel and serial are also equal" << endl;
  else
    cout << "race conditions resulted in ERROR " << endl;
 
  // A final test to see to cache coherence problems  
  // only can do if n==m
  if (n==m){
    // create transpose matrix
    int **B=0;
    new2d(B,m,n);
    for (int i=0;i<n;i++)
      for (int j=0;j<m;j++)
	B[j][i]=A[i][j];
    
    cout << " calculate in parallel mode cache coherence problem. " << endl;
    t1=omp_get_wtime();
#pragma omp parallel for num_threads(NUMT)
    for (int i=0;i<n;i++){
      c[i]=0;
      for (int j=0;j<m;j++) c[i]+=(B[j][i]*x[j]);
    }
    t2=omp_get_wtime();
    cout << "parallel with wrong order  =" << t2 - t1 << " secs" << endl;

    if (compareResults(b,c,n))
      cout << "Third test parallel and serial are also equal" << endl;
    else
      cout << "cache coherence test gives ERROR.. " << endl;
    del2d(B);
  }
  

  del2d(A);
  del1d(x);
  del1d(b);
  del1d(c);
  return 0; 
} 
