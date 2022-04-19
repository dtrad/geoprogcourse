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
    if (a[i]!=b[i]){ equal=false;break;}
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
  int n=5;
  int m=10;
  int **A=0;
  new2d(A,n,m);
  cout << "Using " << omp_get_max_threads() << " threads " << endl;
  int *x=0;
  new1d(x,m);
  int *b=0;
  new1d(b,n);
  int *c=0;
  new1d(c,n);


  // First an example where the threads are initialized once and used several times.
  omp_set_nested(2);
#pragma omp parallel // all threads initialized here once
  {
    printf("thread=%d!\n",omp_get_thread_num());
#pragma omp barrier // wait to all is done. 
#pragma omp single  // need to do something once!
    {
      trackThread(999);
      srand((unsigned)time(NULL));
    }
#pragma omp for     // back to use all threads
      for (int i=0;i<n;i++){
	trackThread(i);
	for (int j=0;j<m;j++){
	  A[i][j]= rand()%9 + 1;
	}
      }
  } // closing threads we initialized earlier


  // lets us do some serial stuff.
  for (int i=0;i<m;i++) x[i]=i;
  printA(A,n,m);
  printV(x,m);
  
  cout << " calculate in serial mode" <<  endl;
  // calculate multiplication in one thread
#pragma omp single
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
#pragma omp for
  for (int i=0;i<n;i++){
    //cout << "thread " << omp_get_thread_num() << " for row " << i << endl;
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



  // A final test to see how much it costs to initialize threads. 
  cout << " calculate in parallel mode a second time. " << endl;
  {
#pragma omp single    
    t1=omp_get_wtime();
#pragma omp for
  for (int i=0;i<n;i++){
    //cout << "thread " << omp_get_thread_num() << " for row " << i << endl;
    c[i]=0;
    for (int j=0;j<m;j++)
      c[i]+=(A[i][j]*x[j]);
  }

#pragma omp single
  t2=omp_get_wtime();
  }
  
  cout << "parallel multiplication =" << t2 - t1 << " secs" << endl;


  if (compareResults(b,c,n))
    cout << "Second test parallel and serial are also equal" << endl;
  else
    cout << "ERROR " << endl;

  del2d(A);
  del1d(x);
  del1d(b);
  del1d(c);
  return 0; 
} 
