#include <stdio.h> 
#include <omp.h>
#include <iostream>

using namespace std;

int slowFunction(int p){
  // use here any slow function returns 1 or 0
  // example check a number is prime or not.
  int isPrime=1;
  for (int i=2;i<p;i++){
    if ( p%i == 0){ isPrime = 0; break;}
  }
  return isPrime;

}



void trackThread(int i){
    cout << "thread " << omp_get_thread_num() << " i=" << i << endl;
    cout.flush();
}

int main()
{
  double elapsedTimeSec=-1.;
  double t1, t2;
  int N=10000;
  int numPrimes = 0;
  
  // Serial run
  t1=omp_get_wtime();
  numPrimes = 0;
  for (int i=1;i<N;i++)  numPrimes += slowFunction(i);
  t2=omp_get_wtime();
  cout << " serial run gives " << numPrimes 
       << " in " << t2-t1 << " secs " << endl;


  // Atomic run
  t1=omp_get_wtime();
  numPrimes = 0;
#pragma omp parallel for 
  for (int i=1;i<N;i++){
    //#pragma omp atomic
    numPrimes += slowFunction(i); 
  }
  t2=omp_get_wtime();
  cout << " Atomic run gives " << numPrimes 
       << " in " << t2-t1 << " secs " << endl;  

  // Critical run 
  t1=omp_get_wtime();
  numPrimes = 0;
#pragma omp parallel for 
  for (int i=1;i<N;i++){
    int temp = slowFunction(i);
    //#pragma omp critical
    numPrimes += temp;
  }
  t2=omp_get_wtime();
  cout << " Critical run gives " << numPrimes 
       << " in " << t2-t1 << " secs " << endl;

  // reduction run 
  t1=omp_get_wtime();
  numPrimes = 0;
#pragma omp parallel for reduction(+:numPrimes)
  for (int i=1;i<N;i++){
    numPrimes += slowFunction(i);
  }
  t2=omp_get_wtime();
  cout << " Reduction run gives " << numPrimes 
       << " in " << t2-t1 << " secs " << endl;



  return 0; 
} 
