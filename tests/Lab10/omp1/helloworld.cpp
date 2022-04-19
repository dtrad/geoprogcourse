#include <stdio.h> 
#include <omp.h>
#include <iostream>
using namespace std;
int main()
{
#pragma omp parallel
  {
    printf("thread=%d!\n",omp_get_thread_num());
  }
  cout << "parallel section is over here " << endl;
  return 0; 
} 
