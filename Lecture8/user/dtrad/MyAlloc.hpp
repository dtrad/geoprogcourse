#ifndef ALLOC
#define ALLOC

#include <iostream>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

/** Delete a singular object, assign NULL to pointer */
template< class A >
inline void MDEL( A& aa )
{
    if( aa != NULL ) {
        delete aa;
        aa = NULL;
    }
}

/** Delete an array pointer, assign NULL to pointer */
template< class A >
inline void MDELARR( A& aa )
{
    if( aa != NULL ) {
        delete [] aa;
        aa = NULL;
    }
}


/** Template function to allocate 1d buffer */
template<class A>
void new1d( A* &p1d, unsigned int size )
{
  try {
        p1d = new A[size];
  }
  catch( ... ) {
    cout << "bad memory alloc for requested size ";
    cout << sizeof(A)*size << ", line " << __LINE__ << "file " << __FILE__ << endl;
    p1d = NULL;
    throw;  // rethrow to allow for handling of exception by caller(s)
  }
}

/** Template function to delete 1d buffer */
template<class A>
void del1d( A* &p1d )
{
    delete []p1d;
    p1d = 0;
}

/** Template function to increase size 1d buffer */
template<class A>
void resize1d( A* &p1d, size_t size )
{
  delete [] p1d;p1d=0;
  try{
    p1d = new A[size];
  }
  catch( ... ){
    cout << "bad memory alloc for requested size ";
    cout << sizeof(A)*size << ", line " << __LINE__ << "file " << __FILE__ << endl;
    p1d = NULL;
    throw;  // rethrow to allow for handling of exception by caller(s)
  } 
}

/** Template function to allocate 2d buffer */
template<class A>
void new2d( A** &p2d, unsigned int size1, unsigned int size2 )
{
  size_t n = 0;
  try {

#if 0
    //no guarantee proper memory alignment ?!
    unsigned int size_ptr_table = size1*sizeof(A*);
    char *buf = new char[size_ptr_table + size1*size2*sizeof(A)];

#else

    size_t size_ptr_table = (size_t)size1*sizeof(A*) / sizeof(A);

    if( (size_t)size1*sizeof(A*) % sizeof(A) != 0 ) { // such as A is double on 32bit machine
      size_ptr_table += 1;
    }

    n = size_ptr_table + (size_t)size1*size2;
    A *buf = new A[n];
#endif

    p2d = (A**)buf;
    p2d[0] = (A*)(buf + size_ptr_table);
    for( unsigned int i = 1; i < size1; ++i )
      p2d[i] = p2d[i-1] + size2;
  }

  catch( ... ) {
    cout << "bad memory alloc for requested size ";
    cout << sizeof(A)*n << ", line " << __LINE__ << "file " << __FILE__ << endl;
    p2d = NULL;
    throw;  // rethrow to allow for handling of exception by caller(s)
  }
}

/** Template function to delete 2d buffer */
template<class A>
void del2d( A** &p2d )
{
    delete []p2d;
    p2d = 0;
}

/** Template function to allocate 3d buffer */
template<class A>
void new3d( A*** &p3d, unsigned int size1, unsigned int size2, unsigned int size3 )
{
    size_t n = 0;
    try {
#if 0
        unsigned int size_ptr_table1 = size1*sizeof(A**);
        unsigned int size_ptr_table2 = size1*size2*sizeof(A*);
        char *buf = new char[size_ptr_table1 + size_ptr_table2 + size1*size2*size3*sizeof(A)];
#else

        size_t size_ptr_table1 = (size_t)size1*sizeof(A**) / sizeof(A);
        if( size1*sizeof(A**) % sizeof(A) != 0 ) {
            size_ptr_table1 += 1;
        }

        size_t size_ptr_table2 = (size_t)size1*size2*sizeof(A*) / sizeof(A);
        if( size1*size2*sizeof(A*) % sizeof(A) != 0 ) {
            size_ptr_table2 += 1;
        }

        n = size_ptr_table1 + size_ptr_table2 + (size_t)size1*size2*size3;
        A *buf = new A[n];

#endif
        p3d = (A***)buf;
        p3d[0] = (A**)(buf+size_ptr_table1);
        p3d[0][0] = (A*)(buf+size_ptr_table1+size_ptr_table2);
        size_t offset = 0;
        for( unsigned int i = 0; i < size1; ++i ) {
            p3d[i] = p3d[0] + (size_t)size2*i;
            for( unsigned int j = 0; j < size2; ++j, offset += size3 ) {
                p3d[i][j] = p3d[0][0] + offset;
            }
        }
    }
    catch( ... /*STD::bad_alloc ba*/ ) {
      cout << "bad memory alloc for requested size ";
      cout << sizeof(A)*n << ", line " << __LINE__ << "file " << __FILE__ << endl;
      p3d = NULL;
      throw;  // rethrow to allow for handling of exception by caller(s)
    }
}

/** Template function to delete 3d buffer */
template<class A>
void del3d( A*** &p3d )
{
    delete []p3d;
    p3d = 0;
}



#endif

