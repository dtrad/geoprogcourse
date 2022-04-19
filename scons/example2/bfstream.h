//==============================================================
// bfstream.h -- Declarations for bfstream (binary file) class
// Time-stamp: <1999-06-29 14:34:55 tswan>
// Copyright (c) 1999 by Tom Swan. All rights reserved.
//==============================================================

#ifndef __bfstream_H
#define __bfstream_H   // Prevent multiple #includes

#include <iostream>
#include <fstream>

// Binary input and output file stream
using namespace std;
class bfstream: public fstream {
public:
  bfstream(const char *fn)
    : fstream(fn, ios::in | ios::out | ios::binary) { }
  void writeBytes(const void *, int);
  void readBytes(void *, int);
  template <class T>
    bfstream & operator<< (const T & data);
  template <class T>
    bfstream & operator>> (T & data);    
};

template <class T>
bfstream & bfstream::operator<< (const T & data)
{
  writeBytes(&data, sizeof(data));
  return *this;
}

template <class T>
bfstream & bfstream::operator>> (T & data)
{
  readBytes(&data, sizeof(data));
  return *this;
}

#endif  // __bfstream_H
