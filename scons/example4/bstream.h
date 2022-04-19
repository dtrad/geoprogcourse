//==============================================================
// bstream.h -- Declarations for bofstream and bifstream classes
// Time-stamp: <1999-06-29 14:09:57 tswan>
// Copyright (c) 1999 by Tom Swan. All rights reserved.
//==============================================================

#ifndef __bstream_H
#define __bstream_H   // Prevent multiple #includes

#include <iostream>
#include <fstream>
using namespace std;
// Binary output file stream

class bofstream: public ofstream {
public:
  bofstream(const char *fn)
    : ofstream(fn, ios::out | ios::binary) { }
  void writeBytes(const void *, int);
  template <class T>
    bofstream & operator<< (const T & data);
};

template <class T>
bofstream & bofstream::operator<< (const T & data)
{
  writeBytes(&data, sizeof(data));
  return *this;
}

// Binary input file stream

class bifstream: public ifstream {
public:
  bifstream(const char *fn)
    : ifstream(fn, ios::in | ios::binary) { }
  void readBytes(void *, int);
  template <class T>
    bifstream & operator>> (T & data);
};

template <class T>
bifstream & bifstream::operator>> (T & data)
{
  readBytes(&data, sizeof(data));
  return *this;
}

#endif  // __bstream_H





