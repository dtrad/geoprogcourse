//==============================================================
// wdouble.cpp -- Writes double values in binary to disk
// Time-stamp: <1999-06-29 11:48:04 tswan>
// To compile:
//   g++ -o wdouble wdouble.cpp
// To run:
//   ./wdouble
// Copyright (c) 1999 by Tom Swan. All rights reserved.
//==============================================================

#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>

#define FILENAME "test.dat"

class bofstream: public ofstream {
public:
  bofstream(const char *fn)
    : ofstream(fn, ios::out | ios::binary) { }
  void writeBytes(const void *, int);
  bofstream & operator<< (double);
};

bofstream & bofstream::operator<< (double d)
{
  writeBytes(&d, sizeof(d));
  return *this;
}

int main()
{
  bofstream bofs(FILENAME);
  if (!bofs) {
    cerr << "Error: unable to write to " << FILENAME << endl;
    exit(1);
  }
  cout << "Writing to " << FILENAME << endl;
  double d = 3.14159;
  bofs << d;
  bofs << d * d;
  bofs << 9.9999999;
  d = 4.7E-8;
  bofs << d;
  return 0;
}

void bofstream::writeBytes(const void *p, int len)
{
  if (!p) return;
  if (len <= 0) return;
  write((char *)p, len);
}
