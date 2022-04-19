//==============================================================
// rdouble.cpp -- Read double values in binary from file
// Time-stamp: <1999-06-29 11:55:41 tswan>
// To compile:
//   g++ -o wdouble wdouble.cpp
//   g++ -o rdouble rdouble.cpp
// To run:
//   ./wdouble
//   ./rdouble
// Copyright (c) 1999 by Tom Swan. All rights reserved.
//==============================================================

#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>

#define FILENAME "test.dat"

class bifstream: public ifstream {
public:
  bifstream(const char *fn)
    : ifstream(fn, ios::in | ios::binary) { }
  void readBytes(void *, int);
  bifstream & operator>> (double &);
};

bifstream & bifstream::operator>> (double &d)
{
  readBytes(&d, sizeof(d));
  return *this;
}

int main()
{
  bifstream bifs(FILENAME);
  if (!bifs) {
    cerr << "Error: unable to open " << FILENAME << endl;
    cerr << "       compile wdouble.cpp and run first\n";
    exit(1);
  }
  double d;
  long count = 0;
  cout.precision(8);
  bifs >> d;
  while (!bifs.eof()) {
    cout << ++count << ": " << d << endl;
    bifs >> d;
  }
  return 0;
}

void bifstream::readBytes(void *p, int len)
{
  if (!p) return;
  if (len <= 0) return;
  read((char *)p, len);
}
