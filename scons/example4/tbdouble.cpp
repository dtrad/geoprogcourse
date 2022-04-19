//==============================================================
// tbdouble.cpp -- Writes and reads double values in a binary file
// Time-stamp: <1999-06-29 15:13:19 tswan>
// To compile:
//   g++ -c bstream.cpp
//   g++ -o tbdouble tbdouble.cpp bstream.o
// To run:
//   ./tbdouble
// Copyright (c) 1999 by Tom Swan. All rights reserved.
//==============================================================

#include <iostream.h>
#include <stdlib.h>    // Need exit()
#include "bstream.h"

#define FILENAME "tbdouble.dat"

int main()
{

// Construct binary output file
  bofstream bofs(FILENAME);
  if (!bofs) {
    cerr << "Error: unable to write to " << FILENAME << endl;
    exit(1);
  }

// Write values and close file
  double d = 3.14159;
  bofs << d << d * d << d * d * d;
  bofs.close();

// Construct binary input file
  bifstream bifs(FILENAME);
  if (!bifs) {
    cerr << "Error: unable to open " << FILENAME << endl;
    exit(2);
  }

// Read and display values from file
  int count = 0;
  cout.precision(8);
  bifs >> d;
  while (!bifs.eof()) {
    cout << ++count << ": " << d << endl;
    bifs >> d;
  }

  return 0;
}

