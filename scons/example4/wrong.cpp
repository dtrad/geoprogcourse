//==============================================================
// wrong.cpp -- Demonstrates incorrect way to write binary data
// Time-stamp: <1999-06-29 11:16:45 tswan>
// To compile:
//   g++ wrong.cpp
// To run:
//   ./a.out
// Copyright (c) 1999 by Tom Swan. All rights reserved.
//==============================================================

// *** Note: Don't use the method illustrated in this program.
//     It demonstrates the WRONG way to write binary data to
//     a file. For the correct solution, see the text and also
//     examine the files bstream.h, bstream.cpp, bfstream.h, 
//     bfstream.cpp, and associated test programs.

#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>

int main()
{
  ifstream xf("test.dat", ios::in);
  if (xf) {
    cerr << "File already exists" << endl;
    cerr << "Delete test.dat and try again" << endl;
    exit(1);
  }
  ofstream ofs("test.dat", ios::out | ios::binary);
  if (ofs) {
    double d = 3.14159;
    ofs << d;       // ???
    ofs << d * d;   // ???
  }
  ofs.close();
  ifstream ifs("test.dat", ios::in | ios::binary);
  if (ifs) {
    double d;
    ifs >> d;    // ???
    cout << "d == " << d << endl;
    ifs >> d;    // ???
    cout << "d * d == " << d << endl;
  }
  return 0;
}
