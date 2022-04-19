//==============================================================
// wchar.cpp -- Writes a text file a char at a time
// Time-stamp: <1999-06-29 10:03:53 tswan>
// To compile:
//   g++ -o wchar wchar.cpp
// To run:
//   ./wchar filename.txt
// Copyright (c) 1999 by Tom Swan. All rights reserved.
//==============================================================

#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
  if (argc <= 1) {
    cerr << "Error: filename missing" << endl;
    exit(1);
  }
  ifstream ifs(argv[1], ios::in);
  if (ifs) {
    cerr << "Error: " << argv[1] << " already exists" << endl;
    cerr << "       Specify a different filename" << endl;
    exit(2);
  }
  ofstream ofs(argv[1], ios::out);
  if (!ofs) {
    cerr << "Error: unable to write to " << argv[1] << endl;
    exit(3);
  }
  ofs << "1: A string" << endl;
  ofs.put('2');
  ofs.put(':');
  ofs.put(' ');
  ofs.put('C').put('h').put('a').put('r').put('s');
  ofs << endl;
  return 0;
}
