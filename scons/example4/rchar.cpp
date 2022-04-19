//==============================================================
// rchar.cpp -- Reads a text file a char at a time
// Time-stamp: <1999-06-25 15:43:43 tswan>
// To compile:
//   g++ -o rchar rchar.cpp
// To run:
//   ./rchar rchar.cpp
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
  if (!ifs) {
    cerr << "Error: unable to open " << argv[1] << endl;
    exit(2);
  }
  char c;
  long nc = 0, nl = 0;
  while (ifs.get(c)) {
    if (c == '\n') 
      nl++;  // Count number of lines
    else 
      nc++;  // Count number of characters
    cout.put(c);
  }
  cout << endl << endl << "Total characters : " << nc;
  cout << endl << "Number of lines  : " << nl << endl;
  return 0;
}
