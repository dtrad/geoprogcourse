//==============================================================
// rline.cpp -- Reads a text file a line at a time
// Time-stamp: <1999-06-25 15:46:04 tswan>
// To compile:
//   g++ -o rline rline.cpp
// To run:
//   ./rline rline.cpp
// Copyright (c) 1999 by Tom Swan. All rights reserved.
//==============================================================

#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include <string.h>

#define BUFLEN 128

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
  char buffer[BUFLEN];
  long nc = 0, nl = 0;
  while (!ifs.eof()) {
    ifs.getline(buffer, sizeof(buffer), '\n');
    if (!(ifs.eof() && strlen(buffer) == 0)) {
      nc += strlen(buffer);
      nl++;
      cout << buffer << endl;
    }
  }
  cout << endl << endl << "Total characters : " << nc;
  cout << endl << "Number of lines  : " << nl << endl;
  return 0;
}
