//==============================================================
// wline.cpp -- Writes a text file a line at a time
// Time-stamp: <1999-06-29 10:47:22 tswan>
// To compile:
//   g++ -o wline wline.cpp
// To run:
//   ./wline filename.txt
// Copyright (c) 1999 by Tom Swan. All rights reserved.
//==============================================================

#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>
#include <string.h>
#include <string>

#define STR "2: Another literal string"

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
  ofs << "1: A literal string" << endl;
  ofs.write(STR, strlen(STR));
  ofs << endl;
  char *c = "String addressed by pointer";
  ofs << "3: " << c << endl;
  string string_object("4: A string object");
  ofs << string_object << endl;
  return 0;
}
