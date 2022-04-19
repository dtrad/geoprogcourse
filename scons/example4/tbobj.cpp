//==============================================================
// tbobj.cpp -- Test binary object file I/O
// Time-stamp: <1999-06-29 15:21:14 tswan>
// To compile:
//   g++ -c tanyclass.cpp
//   g++ -c bstream.cpp
//   g++ -o tbobj tbobj.cpp tanyclass.o bstream.o
// To run:
//   ./a.out
// Copyright (c) 1999 by Tom Swan. All rights reserved.
//==============================================================

#include <iostream.h>
#include <stdlib.h>     // Need exit()
#include "tanyclass.h"
#include "bstream.h"

#define FILENAME "tbobj.dat"

int main()
{
  bofstream bofs(FILENAME);
  if (!bofs) {
    cerr << "Error: unable to write to " << FILENAME << endl;
    exit(1);
  }

  TAnyClass obj;  // Default object

  // Prompt user to enter object values
  cin >> obj;

  // Write object to disk and close file
  cout << "Writing object to disk" << endl;
  bofs << obj;
  bofs.close();

  // Construct binary input file
  cout << "Opening file" << endl;
  bifstream bifs(FILENAME);
  if (!bifs) {
    cerr << "Error: unable to open " << FILENAME << endl;
    exit(2);
  }

  // Read object from file. Use a new TAnyClass
  // object to be sure values are new.
  cout << "Reading object from disk" << endl;
  TAnyClass newObject;
  bifs >> newObject;

  cout << "Object from disk:" << endl;
  cout << newObject; 

  bifs.close();
  return 0;
}
