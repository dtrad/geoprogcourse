//==============================================================
// bfstest.cpp -- Tests the bfstream binary file stream class
// Time-stamp: <1999-06-29 15:20:27 tswan>
// To compile:
//   g++ -c tanyclass.cpp
//   g++ -c bfstream.cpp
//   g++ -o bfstest bfstest.cpp tanyclass.o bfstream.o
// To run:
//   ./bfstest
// Copyright (c) 1999 by Tom Swan. All rights reserved.
//==============================================================

#include <iostream.h>
#include <stdlib.h>    // Need exit()
#include "tanyclass.h"
#include "bfstream.h"

#define FILENAME "abcde.dat"

int main()
{
  TAnyClass a(0,1), b(2,3), c(4,5), d(6,7), e(8,9);
  TAnyClass aa, bb, cc, dd, ee;

  // Construct output file stream object bfs
  bfstream bfs(FILENAME);
  if (!bfs) {
    cerr << "Error: unable to create " << FILENAME << endl;
    exit(1);
  }

  // Write records to file
  cout << "Writing five records to file" << endl;
  bfs << a << b << c << d << e;

  // Reset file to beginning, read and display records
  cout << "Reading five records from file" << endl;
  bfs.seekg(0);
  bfs >> aa >> bb >> cc >> dd >> ee;
  cout << aa << bb << cc << dd << ee;

  // Seek record number 3 and change it
  TAnyClass x(123,456);               // Construct new object
  long rn = 3;                        // Define record number
  bfs.seekp(sizeof(TAnyClass) * rn);  // Seek to record #3
  bfs << x;                           // Write object to file

  // Seek to record number 3 again and read it
  TAnyClass y;                        // Construct empty object
  bfs.seekg(sizeof(TAnyClass) * rn);  // Seek to record #3
  bfs >> y;                           // Read object from file

  cout << endl << "Changed record #3 (4th record) to" << endl;
  cout << y;

  // Make sure other records are undisturbed;
  cout << endl << "Re-reading five records" << endl;
  bfs.seekg(0);
  bfs >> aa >> bb >> cc >> dd >> ee;
  cout << aa << bb << cc << dd << ee;

  // Close file object
  bfs.close();

  return 0;
}

