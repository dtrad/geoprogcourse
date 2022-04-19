//==============================================================
// seekdemo.cpp -- Demonstrates seekg() and seekp() functions
// Time-stamp: <1999-06-30 10:19:17 tswan>
// To compile:
//   g++ -c tanyclass.cpp
//   g++ -c bstream.cpp
//   g++ -o seekdemo seekdemo.cpp
// To run:
//   ./seekdemo
// Copyright (c) 1999 by Tom Swan. All rights reserved.
//==============================================================

#include <iostream>
#include <stdlib.h>     // Need exit()
#include "bfstream.h"
#include "tanyclass.h"

#define FILENAME "./seekdemo.dat"
using namespace std;
int main()
{
  // Construct some objects to store in a file
  //
  TAnyClass a(0,1), b(2,3), c(4,5), d(6,7), e(8,9);

  // Construct objects for tests (could reuse above)
  //
  TAnyClass aa, bb, cc, dd, ee;

  // Create new or overwrite existing file
  //
  bfstream bfs(FILENAME);
  if (!bfs) {
    cerr << "Error: unable to create " << FILENAME << endl;
    exit(1);
  }

  // Write objects to disk
  //
  bfs << a << b << c << d << e;

  // Read objects from disk and display
  //
  bfs.seekg(0);
  bfs >> aa >> bb >> cc >> dd >> ee;
  cout << "Records stored in " << FILENAME << endl;
  cout << aa << bb << cc << dd << ee;

  // Seek to specific record number
  cout << "Seeking to record #3" << endl;
  streampos rn = 3;  // Record number
  bfs.seekg(sizeof(TAnyClass) * rn);
  bfs >> aa;
  cout << aa;

  // Seek back from current position
  cout << "Seeking back two records" << endl;
  rn = 2;  // Number of records to seek backwards
  bfs.seekg(-(sizeof(TAnyClass) * rn), ios::cur);
  bfs >> bb;
  cout << bb;

  // Seek backwards from end of file
  cout << "Seeking back one record from end of file" << endl;
  rn = 1;  // Number of records to seek from end of file
  bfs.seekg(-(sizeof(TAnyClass) * rn), ios::end);
  bfs >> cc;
  cout << cc;

  // Seek forward from beginning of file
  cout << "Seeking four records from beginning of file" << endl;
  rn = 4;  // Number of records to seek from beginning of file
  bfs.seekg((sizeof(TAnyClass) * rn), ios::beg);
  bfs >> dd;
  cout << dd;

  // Prompt user for new record to insert in file
  cout << "Enter values for new record" << endl;
  TAnyClass newObject;
  cin >> newObject;
  rn = 2;  // Record to overwrite
  bfs.seekp(sizeof(TAnyClass) * rn, ios::beg);
  bfs << newObject;
  
  // Re-read changed record to confirm
  //
  cout << "Reading record to confirm insertion" << endl;
  TAnyClass oldObject;  // Use a new object for the test
  rn = 2;
  bfs.seekg(sizeof(TAnyClass) * rn, ios::beg);
  bfs >> oldObject;
  cout << oldObject;  

  bfs.close();
  return 0;
}
