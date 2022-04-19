//==============================================================
// tanyclass.cpp -- Implements the test TAnyClass class
// Time-stamp: <1999-06-29 15:02:48 tswan>
// To compile:
//   g++ -c tanyclass.cpp
// Copyright (c) 1999 by Tom Swan. All rights reserved.
//==============================================================

#include <iostream>
#include "tanyclass.h"

// Implements overloaded ostream output operator
ostream & operator<< (ostream &os, const TAnyClass &c)
{
  os << "-- TAnyClass object --" << endl;
  os << "x == " << c.x << "; ";
  os << "y == " << c.y << endl;
  return os;
}

// Implements overloaded istream input operator
istream & operator>> (istream &is, TAnyClass &c)
{
  cout << " Enter value for X: ";  
  is >> c.x;
  cout << " Enter value for Y: ";
  is >> c.y;
  return is;
}
