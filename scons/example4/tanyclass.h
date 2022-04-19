//==============================================================
// tanyclass.h -- Declarations for TAnyClass test class
// Time-stamp: <1999-06-29 15:03:10 tswan>
// Copyright (c) 1999 by Tom Swan. All rights reserved.
//==============================================================

#ifndef __tanyclass_H
#define __tanyclass_H   // Prevent multiple #includes

#include <iostream>
using namespace std;
class TAnyClass {
private:
  int x;  // Private class data members
  int y;
public:
  friend ostream & operator<< (ostream &, const TAnyClass &);
  friend istream & operator>> (istream &, TAnyClass &);
public:
  TAnyClass(): x(0), y(0) { }
  TAnyClass(int X, int Y): x(X), y(Y) { }
};

#endif  // __tanyclass_H
