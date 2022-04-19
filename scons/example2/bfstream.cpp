//==============================================================
// bfstream.cpp -- Implements the bfstream (binary file) class
// Time-stamp: <1999-06-29 14:49:23 tswan>
// To compile:
//   g++ -c bfstream.cpp
// Copyright (c) 1999 by Tom Swan. All rights reserved.
//==============================================================

#include "bfstream.h"

void bfstream::writeBytes(const void *p, int len)
{
  if (!p) return;
  if (len <= 0) return;
  write((char *)p, len);
}

void bfstream::readBytes(void *p, int len)
{
  if (!p) return;
  if (len <= 0) return;
  read((char *)p, len);
}
