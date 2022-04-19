//==============================================================
// bstream.cpp -- Implements the bofstream and bifstream classes
// Time-stamp: <1999-06-29 14:19:46 tswan>
// To compile:
//   g++ -c bstream.cpp
// Copyright (c) 1999 by Tom Swan. All rights reserved.
//==============================================================

#include "bstream.h"

void bofstream::writeBytes(const void *p, int len)
{
  if (!p) return;
  if (len <= 0) return;
  write((char *)p, len);
}

void bifstream::readBytes(void *p, int len)
{
  if (!p) return;
  if (len <= 0) return;
  read((char *)p, len);
}
