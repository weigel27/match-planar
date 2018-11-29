// -*- C++ -*-
// $Header: /home/weigel/cvs/simtools/my_throw/my_throw.h,v 1.8 2008-08-11 12:19:47 weigel Exp $

#include <iostream>
#include <stdlib.h>
#include "mystring.h"

#ifndef _MY_THROW_H_
#define _MY_THROW_H_

//==============================================================================

inline void my_throw(const char* s, bool abort = true)
{
  cerr<<s;
  if(abort) exit(1);
}

inline void my_throw(const mystring& s, bool abort = true)
{
  cerr<<s;
  if(abort) exit(1);
}

//==============================================================================

#endif
