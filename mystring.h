// -*- C++ -*-
// $Header: /home/weigel/cvs/simtools/mystring/mystring.h,v 1.10 2008-08-11 12:19:47 weigel Exp $

#include <stdio.h>
#include <string>

using namespace std;

#ifndef _MYSTRING_H
#define _MYSTRING_H

class mystring : public string {
public:
  mystring() : string() { }
  mystring(const string &other) : string(other) { }
  mystring(const char *other) : string(other) { }
  mystring(int d) { char buf[20]; sprintf(buf,"%d",d); *this = buf; }
  mystring(unsigned int d) { char buf[20]; sprintf(buf,"%d",d); *this = buf; }
  mystring(char c) { char buf[2]; sprintf(buf,"%s",&c); *this = buf; }
  mystring(double d) { char buf[20]; sprintf(buf,"%f",d); *this = buf; }

  template<class A> mystring& operator<<(const A& a) { *this += mystring(a); return *this; }
};

//template<class A> mystring& operator<<(mystring& s, const A& a) { s += mystring(a); return s; }

#endif
