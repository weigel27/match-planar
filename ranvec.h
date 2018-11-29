// -*- C++ -*-
// $Header: /home/weigel/cvs/simtools/ranvec/ranvec.h,v 1.12 2008-08-11 12:19:47 weigel Exp $

#include <cmath>
#include <string>

using namespace std;

#ifndef _RANVEC_H
#define _RANVEC_H

//==============================================================================
// interface class

class ranvec {
protected:
  double logfactorial(int);
public:
  virtual double floating() = 0;
  inline double gaussian(double, double);
  inline int poisson(double);
  inline int poisson_small(double);
  int poisson_large(double);
};

//==============================================================================
// improved R250 random number generator

class r250 : public ranvec {
  static const int nrand, biginteger, nwarm, bigmagic1;
  static const int smallmagic1, bigmagic2, smallmagic2, nbit;
  static const double bigfloat, multiply, factor;

  int count;
  int *rand_w_array1;
  int *rand_w_array2;
  double *random_numbers;

  void fill_arrays();
public:
  r250() { }
  r250(int iseed) { setup(iseed); } 
  ~r250();

  void setup(int iseed);
  inline int integer(int);
  inline unsigned int integer(unsigned int);
  inline double floating();
  
  void write(const char*);
  void read(const char*);
  void write(const string &s) { write(s.c_str()); }
  void read(const string &s) { read(s.c_str()); }
};

//==============================================================================
// Mersenne twister

class mersenne : public ranvec {
  static const int N = 624;
  static const int M = 397;;
  static const unsigned long MATRIX_A;   /* constant vector a */
  static const unsigned long UPPER_MASK; /* most significant w-r bits */
  static const unsigned long LOWER_MASK; /* least significant r bits */
  unsigned long mt[N]; /* the array for the state vector  */
  int mti; /* mti==N+1 means mt[N] is not initialized */
  
  void init_genrand(unsigned long s);
  unsigned long genrand_int32();
public:
  mersenne() { mti = N+1; }
  mersenne(unsigned long init[], int length) { setup(init, length); } 
  ~mersenne() { }

  void setup(unsigned long init[], int length);
  inline int integer(int);
  inline unsigned int integer(unsigned int);
  inline double floating();
  
  void write(const char*);
  void read(const char*);
  void write(const string &s) { write(s.c_str()); }
  void read(const string &s) { read(s.c_str()); }
};

//==============================================================================
// inline methods for interface

inline double ranvec::gaussian(double mean, double sigma)
{
  return(sigma*sqrt(-2*log(floating()))*cos(2*M_PI*floating())+mean);
}

//==============================================================================

inline int ranvec::poisson(double lambda)
{
  return (lambda < 30.0) ? poisson_small(lambda) : poisson_large(lambda);	
}

//==============================================================================

inline int ranvec::poisson_small(double lambda)
{
  // Algorithm due to Donald Knuth, 1969.
  double p = 1.0, L = exp(-lambda);
  int k = 0;
  do
  {
    k++;
    p *= floating();
  } while (p > L);
  return k - 1;
}


//==============================================================================
// inline methods for r250

//==============================================================================

inline double r250::floating()
{
  do {
    count++; 
    if(count>=nrand) {
      fill_arrays();
      count=0;
    }
    if(random_numbers[count]>0 && random_numbers[count]<1) break;
  } while(1);

  return(random_numbers[count]);

}

//==============================================================================

inline int r250::integer(int max) {
  return(static_cast<int>(floating()*max));
}

//==============================================================================

inline unsigned int r250::integer(unsigned int max) {
  return(static_cast<unsigned int>(floating()*max));
}

//==============================================================================
//==============================================================================
// inline methods for mersenne

//==============================================================================

inline double mersenne::floating()
{
  return genrand_int32()*(1.0/4294967296.0);
}

//==============================================================================

inline int mersenne::integer(int max) {
  return(static_cast<int>(floating()*max));
}

//==============================================================================

inline unsigned int mersenne::integer(unsigned int max) {
  return(static_cast<unsigned int>(floating()*max));
}

//==============================================================================

#ifdef _INCLUDE_ALL
#include "ranvec.cc"
#endif

#endif
