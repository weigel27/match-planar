// -*- C++ -*-
// $Header: /home/weigel/cvs/simtools/alloc/ccalloc.h,v 1.6 2008-08-11 12:19:46 weigel Exp $

#ifndef _CCALLOC_H
#define _CCALLOC_H

#include <string.h>
#include "my_throw.h"

/*!
  \file Heap allocation for C++ progs
*/

//==============================================================================
//! Generic programming template to generate expressions such as '***T'
//! with variable number of '*'

template<typename T, unsigned int n> struct _add_pointer {
  typedef typename _add_pointer<T, n-1>::type* type;
};

template<typename T> struct _add_pointer<T,0> {
  typedef T type;
};

//! Generic programming template to remove a number of pointers from a type
template<typename T, unsigned int n> struct _remove_pointer {
  typedef T type;
};

template<typename T, unsigned int n> struct _remove_pointer<T*, n> {
  typedef typename _remove_pointer<T, n-1>::type type;
};

//==============================================================================
//! Allocation

template<typename T, unsigned int d> typename _add_pointer<T, d>::type cc_alloc(const unsigned int n[]);

//! Unfortunately, partial specialization of function templates is not allowed
//! but needed here to avoid infinite regress. Thus, a dummy function has to be
//! fully specialized for all intended types
template<> _add_pointer<unsigned int, 0>::type cc_alloc<unsigned int,0>(const unsigned int n[]);
template<> _add_pointer<short int, 0>::type cc_alloc<short int,0>(const unsigned int n[]);
template<> _add_pointer<unsigned short int, 0>::type cc_alloc<unsigned short int,0>(const unsigned int n[]);
template<> _add_pointer<long int, 0>::type cc_alloc<long int,0>(const unsigned int n[]);
template<> _add_pointer<long long int, 0>::type cc_alloc<long long int,0>(const unsigned int n[]);
template<> _add_pointer<char, 0>::type cc_alloc<char,0>(const unsigned int n[]);
template<> _add_pointer<unsigned char, 0>::type cc_alloc<unsigned char,0>(const unsigned int n[]);
template<> _add_pointer<int, 0>::type cc_alloc<int,0>(const unsigned int n[]);
template<> _add_pointer<float, 0>::type cc_alloc<float,0>(const unsigned int n[]);
template<> _add_pointer<double, 0>::type cc_alloc<double,0>(const unsigned int n[]);
template<> _add_pointer<long double, 0>::type cc_alloc<long double,0>(const unsigned int n[]);

//! The real allocation template
template<typename T, unsigned int d> typename _add_pointer<T, d>::type cc_alloc(const unsigned int n[])
{
  typedef typename _add_pointer<T, d>::type pointer_t;
  typedef typename _add_pointer<T, d-1>::type removed_pointer_t;
  pointer_t pointer;

  //if(n[d-1] < 1) my_throw("cc_alloc(): illegal array dimension\n");
  //pointer = (pointer_t)calloc(n[d-1], sizeof(typename _remove_pointer<pointer_t,1>::type));
  pointer = (pointer_t)calloc(n[d-1], sizeof(removed_pointer_t));
  if(!pointer) my_throw("cc_alloc(): out of memory\n");
  //fprintf(stdout,"allocated %d elements of size %d at address %p\n",n[d-1],sizeof(removed_pointer_t),pointer);

  if(d > 1)
    for(unsigned int i = 0; i < n[d-1]; ++i) pointer[i] = cc_alloc<T, d-1>(n);

  return(pointer);
}

//! Specializations for 1- and 2-dimensional arrays
//template<typename T> typename _add_pointer<T, 1>::type cc_alloc(unsigned int n) { return cc_alloc<T, 1>((const unsigned int[]){n}); }
//template<typename T> typename _add_pointer<T, 2>::type cc_alloc(unsigned int n1, unsigned int n2) { return cc_alloc<T, 2>((const unsigned int[]){n1, n2}); }

//==============================================================================
//! Deallocation

template<typename T, unsigned int d> void cc_free(typename _add_pointer<T, d>::type p, const unsigned int n[]);

template<> void cc_free<unsigned int,0>(_add_pointer<unsigned int, 0>::type p, const unsigned int n[]);
template<> void cc_free<short int,0>(_add_pointer<short int, 0>::type p, const unsigned int n[]);
template<> void cc_free<unsigned short int,0>(_add_pointer<unsigned short int, 0>::type p, const unsigned int n[]);
template<> void cc_free<long int,0>(_add_pointer<long int, 0>::type p, const unsigned int n[]);
template<> void cc_free<long long int,0>(_add_pointer<long long int, 0>::type p, const unsigned int n[]);
template<> void cc_free<char,0>(_add_pointer<char, 0>::type p, const unsigned int n[]);
template<> void cc_free<unsigned char,0>(_add_pointer<unsigned char, 0>::type p, const unsigned int n[]);
template<> void cc_free<int,0>(_add_pointer<int, 0>::type p, const unsigned int n[]);
template<> void cc_free<float,0>(_add_pointer<float, 0>::type p, const unsigned int n[]);
template<> void cc_free<double,0>(_add_pointer<double, 0>::type p, const unsigned int n[]);
template<> void cc_free<long double,0>(_add_pointer<long double, 0>::type p, const unsigned int n[]);

template<typename T, unsigned int d> void cc_free(typename _add_pointer<T, d>::type p, const unsigned int n[])
{
  if(d > 1)
    for(unsigned int i = 0; i < n[d-1]; ++i) cc_free<T, d-1>(p[i], n);

  free(p);
}

//! Specializations for 1- and 2-dimensional arrays
//template<typename T> void cc_free(typename _add_pointer<T, 1>::type p, unsigned int n) { cc_free<T, 1>(p, (const unsigned int[]){n}); }
//template<typename T> void cc_free(typename _add_pointer<T, 2>::type p, unsigned int n1, unsigned int n2) { cc_free<T, 2>(p, (const unsigned int[]){n1, n2}); }

//==============================================================================
//! Content deletion

template<typename T, unsigned int d> void cc_clear(typename _add_pointer<T, d>::type p, const unsigned int n[]);

template<> void cc_clear<unsigned int,0>(_add_pointer<unsigned int, 0>::type p, const unsigned int n[]);
template<> void cc_clear<short int,0>(_add_pointer<short int, 0>::type p, const unsigned int n[]);
template<> void cc_clear<unsigned short int,0>(_add_pointer<unsigned short int, 0>::type p, const unsigned int n[]);
template<> void cc_clear<long int,0>(_add_pointer<long int, 0>::type p, const unsigned int n[]);
template<> void cc_clear<long long int,0>(_add_pointer<long long int, 0>::type p, const unsigned int n[]);
template<> void cc_clear<char,0>(_add_pointer<char, 0>::type p, const unsigned int n[]);
template<> void cc_clear<unsigned char,0>(_add_pointer<unsigned char, 0>::type p, const unsigned int n[]);
template<> void cc_clear<int,0>(_add_pointer<int, 0>::type p, const unsigned int n[]);
template<> void cc_clear<float,0>(_add_pointer<float, 0>::type p, const unsigned int n[]);
template<> void cc_clear<double,0>(_add_pointer<double, 0>::type p, const unsigned int n[]);
template<> void cc_clear<long double,0>(_add_pointer<long double, 0>::type p, const unsigned int n[]);

template<typename T, unsigned int d> void cc_clear(typename _add_pointer<T, d>::type p, const unsigned int n[])
{
  if(d > 1)
    for(unsigned int i = 0; i < n[d-1]; ++i) cc_clear<T, d-1>(p[i], n);
  else  
    bzero(p,sizeof(*p)*n[d-1]);  
}

//! Specializations for 1- and 2-dimensional arrays
//template<typename T> void cc_clear(typename _add_pointer<T, 1>::type p, unsigned int n) { cc_clear<T, 1>(p, (const unsigned int[]){n}); }
//template<typename T> void cc_clear(typename _add_pointer<T, 2>::type p, unsigned int n1, unsigned int n2) { cc_clear<T, 2>(p, (const unsigned int[]){n1, n2}); }

//==============================================================================
//! Copying

template<typename T, unsigned int d> void cc_copy(typename _add_pointer<T, d>::type src, typename _add_pointer<T, d>::type dst, const unsigned int n[]);

template<> void cc_copy<unsigned int, 0>(_add_pointer<unsigned int, 0>::type src, _add_pointer<unsigned int, 0>::type dst, const unsigned int n[]);
template<> void cc_copy<short int, 0>(_add_pointer<short int, 0>::type src, _add_pointer<short int, 0>::type dst, const unsigned int n[]);
template<> void cc_copy<unsigned short int, 0>(_add_pointer<unsigned short int, 0>::type src, _add_pointer<unsigned short int, 0>::type dst, const unsigned int n[]);
template<> void cc_copy<long int, 0>(_add_pointer<long int, 0>::type src, _add_pointer<long int, 0>::type dst, const unsigned int n[]);
template<> void cc_copy<long long int, 0>(_add_pointer<long long int, 0>::type src, _add_pointer<long long int, 0>::type dst, const unsigned int n[]);
template<> void cc_copy<char, 0>(_add_pointer<char, 0>::type src, _add_pointer<char, 0>::type dst, const unsigned int n[]);
template<> void cc_copy<unsigned char, 0>(_add_pointer<unsigned char, 0>::type src, _add_pointer<unsigned char, 0>::type dst, const unsigned int n[]);
template<> void cc_copy<int, 0>(_add_pointer<int, 0>::type src, _add_pointer<int, 0>::type dst, const unsigned int n[]);
template<> void cc_copy<float, 0>(_add_pointer<float, 0>::type src, _add_pointer<float, 0>::type dst, const unsigned int n[]);
template<> void cc_copy<double, 0>(_add_pointer<double, 0>::type src, _add_pointer<double, 0>::type dst, const unsigned int n[]);
template<> void cc_copy<long double, 0>(_add_pointer<long double, 0>::type src, _add_pointer<long double, 0>::type dst, const unsigned int n[]);

template<typename T, unsigned int d> void cc_copy(typename _add_pointer<T, d>::type src, typename _add_pointer<T, d>::type dst, const unsigned int n[])
{
  if(d > 1)
    for(unsigned int i = 0; i < n[d-1]; ++i) cc_copy<T, d-1>(src[i], dst[i], n);
  else
    memcpy(dst,src,sizeof(*src)*n[d-1]);
}

//! Specializations for 1- and 2-dimensional arrays
//template<typename T> void cc_copy(typename _add_pointer<T, 1>::type src, typename _add_pointer<T, 1>::type dst, unsigned int n) {cc_copy<T, 1>(src, dst, (const unsigned int[]){n}); }
//template<typename T> void cc_copy(typename _add_pointer<T, 2>::type src, typename _add_pointer<T, 2>::type dst, unsigned int n1, unsigned int n2) {cc_copy<T, 1>(src, dst, (const unsigned int[]){n1, n2}); }

//==============================================================================

#endif
