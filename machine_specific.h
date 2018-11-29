#ifndef _MACHINE_SPECIFIC_H
#define _MACHINE_SPECIFIC_H

//template<class T> double fabs(T val) { return fabs((double)val); }

#ifndef M_PI
#define M_PI    3.14159265358979323846
#endif

#define GCC_VERSION (__GNUC__ * 10000 \
                     + __GNUC_MINOR__ * 100 \
                     + __GNUC_PATCHLEVEL__)

#define MAXFLOAT 3.40282347e+38F

#ifdef __SR8000
  #define ios_out ios_base::out
  #define ios_in ios_base::in
  #define ios_app ios_base::app
  #define ios_binary ios_base::binary
  #define ios_beg ios_base::beg
  #define ios_end ios_base::end
  #define ios_cur ios_base::cur
  #define ios_scientific ios_base::scientific
  #define ios_floatfield ios_base::floatfield
#elif __GNUC__ == 3 || __INTEL_COMPILER
  #define ios_out ios_base::out
  #define ios_in ios_base::in
  #define ios_app ios_base::app
  #define ios_binary ios_base::binary
  #define ios_beg ios_base::beg
  #define ios_end ios_base::end
  #define ios_cur ios_base::cur
  #define ios_scientific ios_base::scientific
  #define ios_floatfield ios_base::floatfield
#else
  #define ios_out ios::out
  #define ios_in ios::in
  #define ios_app ios::app
  #define ios_binary ios::binary
  #define ios_beg ios::beg
  #define ios_end ios::end
  #define ios_cur ios::cur
  #define ios_scientific ios::scientific
  #define ios_floatfield ios::floatfield
#endif

#define BOOST_EXCEPTION_DISABLE

#ifdef __cplusplus
#include <iostream>
#if defined(_CRAY)
  typedef ios::open_mode ios_open_mode;
  typedef ios::seek_dir ios_seek_dir;
  namespace std { };
#elif defined(__SR8000)
  using namespace std;
  typedef ios_base::open_mode ios_open_mode;
  typedef ios_base::seekdir ios_seek_dir;
#elif __GNUC__ >= 3
  using namespace std;
  typedef std::ios_base::openmode ios_open_mode;
  typedef std::ios_base::seekdir ios_seek_dir;
#else
  typedef ios::openmode ios_open_mode;
  typedef ios::seekdir ios_seek_dir;
#endif
#endif

#endif
