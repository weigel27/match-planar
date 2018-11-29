// -*- C++ -*-
// $Header: /home/weigel/cvs/simtools/parse_opts/parse_opts.h,v 1.21 2009-02-19 12:47:14 weigel Exp $ 

#include <machine_specific.h>
#include <cmath>
#include <limits>

using namespace std;

#include "mystring.h"
#include <strings.h>
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <stdlib.h>

#ifndef _PARSE_OPTS_H_
#define _PARSE_OPTS_H_

//==============================================================================

class Options_Parser;

class Simple_Option {
protected:
  char c;
  bool _isset;
  bool *result;
  string use;
  bool must;
  unsigned int has_args;
  vector<Simple_Option*> include;
  vector<Simple_Option*> exclude;
  string name;
  unsigned int print_length;
  bool ishelp;
  
  string argument;
  vector<string> arguments;
  int length;
  char other_num;
  char separator;
  
  Simple_Option(char cc) : c(cc), _isset(false), must(false), has_args(0), ishelp(false) { }
  ~Simple_Option() { }

  friend class Options_Parser;
};

//==============================================================================

class Options_Parser {
  map<char,Simple_Option*> regoptions;
  string purpose, usage, files;
  string my_name;
  string optstring;
  bool parsed;
  
  void parse_inclexcl(const char *str, vector<Simple_Option*> &arr);
  void print_usage(bool ex = true);
  void die(mystring mess) { cerr<<my_name<<": "<<mess<<endl; exit(1); }
  void split_string(char sep, string &orig, vector<string> &dest);
public:
  Options_Parser(const char *p, const char* f) { purpose = p; files = f; parsed = false; }
  ~Options_Parser();

  typedef map<char,Simple_Option*>::iterator MI;
  string cmdline;
  
  //int first_non_opt() const { extern int optind; return optind; }

  void option(char cc, bool *res, const char* u, bool mst=false, const char *incl=NULL, const char *excl=NULL);
  void help_option(char cc);
  void option(char cc, bool *res, bool mandatory, const char *argname, const char* u,
	      bool mst=false, const char *incl=NULL, const char *excl=NULL);
  void option(char cc, bool *res, bool mandatory, int num, const char sep,
	      const char *argname, const char* u, bool mst=false, const char *incl=NULL,
	      const char *excl=NULL);
  void option(char cc, bool *res, bool mandatory, char *other_num, const char sep,
	      const char *argname, const char* u, bool mst=false, const char *incl=NULL,
	      const char *excl=NULL);
  
  void parse(int *argc, char **argv);
  void convert(char cc, int &res, int min = numeric_limits<int>::min(), int max = numeric_limits<int>::max());
  void convert(char cc, vector<int> &res, int min = numeric_limits<int>::min(), int max = numeric_limits<int>::max());
  // FIXME: g++ doesn't allow to use +/- HUGE_VAL as default parameter
  void convert(char cc, double &res, double min = numeric_limits<int>::min(), double max = numeric_limits<int>::max());
  void convert(char cc, vector<double> &res, double min = numeric_limits<int>::min(), double max = numeric_limits<int>::max());
  void convert(char cc, string &res);
  void convert(char cc, vector<string> &res);
  void convert(char cc, char &res, const char* allowed);
};

//==============================================================================

#ifdef _INCLUDE_ALL
#include "parse_opts.cc"
#endif

#endif
