// -*- C++ -*-
// $Header: /home/weigel/cvs/simtools/parse_opts/parse_opts.cc,v 1.21 2009-02-19 12:47:14 weigel Exp $ 

#include <machine_specific.h>
#include "my_throw.h"
#include "parse_opts.h"
#include "mystring.h"
#include <cmath>
#include <sstream>
#include <unistd.h>
#include <ctype.h>
#include <stdlib.h>
#include <libgen.h>
#include <string.h>

//==============================================================================

void Options_Parser::print_usage(bool ex)
{
  if(usage.length()) {
    cerr<<usage;
    if(ex) exit(1);
    return;
  }

  Simple_Option* opt;
  string offset(8+my_name.length(),' ');
  string line;
  unsigned int maxlength = 0;
  usage += my_name+" - "+purpose+"\n";
  line  += "usage: "+my_name+" ";
  for(MI i = regoptions.begin(); i != regoptions.end(); ++i) {
    opt = (*i).second;
    if(opt->must == false) line += "[";
    line += string("-")+(*i).first;
    if(opt->has_args) {
      maxlength = ((opt->name).length()+2*(opt->has_args%2) > maxlength) ? 
	(opt->name).length()+2*(opt->has_args%2) : maxlength;
      if(opt->has_args % 2) line += "[";
      else line+= " ";
      line += opt->name;
      if(opt->has_args % 2) line += "]";
    }
    if(opt->must == false) line += "]";
    line += " ";
    if(line.length() > 70 && i != regoptions.end()) {
      usage += line + "\n";
      line = offset;
    }
  }
  usage += line + files + "\n";
  for(MI i = regoptions.begin(); i != regoptions.end(); ++i) {
    opt = (*i).second;
    if(opt->ishelp) continue;
    line = string(" -") + (*i).first;
    if(opt->has_args) {
      if(opt->has_args % 2) line += "[";
      else line += " ";
      line += opt->name;
      if(opt->has_args % 2) line += "] ";
    }
    line += string(maxlength-2*(opt->has_args%2)-(opt->name).length()-(opt->has_args!=0)+4,' ');
    for(unsigned int j=0; j < (opt->use).length(); ++j)
      if((opt->use)[j] == '\n' || (isspace(opt->use[j]) && line.length() > 70)) {
	usage += line + "\n";
	line = string(maxlength+8,' ');
      }
      else line += (opt->use)[j];
    usage += line + "\n";
  }
  cerr<<usage;
  if(ex) exit(1);
}

//==============================================================================

void Options_Parser::parse_inclexcl(const char *str, vector<Simple_Option*> &arr)
{
  if(!str) return;
  arr.clear();
  for(const char *c = str; *c != '\0'; ++c) {
    if(!regoptions[*c]) regoptions[*c] = new Simple_Option(*c);
    arr.push_back(regoptions[*c]);
  }
}

//==============================================================================

void Options_Parser::option(char cc, bool* res, const char *u, bool mst, const char *incl, const char *excl)
{
  if(!regoptions[cc]) regoptions[cc] = new Simple_Option(cc);
  optstring += cc;
  regoptions[cc]->result = res;
  regoptions[cc]->use = u;
  regoptions[cc]->must = mst;
  regoptions[cc]->has_args = 0;
  regoptions[cc]->print_length = 3;
  parse_inclexcl(incl,regoptions[cc]->include);
  parse_inclexcl(excl,regoptions[cc]->exclude);
}

//==============================================================================

void Options_Parser::option(char cc, bool *res, bool mandatory, const char *argname, const char* u,
			    bool mst, const char *incl, const char *excl)
{
  if(parsed) my_throw("Options_Parser: options have already been parsed\n");
  option(cc,res,u,mst,incl,excl);
  regoptions[cc]->has_args = mandatory + 1;
  if(mandatory) optstring += ":";
  else optstring += "::";
  regoptions[cc]->print_length = 4 + strlen(argname) + (1-mandatory)*4;
  regoptions[cc]->name = argname;
}

//==============================================================================

void Options_Parser::option(char cc, bool *res, bool mandatory, int num, const char sep, const char *argname,
			    const char* u, bool mst, const char *incl, const char *excl)
{
  option(cc,res,u,mst,incl,excl);
  regoptions[cc]->has_args = mandatory + 3;
  if(mandatory) optstring += ":";
  else optstring += "::";
  regoptions[cc]->length = num;
  regoptions[cc]->separator = sep;
  regoptions[cc]->name = argname;
}

//==============================================================================

void Options_Parser::option(char cc, bool *res, bool mandatory, char *other_num, const char sep, const char* argname,
			    const char* u, bool mst, const char *incl, const char *excl)
{
  option(cc,res,u,mst,incl,excl);
  regoptions[cc]->has_args = mandatory + 3;
  if(mandatory) optstring += ":";
  else optstring += "::";
  regoptions[cc]->length = 0;
  regoptions[cc]->other_num = *other_num;
  if(!regoptions[*other_num]) my_throw("Options_Parser: 'other_num' option not registered\n");
  if(regoptions[*other_num]->has_args > 2) my_throw("Options_Parser: 'other_num' must be option with scalar arg\n");
  regoptions[cc]->separator = sep;
  regoptions[cc]->name = argname;
}

//==============================================================================

void Options_Parser::help_option(char cc)
{
  if(!regoptions[cc]) regoptions[cc] = new Simple_Option(cc);
  optstring += cc;
  regoptions[cc]->ishelp = true;
}

//==============================================================================

void Options_Parser::split_string(char sep, string &orig, vector<string> &dest)
{
  unsigned int pos = 0;
  string::const_iterator i = orig.begin();
  for(;;) {
    for(; i != orig.end() && *i != sep; ++i) dest[pos] += *i;
    if(i == orig.end()) break;
    ++i; ++pos;
    dest.resize(pos+1);
  }
}

//==============================================================================

void Options_Parser::parse(int *argc, char **argv)
{
  if(parsed) return;

  extern char *optarg;
  Simple_Option* opt;
  char result;

  for(int i = 0; i < *argc; ++i) cmdline += argv[i]+string(" ");
  
  my_name = basename(argv[0]);

  result = getopt(*argc,argv,optstring.c_str());
  while(result != EOF) {
    if(result == '?' || regoptions[result]->ishelp) print_usage();
    opt = regoptions[result];
    *(opt->result) = true;
    opt->_isset = true;
    if(optarg != NULL) {
      if(opt->has_args < 3) (opt->argument) = string(optarg);
      else ((opt->arguments).push_back(string(optarg)));
    }
    result = getopt(*argc,argv,optstring.c_str());
  }

  int mynum = 0;
  string mystr;
  
  for(MI i = regoptions.begin(); i != regoptions.end(); ++i) {
    opt = (*i).second;
    if(opt->_isset) {
      for(vector<Simple_Option*>::iterator j=(opt->include).begin(); j != (opt->include).end(); ++j)
	if(!((*j)->_isset)) die(string("option ")+(*i).first+" implies option "+(*j)->c);
      for(vector<Simple_Option*>::iterator j=(opt->exclude).begin(); j != (opt->exclude).end(); ++j)
	if(((*j)->_isset)) die(string("option ")+(*i).first+" excludes option "+(*j)->c);
      if(opt->has_args > 2) {
	if(opt->length == 0) {
	  if(!regoptions[opt->other_num]->_isset) die(string("option ")+(*i).first+" implies option "+regoptions[opt->other_num]->c);
	  mynum = atoi((regoptions[opt->other_num]->argument).c_str());
	  if(mynum <= 0) die(string("option ")+opt->other_num+" requires positive integer argument");
	}
	else mynum = opt->length;
	mystr = (opt->arguments)[0];
	(opt->arguments)[0] = "";
	split_string(opt->separator,mystr,opt->arguments);
	if(opt->length != -1 && ((int)(opt->arguments).size() < mynum 
				 || (int)(opt->arguments).size() > mynum))
	  die(mystring("option ")+(*i).first+" requires "+mystring(mynum)+" arguments");
      }
    }
    else if(opt->must) die(string("option ")+(*i).first+" must be given");
  }
  parsed = true;
}

//==============================================================================

void Options_Parser::convert(char cc, int &res, int min, int max)
{
  if(!parsed) my_throw("Options_Parser: must parse options before converting\n");
  if(!regoptions[cc]) my_throw("Options_Parser: unregistered option\n");
  if(!regoptions[cc]->_isset) return;
  string *s = &(regoptions[cc]->argument);
  if(s->empty()) return;
  char *endptr;
  const char *nptr = s->c_str();
  res = strtol(nptr, &endptr, 10);
  if(*nptr == '\0' || *endptr != '\0') die(string("option ")+cc+" requires integer argument");
  if(res < min || res > max) die(string("argument out of range for option ")+cc);
}

//==============================================================================

void Options_Parser::convert(char cc, vector<int> &res, int min, int max)
{
  if(!parsed) my_throw("Options_Parser: must parse options before converting\n");
  if(!regoptions[cc]) my_throw("Options_Parser: unregistered option\n");
  if(!regoptions[cc]->_isset) return;
  vector<string> *s = &(regoptions[cc]->arguments);
  if(s->empty()) return;
  res.clear();
  unsigned int count = 0;
  for(vector<string>::iterator i = s->begin(); i != s->end(); ++i, ++count) {
    const char *nptr = i->c_str();
    char *endptr;
    res.push_back(strtol(nptr, &endptr, 10));
    if(*nptr == '\0' || *endptr != '\0') die(string("option ")+cc+" requires integer argument");
    if(res.back() < min || res.back() > max) die(string("argument out of range for option ")+cc);
  }
}

//==============================================================================

void Options_Parser::convert(char cc, double &res, double min, double max)
{
  if(!parsed) my_throw("Options_Parser: must parse options before converting\n");
  if(!regoptions[cc]) my_throw("Options_Parser: unregistered option\n");
  if(!regoptions[cc]->_isset) return;
  string *s = &(regoptions[cc]->argument);
  if(s->empty()) return;
  char *endptr;
  const char *nptr = s->c_str();
  res = strtod(nptr, &endptr);
  if(*nptr == '\0' || *endptr != '\0') die(string("option ")+cc+" requires float argument");
  if(res < min || res > max) die(string("argument out of range for option ")+cc);
}

//==============================================================================

void Options_Parser::convert(char cc, vector<double> &res, double min, double max)
{
  if(!parsed) my_throw("Options_Parser: must parse options before converting\n");
  if(!regoptions[cc]) my_throw("Options_Parser: unregistered option\n");
  if(!regoptions[cc]->_isset) return;
  vector<string> *s = &(regoptions[cc]->arguments);
  if(s->empty()) return;
  res.clear();
  for(vector<string>::iterator i = s->begin(); i != s->end(); ++i) {
    const char *nptr = i->c_str();
    char *endptr;
    res.push_back(strtod(nptr, &endptr));
    if(*nptr == '\0' || *endptr != '\0') die(string("option ")+cc+" requires float argument");
    if(res.back() < min || res.back() > max) die(string("argument out of range for option ")+cc);
  }
}

//==============================================================================

void Options_Parser::convert(char cc, char &res, const char* allowed)
{
  if(!parsed) my_throw("Options_Parser: must parse options before converting\n");
  if(!regoptions[cc]) my_throw("Options_Parser: unregistered option\n");
  if(!regoptions[cc]->_isset) return;
  string *s = &(regoptions[cc]->argument);
  if(s->empty()) return;
  if(s->size() > 1) {
    cout<<"found argument '"<<*s<<"'"<<endl;
    die(string("option ")+cc+" requires char argument");
  }
  res = (*s)[0];
  const char *i = allowed;
  for(; *i != '\0' && *i != res; ++i);
  if(*i == '\0') die(string("bad char argument for option ")+cc);
}

//==============================================================================

void Options_Parser::convert(char cc, string &res)
{
  if(!parsed) my_throw("Options_Parser: must parse options before converting\n");
  if(!regoptions[cc]) my_throw("Options_Parser: unregistered option\n");
  if(!regoptions[cc]->_isset) return;
  res = regoptions[cc]->argument;
}

//==============================================================================

void Options_Parser::convert(char cc, vector<string> &res)
{
  if(!parsed) my_throw("Options_Parser: must parse options before converting\n");
  if(!regoptions[cc]) my_throw("Options_Parser: unregistered option\n");
  if(!regoptions[cc]->_isset) return;
  res = regoptions[cc]->arguments;
}

//==============================================================================

Options_Parser::~Options_Parser()
{
  for(MI p = regoptions.begin(); p != regoptions.end(); ++p) delete p->second;
}

//==============================================================================
