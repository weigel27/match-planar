// -*- C++ -*-
// $Header: /home/weigel/cvs/progs/simulation/matching/options.cc,v 1.14 2007-06-15 21:11:00 weigel Exp $

#include "match.h"
#include "parse_opts.h"

_params params;

void parse_options(int *argc, char *argv[])
{
  // default values
  params.neighbors = numeric_limits<int>::max();
  params.frustrated = numeric_limits<int>::max();
  params.max_dijkstra_cutoff = 10;
  params.min_dijkstra_cutoff = 10;
  
  vector<double> range1, range2;
  range1.resize(2); range2.resize(2);
  range2[0] = params.min_dijkstra_cutoff;
  range2[1] = params.max_dijkstra_cutoff;
  bool dummy;

  Options_Parser parser("Ising model ground states on planar graphs", "");

  parser.option('d',&dummy,true,2,':',"MIN:MAX", "range of weight cut-offs for Dijkstra");
  parser.option('n',&dummy,true,"NEIGHB", "maximum number of neighbors in Dijkstra");
  parser.option('f',&dummy,true,"FRUST", "maximum number of frustrated neighbors in Dijkstra");
  parser.option('o',&dummy,true,"PATH","filename of output spin configuration",true);  
  parser.option('G',&dummy,true,"FILE","filename of graph to use",true);

  parser.help_option('h');
  parser.help_option('?');
  parser.parse(argc,argv);

  parser.convert('n',params.neighbors,1);
  parser.convert('f',params.frustrated,1);
  parser.convert('d',range2);
  parser.convert('o',params.outprefix);  
  parser.convert('G',params.graphpath);

  params.min_dijkstra_cutoff = range2[0];
  params.max_dijkstra_cutoff = range2[1];
}

istream& operator>>(istream& s, _params& p)
{
  s>>p.neighbors>>p.frustrated;
  s>>p.min_dijkstra_cutoff>>p.max_dijkstra_cutoff;
  s>>p.outprefix>>p.graphpath;
  
  return s;
}

ostream& operator<<(ostream& s, _params& p)
{
  s<<p.neighbors<<" "<<p.frustrated<<" ";
  s<<p.min_dijkstra_cutoff<<" "<<p.max_dijkstra_cutoff<<" ";
  s<<p.outprefix<<" "<<p.graphpath<<endl;

  return s;
}
