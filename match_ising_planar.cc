// -*- C++ -*-
// $Header: /home/weigel/cvs/progs/simulation/matching/match_ising_DW.cc,v 1.7 2011-08-23 12:29:30 weigel Exp $

#include <limits>
#include <map>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/graphviz.hpp>
#include <readline/readline.h>
#include <sys/stat.h>
#include <sys/times.h>
#include <unistd.h>
#include <iomanip>
#include "stdlib.h"
#include "match.h"
#include "dijkstra.h"
#include "ccalloc.h"
#include "mystring.h"
#include "zfstream.h"
#include "ranvec.h"

// content types
typedef short int Vertex_Content;
typedef float Bond_Content;

//==============================================================================
//! Does the actual ground-state computation
/*! Given a starting configuration, subsequent
  applications of an externally supplied matching
  algorithm are used to potentially find (one of)
  the ground states of the planar XY model.
*/

template<typename Graph, typename Dual_Graph, typename ContentMap, typename BondContentMap, typename DualBondContentMap,
	 typename EdgeToDualMap, typename DualToEdgeMap, typename FrustToPlaqMap, typename PlaqToFrustMap,
	 typename DistanceMatrix, typename PredecessorMatrix, typename DVertexIndexMap, typename BrokenLabelMap,
	 typename FrustratedLabelMap>
bool match(const Graph& g, const Dual_Graph& dg, ContentMap& vertexcontent, const BondContentMap& bondcontent,
	   DualBondContentMap& dual_bondcontent, const EdgeToDualMap& edge_to_dual, const DualToEdgeMap& dual_to_edge,
	   const FrustToPlaqMap& frust_to_plaq, const PlaqToFrustMap& plaq_to_frust, DistanceMatrix& distance_matrix,
	   PredecessorMatrix& predecessor_matrix, const DVertexIndexMap& dvertexindex, BrokenLabelMap& broken_label,
	   FrustratedLabelMap& frustrated_label, unsigned int num_frustrated,
	   float cutoff, int neighbors, int frustrated)
{
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Dual_Graph>::edge_descriptor d_edge_descriptor;
  typedef typename graph_traits<Dual_Graph>::edge_iterator d_edge_iterator;
  typedef typename graph_traits<Dual_Graph>::out_edge_iterator d_out_edge_iterator;
  typedef typename graph_traits<Dual_Graph>::vertex_iterator d_vertex_iterator;
  typedef typename property_traits<ContentMap>::value_type Vertex_Content;

  //const unsigned int multiplier = 999999;
  const unsigned int multiplier = 89999999;
  //const unsigned int multiplier = 899999999;
  // shift to ensure that edge weights are positive for Dijkstra
  
  vertex_iterator vi, vend;

  // FIXME: iteration is now over in- and out-egdes
  // which involves twice the necessary computations

  edge_iterator ei, eend;
  for(tie(ei, eend) = edges(g); ei != eend; ++ei)
    if(bondcontent[*ei] != 0) dual_bondcontent[edge_to_dual[*ei]] =  fabs(bondcontent[*ei]);
  
  // prepare to compute shortest paths in Dual_Graph to
  // set up the matching graph
    
  d_vertex_iterator dvi, dvend;

  // compute shortest paths between frustrated plaquettes using
  // Dijkstra's algorithm

  //unsigned int count = 0;    
  for(unsigned int d = 0; d < num_frustrated; ++d) {
    my_dijkstra_shortest_paths(dg, frust_to_plaq[d],
			       make_iterator_property_map(&predecessor_matrix[d][0], dvertexindex),
			       make_iterator_property_map(&distance_matrix[d][0], dvertexindex),
			       get(edge_content, dg), frustrated_label, cutoff, neighbors, frustrated);
    //for(unsigned int i = 0; i < num_vertices(dg); ++i) if(source(predecessor_matrix[d][i], dg) == i) ++count;
  }
  //cout<<count/((double)(num_vertices(dg)*num_frustrated))*100<<" % unsused"<<endl;

  // output plaquette graph to hand it down to the blossom algorithm
  // use plain integer indices here to allow for "j < i"

  // find minimum and maximum edge weights for integer re-scaling
    
  float maxweight = numeric_limits<float>::min(), minweight = numeric_limits<float>::max();
  unsigned true_bonds = 0; 

  for(unsigned int i = 0; i < num_frustrated; ++i) {
    //cout<<frust_to_plaq[i]<<": ";
    for(unsigned int j = 0; j < i; ++j) {
      distance_matrix[i][j] = distance_matrix[i][frust_to_plaq[j]];
      if(distance_matrix[i][j] == MAXFLOAT) continue;
      ++true_bonds;
      maxweight = (distance_matrix[i][frust_to_plaq[j]] > maxweight) ? distance_matrix[i][frust_to_plaq[j]] : maxweight;
      minweight = (distance_matrix[i][frust_to_plaq[j]] < minweight) ? distance_matrix[i][frust_to_plaq[j]] : minweight;
      //cout<<frust_to_plaq[j]<<" ("<<distance_matrix[i][j]<<") ";
    }
    //cout<<endl;
  }
  //cout<<endl;

#ifndef NDEBUG    
  cout<<"using "<<true_bonds<<" bonds for matching"<<endl;
  cout<<"minimum weight: "<<minweight<<", maximum weight: "<<maxweight<<endl;
#endif    

  float maxdiff = 0;

  // maxweight = minweight can essentially only occur for extremely small lattices,
  // but just to be sure ...

  float divider = (maxweight-minweight) ? maxweight-minweight : 1;
    
  // re-scale the edge weights to an integer regime
    
  for(unsigned int i = 0; i < num_frustrated; ++i)
    for(unsigned int j = 0; j < i; ++j) {
      if(distance_matrix[i][j] == MAXFLOAT) continue;
      int discrete = static_cast<int>((distance_matrix[i][j]-minweight)/divider*multiplier);
      float diff = (distance_matrix[i][j]-(((float)(discrete))/multiplier*divider+minweight))/distance_matrix[i][j];
      maxdiff = (diff > maxdiff) ? diff : maxdiff;
      //cout<<"original: "<<distance_matrix[i][j]<<", discrete: "<<discrete<<endl;
      distance_matrix[i][j] = discrete;
    }
    
#ifndef NDEBUG
  cout<<"maximum relative deviation: "<<maxdiff<<endl;
#endif    

  // do the actual matching using the Blossom IV implementation

  perfect_match(num_frustrated, true_bonds, distance_matrix, 0, 0, 0, 1);

#ifndef NDEBUG
  cout<<"done matching"<<endl;
#endif 

  // using the perfect matching, label all the edges of the
  // plaquette graph as broken or satisfied
  
  d_edge_iterator dei, deend;
    
  for(tie(dei, deend) = edges(dg); dei != deend; ++dei) broken_label[*dei] = false;

  unsigned int matchcount = 0, brokenedges_count = 0;

  /*
  cout<<"distance matrix:"<<endl;
  for(unsigned int i = 0; i < num_frustrated; ++i) {
    for(unsigned int j = 0; j < i; ++j) {
      cout<<setw(10)<<-distance_matrix[i][j]*divider/multiplier+minweight;
    }
    cout<<endl;
  }
  cout<<endl;
  */  

  for(unsigned int i = 0; i < num_frustrated; ++i) {
    for(unsigned int j = 0; j < i; ++j) {
      float len = distance_matrix[i][j];
      if(len < 0) { // i and j are matched
	//cout<<frust_to_plaq[i]<<" -> "<<frust_to_plaq[j]<<": "<<-len*divider/multiplier+minweight;
	//cout<<frust_to_plaq[i]<<" -> "<<frust_to_plaq[j]<<": "<<-len;
	++matchcount;
	unsigned int src = i, trg = j;
	unsigned int cur = frust_to_plaq[trg], pred;
	unsigned int pathlength = 0;
	// walk up the predecessor tree to the source
	//cout<<", path: ";
	while((pred = predecessor_matrix[src][cur]) != cur) {
	  //cout<<cur;
	  ++brokenedges_count;
	  ++pathlength;
	  // caution: we can have parallel edges!
	  d_edge_descriptor e = smallest_weight_edge(cur, pred, dual_bondcontent, dg);
	  //cout<<" ("<<dual_bondcontent[e]<<") - ";
	  broken_label[e] = true;
	  // get correct reverse edge by mapping to the original graph and back
	  broken_label[edge_to_dual[reverse_edge(dual_to_edge[e], g)]] = true;
	  cur = pred;
	}
	if(cur != (unsigned int)frust_to_plaq[src]) my_throw("oops, broken predecessor tree\n");
	if(!pathlength) my_throw("oops, path of zero length encountered\n");
	//cout<<cur<<endl;
      }
    }
  }

#ifndef NDEBUG
  cout<<"done setting broken labels"<<endl;
#endif 

  d_out_edge_iterator doi, doend;

  for(tie(dvi, dvend) = vertices(dg); dvi != dvend; ++dvi) {
    unsigned int count_frust = 0;
    for(tie(doi, doend) = out_edges(*dvi, dg); doi != doend; ++doi)
      //count_frust += (broken_label[*doi] and dual_bondcontent[*doi] != 0);
      count_frust += (broken_label[*doi]);
    if(count_frust%2 != frustrated_label[*dvi]) {
#ifndef QUIET
      cerr<<"#### oops, inconsistency due to crossing matching paths"<<endl;
      cerr<<"#### label is "<<count_frust%2<<", should be "<<frustrated_label[*dvi]<<endl;
      cerr<<"plaquette "<<*dvi<<": ";
      for(tie(doi, doend) = out_edges(*dvi, dg); doi != doend; ++doi) {
	cerr<<bondcontent[dual_to_edge[*doi]]<<" ";
      }
      cerr<<endl;
#endif
      return false;
    }
  }
    
#ifndef NDEBUG
  cout<<"found "<<brokenedges_count<<" broken edges"<<endl;
  cout<<"found "<<matchcount<<" matching edges"<<endl;
#endif

  // construct the new spin configuration according to
  // the "broken" or "satisfied" labels, using a
  // breadth-first walk through the graph

  //nonzero_bonds<typename property_map<Graph, edge_content_t>::type> nonzero(get(edge_content, g));
  //nonzero_bonds<typename property_map<Graph, edge_content_t>::type> nonzero(bondcontent);
  //filtered_graph<Graph, nonzero_bonds<typename property_map<Graph, edge_content_t>::type> > fg(g, nonzero);

  //for(tie(dvi, dvend) = vertices(dg); dvi != dvend; ++dvi)
  short int r = 1;
  //breadth_first_search(fg, 0, visitor(reconstruct(edge_to_dual, vertexcontent, bondcontent, broken_label, r)));
  breadth_first_search(g, 0, visitor(reconstruct(edge_to_dual, vertexcontent, bondcontent, broken_label, r)));
    
  return reconstruct_visitor<const EdgeToDualMap, ContentMap, const BondContentMap, BrokenLabelMap, short int>::success;
}

//==============================================================================
//! Wrapper function for the "local" minimization procedures

template<typename Graph, typename Dual_Graph, typename ContentMap, typename BondContentMap, typename DualBondContentMap,
	 typename EdgeToDualMap, typename DualToEdgeMap, typename FrustToPlaqMap, typename PlaqToFrustMap,
	 typename DistanceMatrix, typename PredecessorMatrix, typename DVertexIndexMap, typename BrokenLabelMap,
	 typename FrustratedLabelMap>
bool perform_minimization(const Graph& g, const Dual_Graph& dg, ContentMap& vertexcontent, const BondContentMap& bondcontent,
			  DualBondContentMap& dual_bondcontent, const EdgeToDualMap& edge_to_dual, const DualToEdgeMap& dual_to_edge,
			  const FrustToPlaqMap& frust_to_plaq, const PlaqToFrustMap& plaq_to_frust, DistanceMatrix& distance_matrix,
			  PredecessorMatrix& predecessor_matrix, const DVertexIndexMap& dvertexindex, BrokenLabelMap& broken_label,
			  FrustratedLabelMap& frustrated_label, unsigned int num_frustrated, int neighbors,
			  int frustrated, float cutoff)
{
  bool success = match(g, dg, vertexcontent, bondcontent, dual_bondcontent, edge_to_dual, dual_to_edge, frust_to_plaq,
		       plaq_to_frust, distance_matrix, predecessor_matrix, dvertexindex, broken_label, frustrated_label,
		       num_frustrated, cutoff, neighbors, frustrated);

  return success;
}

//==============================================================================

template<typename Graph> void read_graph(Graph& g, const mystring& file)
{
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename property_map<Graph, edge_content_t>::type bondcontent_map_t;
  
  bondcontent_map_t bondcontent_map = get(edge_content, g);
  
  ifstream in;
  in.open(file.c_str(),ios_in);
  
  unsigned int vol, source, target;

  in>>vol;
  for(unsigned int i = 0; i < vol; ++i) add_vertex(g);

  Bond_Content bc;
  while(in>>source>>target>>bc) {
    edge_descriptor e;
    bool accept;
    tie(e,accept) = edge(vertex(source-1, g), vertex(target-1, g), g);
    if(accept) my_throw("read_graph(): multiple edge detected\n");
    tie(e,accept) = add_edge(vertex(source-1, g), vertex(target-1, g), g);
    if(!accept) my_throw("read_graph(): could not insert edge\n");
    bondcontent_map(e) = bc;
  }

  in.close();
}

//==============================================================================

template<typename Graph> void read_graph_leda(Graph& g, const mystring& file)
{
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename property_map<Graph, edge_content_t>::type bondcontent_map_t;
  
  bondcontent_map_t bondcontent_map = get(edge_content, g);
  
  gzifstream in;
  in.open(file.c_str(),ios_in);
  
  unsigned int vol, source, target, dummy, reverse;

  string sn, sl;
  
  in.ignore(50,'\n');
  in>>sn>>sl;
  in>>vol;
  for(unsigned int i = 0; i < vol; ++i) add_vertex(g);
  for(unsigned int k = 0; k < vol; ++k) { in.ignore(10,'{'); in.ignore(10,'|'); }
  in>>dummy;
  
  Bond_Content bc;
  while(!in.eof()) {
    in>>source>>target>>reverse;
    if(!in.eof()) {
      in.ignore(10,'{');
      in>>bc;
      in.ignore(10,'|');
      edge_descriptor e;
      bool accept;
      tie(e,accept) = edge(vertex(source-1, g), vertex(target-1, g), g);
      if(accept) my_throw("read_graph(): multiple edge detected\n");
      tie(e,accept) = add_edge(vertex(source-1, g), vertex(target-1, g), g);
      if(!accept) my_throw("read_graph(): could not insert edge\n");
      bondcontent_map(e) = bc;
    }
  }

  in.close();
}

//==============================================================================

int main(int argc, char *argv[])
{
  parse_options(&argc,argv);
  struct stat buf;
  if(stat((mystring(params.outprefix)<<".params").c_str(),&buf) == 0) exit(1);
  
  typedef typename adjacency_list_traits<vecS, vecS, directedS>::edge_descriptor edgedesc;
  typedef adjacency_list<vecS, vecS, directedS, no_property,
			 property<edge_content_t, Bond_Content, property<edge_dual_t, edgedesc> > > Graph;

  Graph g;
  read_graph(g, params.graphpath);
  //bool is_planar = boyer_myrvold_planarity_test(g);
  //if(!is_planar) my_throw("input graph not planar!\n");
  unsigned int vol = num_vertices(g);
  Vertex_Content *content_array = new Vertex_Content[vol];
  for(unsigned int k = 0; k < vol; ++k) content_array[k] = 1;
  
  // define dual lattice
  typedef adjacency_list<vecS, vecS, directedS,
    property<vertex_frustrated_t, bool>,
    property<edge_content_t, double,
    property<edge_dual_t, edgedesc,
    property<edge_broken_t, bool> > > > Dual_Graph;

  Dual_Graph dg;

  // define access to graph content
  typedef property_map<Graph, vertex_index_t>::type vertex_index_map_t;
  typedef property_map<Dual_Graph, vertex_index_t>::type dual_vertexindex_t;
  typedef iterator_property_map<Vertex_Content*, vertex_index_map_t> vertexcontent_t;
  typedef property_map<Graph, edge_content_t>::type bondcontent_t;
  typedef property_map<Dual_Graph, edge_content_t>::type dual_bondcontent_t;
  typedef property_map<Graph, edge_dual_t>::type edge_to_dual_t;
  typedef property_map<Dual_Graph, edge_dual_t>::type dual_to_edge_t;
  typedef property_map<Dual_Graph, vertex_frustrated_t>::type frustrated_t;
  typedef property_map<Dual_Graph, edge_broken_t>::type broken_t;
  typedef graph_traits<Dual_Graph>::vertex_descriptor d_vertex_descriptor;
  typedef graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef graph_traits<Graph>::vertex_descriptor vertex_descriptor;
 
  dual_vertexindex_t dvertexindex = get(vertex_index, dg);

  // property maps for graph content
  vertex_index_map_t vertexindex = get(vertex_index, g);
  vertexcontent_t vertexcontent = make_iterator_property_map(&content_array[0], vertexindex);
  bondcontent_t bondcontent = get(edge_content, g);
  dual_bondcontent_t dual_bondcontent = get(edge_content, dg);
  edge_to_dual_t edge_to_dual = get(edge_dual, g);
  dual_to_edge_t dual_to_edge = get(edge_dual, dg);
  frustrated_t frustrated_label = get(vertex_frustrated, dg);
  broken_t broken_label = get(edge_broken, dg);

  // remove zero coupling edges
  zero_bonds<bondcontent_t> zero(get(edge_content, g));
  remove_edge_if(zero, g);
  
  // allocate memory for the matching

  unsigned int num_frustrated = dual_graph(g, dg);

  /*
  ofstream vout;
  vout.open("graph.dot",ios::out);
  write_graphviz(vout, dg);
  vout.close();
  */ 

  unsigned int dvol = num_vertices(dg);
  cout<<vol<<" sites, "<<dvol<<" plaquettes, "<<num_frustrated<<" frustrated plaquettes"<<endl;
  
  unsigned int fulldims[2] = {dvol, dvol};
    
  unsigned int **predecessor_matrix = cc_alloc<unsigned int, 2>(fulldims);
  float **distance_matrix = cc_alloc<float, 2>(fulldims);
  
  int* frust_to_plaq = cc_alloc<int, 1>(fulldims);
  int* plaq_to_frust = cc_alloc<int, 1>(fulldims);
  
  //r250 rannum;
  //rannum.setup(params.seed);
  mersenne rannum;
  unsigned long init_key[4]={0x123, 0x231, 0x235, 0x442}, key_length=4;
  rannum.setup(init_key, key_length);

  for(unsigned int i = 0; i < dvol; ++i) plaq_to_frust[i] = -1;

  unsigned int index = 0;
  for(unsigned int i = 0; i < num_vertices(dg); ++i)
    if(frustrated_label[i]) {
      frust_to_plaq[index] = i;
      plaq_to_frust[i] = index;
      ++index;
    }

  /*
  double en = 0;
  BGL_FORALL_EDGES(e, g, Graph)
    {
      en += bondcontent[e]*vertexcontent[source(e,g)]*vertexcontent[target(e,g)];
    }
  en /= 2;

  cout<<"energy before: "<<en<<endl;
  */  

  perform_minimization(g, dg, vertexcontent, bondcontent, dual_bondcontent, edge_to_dual, dual_to_edge, frust_to_plaq,
		       plaq_to_frust, distance_matrix, predecessor_matrix, dvertexindex, broken_label, frustrated_label,
		       num_frustrated, params.neighbors, params.frustrated, params.min_dijkstra_cutoff);
  double en = 0;
  BGL_FORALL_EDGES(e, g, Graph)
    {
      en += bondcontent[e]*vertexcontent[source(e,g)]*vertexcontent[target(e,g)];
    }
  en /= 2;

  ofstream out;
  out.open(params.outprefix.c_str(),ios_out);
  for(unsigned int i = 0; i < vol; ++i) out<<vertexcontent[i]<<endl;
  cout<<"ground-state energy: "<<-en<<endl;
  out.close();
  
  //cout<<"energy after: "<<en<<endl;
  
  // de-allocate memory
  cc_free<float, 2>(distance_matrix, fulldims);
  cc_free<unsigned int, 2>(predecessor_matrix, fulldims);
  cc_free<int, 1>(frust_to_plaq, fulldims);
  cc_free<int, 1>(plaq_to_frust, fulldims);
  
  return EXIT_SUCCESS;
}

//==============================================================================
