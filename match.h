// -*- C++ -*-
// $Header: /home/weigel/cvs/progs/simulation/matching/match.h,v 1.25 2011-08-23 12:29:30 weigel Exp $

/*! \file
   \brief Ground-state calculations for planar XY
   models via matching, definitions and helper functions
*/

#ifndef _MATCH_H
#define _MATCH_H

#include "machine_specific.h"

#include <vector>
#include <unistd.h>
#include <iostream>
#include <map>
#include <limits>
#include <algorithm>
#include <cmath>
#include <boost/graph/breadth_first_search.hpp>

#include "my_throw.h"
#include "bgl_properties.h"

#define NDEBUG
//#define QUIET

extern "C" {
  int perfect_match (int ncount, int ecount, float **dists,
		     int just_fractional, int no_fractional, int use_all_trees,
		     int partialprice);
} 

//! new property for the frustration label of the dual graph vertices
enum vertex_frustrated_t { vertex_frustrated };
namespace boost {
  BOOST_INSTALL_PROPERTY(vertex, frustrated);
}

//! new property for the broken label of the dual graph edges
enum edge_broken_t { edge_broken };
namespace boost {
  BOOST_INSTALL_PROPERTY(edge, broken);
}

using namespace boost;
using namespace std;

//==============================================================================
// struct of simulation parameters

struct _params {
  int siblings;
  int iterations;
  int q;
  int neighbors;
  int frustrated;
  int local;
  int total;
  double max_exc_cutoff;
  double min_exc_cutoff;
  double max_dijkstra_cutoff;
  double min_dijkstra_cutoff;
  double crossfactor;
  double mutation_rate;
  double J_strong;
  double J_strong_threshold;
  double maxover;
  double randomize_start;
  double randomize_factor;
  double noise;
  string graphpath;
  string outprefix;
  string reference_conf;
  string loadstate;
  bool all;
  bool spin_DW;
  bool chiral_DW;
  bool null_DW;
  bool random_DW;
  unsigned int picky;
  int gradient;
  int seed;
  int target_frustrated;
};

extern _params params;
void parse_options(int *argc, char *argv[]);

istream& operator>>(istream& s, _params& p);
ostream& operator<<(ostream& s, _params& p);

//==============================================================================
//! Cyclic successor edge
/*! Return the edge_descriptor of the cyclic successor of an edge.
  \param e edge_descriptor of the considered edge
  \param g the graph to be considered
  \warning pretty slow due to bad hacking style
  \warning does not work for graphs with parallel edges
  \todo Write a new BGL container for adjacency_list with a cyclic
  edge list and corresponding front-end functions for cyclic edge
  traversal.
*/

template<class Graph> typename graph_traits<Graph>::edge_descriptor
edge_succ(const typename graph_traits<Graph>::edge_descriptor& e, const Graph& g)
{
  typename graph_traits<Graph>::out_edge_iterator beg, end, cur;

  for(tie(cur, end) = out_edges(source(e,g), g); *cur != e; ++cur);
  tie(beg, end) = out_edges(source(e,g), g);
  if(++cur != end) return *cur;
  else return *beg;
} 

//==============================================================================
//! Reverse an edge
/*
  \param e edge descriptor for the edge to be reversed
  \param g the graph to be considered
  \warning only makes sense for a bidirectional graph
  \warning does not work for graphs with parallel edges
 */

template<class Graph> typename graph_traits<Graph>::edge_descriptor
reverse_edge(const typename graph_traits<Graph>::edge_descriptor& e, const Graph& g)
{
  //if(edge(target(e,g), source(e,g), g).second == false) {
  //  cout<<"link "<<source(e,g)<<" <-> "<<target(e,g)<<endl;
  //  my_throw("corrupt link reversal information\n");
  //}
  return edge(target(e,g), source(e,g), g).first;
}

//==============================================================================

template<class Descr>
struct _my_less : public binary_function<Descr,Descr,bool>
{
  bool operator()(const Descr& a, const Descr& b) const {
    return (a.m_source < b.m_source || (!(b.m_source < a.m_source) && a.m_target < b.m_target));
  }
};

//==============================================================================
//! Helper type for loop info

template<typename Descriptor> struct _edge_info_t {
  Descriptor dual;
  unsigned char label;
};

//==============================================================================
//! Compute the dual graph
/*
  \param g reference to the original graph
  \param dg reference to the dualized graph
 */

template<class Graph, class Dual_Graph> unsigned int dual_graph(Graph& g, Dual_Graph& dg)
{
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
  typedef typename graph_traits<Dual_Graph>::vertex_descriptor d_vertex_descriptor;
  typedef typename graph_traits<Dual_Graph>::edge_descriptor d_edge_descriptor;
  typedef typename graph_traits<Dual_Graph>::vertex_iterator d_vertex_iterator;

  // property maps: edge_content and edge to dual edge mappings
  typedef property_map<Graph, edge_content_t> bondcontent_t;
  typedef property_map<Dual_Graph, edge_content_t> dual_bondcontent_t;
  typedef property_map<Graph, edge_dual_t> edge_to_dual_t;
  typedef property_map<Dual_Graph, edge_dual_t> to_dual_edge_t;
  typedef property_map<Dual_Graph, vertex_frustrated_t> frustrated_t;

  typename bondcontent_t::type bondcontent = get(edge_content, g);
  typename dual_bondcontent_t::type dual_bondcontent = get(edge_content, dg);
  typename edge_to_dual_t::type edge_to_dual = get(edge_dual, g);
  typename to_dual_edge_t::type dual_to_edge = get(edge_dual, dg);
  typename frustrated_t::type frustrated = get(vertex_frustrated, dg);

  // maps for loop construction and bond setting
  map<edge_descriptor, _edge_info_t<d_vertex_descriptor>, _my_less<edge_descriptor>,
    allocator<pair<const edge_descriptor, _edge_info_t<d_vertex_descriptor> > > > loop;
  //map<edge_descriptor, _edge_info_t<d_vertex_descriptor>, _my_less<edge_descriptor> > loop;

  edge_iterator ei, eend;

  // walk through graph and label edges according to
  // loop "numbers", insert vertices in dual graph
  for(tie(ei, eend) = edges(g); ei != eend; ++ei) {
    // ignore open boundary edges
    if(bondcontent[*ei] == 0) loop[*ei].label = 2;
    if(loop[*ei].label == 0) { // new loop or dual vertex
      d_vertex_descriptor u = add_vertex(dg);
      edge_descriptor pos = reverse_edge(*ei, g);
      do {
	// mark bond
	loop[reverse_edge(pos, g)].dual = u;
	loop[reverse_edge(pos, g)].label = 1;
	// walk backwards along the loop
	// open boundary edges are ignored
      	do { pos = edge_succ(pos, g); } while(bondcontent[pos] == 0);
	//pos = edge_succ(pos, g);
	pos = reverse_edge(pos, g);
      } while(pos != reverse_edge(*ei, g));
    }
  }

  unsigned int count = 0, frust = 0;
  d_edge_descriptor de;
  bool accept;

  unsigned count_neg = 0;
  // walk through the graph a second time to set the dual bonds
  // such that they can be cyclically ordered
  for(tie(ei, eend) = edges(g); ei != eend; ++ei) {
    if(loop[*ei].label < 2) { // not yet considered loop or dual vertex
      edge_descriptor pos = reverse_edge(*ei, g);
      frust = 0;
      do {
	edge_descriptor orig = reverse_edge(pos, g);
	// insert new edge
	tie(de,accept) = add_edge(loop[orig].dual, loop[pos].dual, dg);
	// set up edge to dual edge (and reverse) mapping
	edge_to_dual[orig] = de; dual_to_edge[de] = orig;
	// copy edge content
	dual_bondcontent[de] = bondcontent[orig];
	// mark the edge as done
	loop[reverse_edge(pos, g)].label = 2;
	// count number of negative bonds
	frust += (bondcontent[orig] < 0);
        count_neg += (bondcontent[orig] < 0);
	// walk backwards along the loop
	// open boundary edges are ignored
	do { pos = edge_succ(pos, g); } while(bondcontent[pos] == 0);
	//pos = edge_succ(pos, g);
	pos = reverse_edge(pos, g);
      } while(pos != reverse_edge(*ei, g));
      // frustrated plaquettes have an odd number of negative bonds
      if(frust%2) {
	frustrated[loop[*ei].dual] = true;
	++count;
      }
      else frustrated[loop[*ei].dual] = false;
    }
  }

  // for open boundary conditions:
  // assume that the first plaquette found is the "outer" plaquette
  // and define the frustration of the outer vertex such that the
  // total number of frustrated plaquettes is even
  // (for periodic bc, count is always even)

  if(count%2) {
    d_vertex_iterator vbegin, vend;
    tie(vbegin, vend) = vertices(dg);
    frustrated[*vbegin] = !frustrated[*vbegin];
    cout<<"changed frustration status of outer plaquette"<<endl;
  }

#ifndef QUIET
  cout<<"found "<<count<<" frustrated plaquettes"<<endl;
#endif
  //cout<<count_neg/2<<"\t"<<count;
  return count;
}

//==============================================================================
//! Helper class for filtered graph with non-zero edge weights

template <typename EdgeWeightMap>
struct nonzero_bonds {
  nonzero_bonds() { }
  nonzero_bonds(const EdgeWeightMap& weight) : m_weight(weight) { }
  template <typename Edge>
  bool operator()(const Edge& e) const {
    return (get(m_weight, e) != 0);
    //return true;
  }
  EdgeWeightMap m_weight;
};

template <typename EdgeWeightMap>
struct zero_bonds {
  zero_bonds() { }
  zero_bonds(const EdgeWeightMap& weight) : m_weight(weight) { }
  template <typename Edge>
  bool operator()(const Edge& e) const {
    return (get(m_weight, e) == 0);
  }
  EdgeWeightMap m_weight;
};

template <typename EdgeWeightMap>
struct all_bonds {
  all_bonds() { }
  all_bonds(const EdgeWeightMap& weight) { }
  template <typename Edge>
  bool operator()(const Edge& e) const {
    return true;
  }
};

//==============================================================================
//! Breadth-first search visitor used for constructing the spin configuration
/*!
  Flips discovered spins when the corresponding edge broken label
  does not match spin content
*/

template<typename edge_to_dual_t, typename vertexcontent_t, typename bondcontent_t, typename broken_t, typename vector_t>
class reconstruct_visitor : public default_bfs_visitor {
public:
  //! static functions for statistics                                                                                                                                            
  static bool success;
  static unsigned int error_count, run_count;
  static void reset() { error_count = run_count = 0; }
  static double error_rate() { return error_count/((double)run_count); }

  reconstruct_visitor(edge_to_dual_t& etd, vertexcontent_t& vc, bondcontent_t& bc,  broken_t& br, vector_t& rr)
    : edge_to_dual(etd), vertexcontent(vc), bondcontent(bc), broken_label(br), r(rr)
  { success = true; real = false; flip_count = 0; ++run_count; }

  ~reconstruct_visitor() {
#ifndef NDEBUG
    if(flip_count) cout<<"flipped "<<flip_count<<" spins"<<endl;
#endif 
    if(real and !success) ++error_count;
  }

  template<typename Edge, typename Graph> void tree_edge(Edge e, const Graph& g) {
    typename graph_traits<Graph>::vertex_descriptor  u = source(e, g), v = target(e, g);
    //if(bondcontent[e] == 0) return;
    double en = bondcontent[e]*(vertexcontent[u]*r)*(vertexcontent[v]*r);
    bool broken = broken_label[edge_to_dual[e]];
    // spin should be flipped if broken label does not fit vertex content
    if((broken and en > 0) or (!broken and en < 0)) {
      vertexcontent[v] -= 2*(vertexcontent[v]*r)*r;
      ++flip_count;
    }
    real = true;
  }

  template<typename Edge, typename Graph> void non_tree_edge(Edge e, const Graph& g) {
    typename graph_traits<Graph>::vertex_descriptor  u = source(e, g), v = target(e, g);
    //if(bondcontent[e] == 0) return;
    double en = bondcontent[e]*(vertexcontent[u]*r)*(vertexcontent[v]*r);
    bool broken = broken_label[edge_to_dual[e]];
    // broken labels are not satisfiable
    if((broken and en > 0) or (!broken and en < 0)) {
      //if((broken and en > 1e-5) or (!broken and en < -1e-5)) {
      //vertexcontent[v] -= 2*(vertexcontent[v]*r)*r;
      //my_throw("reconstruct_visitor::non_tree_edge(): bad \"broken\" labels\n");
#ifndef QUIET
      if(success) cerr<<"#### reconstruct_visitor::non_tree_edge(): bad \"broken\" labels\n";
#endif
      success = false;
    }
  }

private:
  bool real;
  unsigned int flip_count;
  edge_to_dual_t& edge_to_dual;
  vertexcontent_t& vertexcontent;
  bondcontent_t& bondcontent;
  broken_t& broken_label;
  vector_t& r;
};

template<typename edge_to_dual_t, typename vertexcontent_t, typename bondcontent_t, typename broken_t, typename vector_t>
bool reconstruct_visitor<edge_to_dual_t, vertexcontent_t, bondcontent_t, broken_t, vector_t>::success = true;
template<typename edge_to_dual_t, typename vertexcontent_t, typename bondcontent_t, typename broken_t, typename vector_t>
unsigned int reconstruct_visitor<edge_to_dual_t, vertexcontent_t, bondcontent_t, broken_t, vector_t>::run_count = 0;
template<typename edge_to_dual_t, typename vertexcontent_t, typename bondcontent_t, typename broken_t, typename vector_t>
unsigned int reconstruct_visitor<edge_to_dual_t, vertexcontent_t, bondcontent_t, broken_t, vector_t>::error_count = 0;

//==============================================================================
//! Convenience function for reconstruct_visitor

template<typename edge_to_dual_t, typename vertexcontent_t, typename bondcontent_t, typename broken_t, typename vector_t>
reconstruct_visitor<edge_to_dual_t, vertexcontent_t, bondcontent_t, broken_t, vector_t>
reconstruct(edge_to_dual_t& etd, vertexcontent_t& vc, bondcontent_t& bc, broken_t& br, vector_t& rr)
{
  return reconstruct_visitor<edge_to_dual_t, vertexcontent_t, bondcontent_t, broken_t, vector_t>(etd, vc, bc, br, rr);
}

//==============================================================================
//! Prints dual graph with content

template<typename Dual_Graph, typename bondcontent_t, typename frustrated_t, typename broken_t> void
print_dual_graph(Dual_Graph& dg, bondcontent_t& dual_bondcontent,
		 frustrated_t& frustrated, broken_t& broken_label)
{
  typedef typename graph_traits<Dual_Graph>::vertex_iterator d_vertex_iterator;
  typedef typename graph_traits<Dual_Graph>::edge_iterator d_edge_iterator;
  typedef typename graph_traits<Dual_Graph>::out_edge_iterator d_out_edge_iterator;
  d_vertex_iterator dvi, dvend;
  d_out_edge_iterator dei, deend;

  cout<<"Dual graph ============================="<<endl;
  for(tie(dvi, dvend) = vertices(dg); dvi != dvend; ++dvi) {
    cout<<*dvi<<"["<<frustrated[*dvi]<<"]"<<":";
    for(tie(dei, deend) = out_edges(*dvi, dg); dei != deend; ++dei) {
      cout<<" ";
      cout.width(3);
      cout<<*dei<<"[";
      cout.width(2);
      cout<<dual_bondcontent[*dei]<<",";
      cout<<broken_label[*dei]<<"]";
    }
    cout<<endl;
  }
  cout<<"========================================"<<endl;
}

//==============================================================================
//! Prints original graph with content

template<typename Graph, typename vertexcontent_t, typename bondcontent_t> void
print_orig_graph(const Graph& g, const vertexcontent_t& vertexcontent, const bondcontent_t& bondcontent)
{
  typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iterator;
  vertex_iterator vi, vend;
  out_edge_iterator ei, eend;

  cout<<"Graph =================================="<<endl;
  for(tie(vi, vend) = vertices(g); vi != vend; ++vi) {
    cout<<*vi<<"["<<vertexcontent[*vi]<<"]"<<":";
    for(tie(ei, eend) = out_edges(*vi, g); ei != eend; ++ei) {
      cout<<" ";
      cout.width(3);
      cout<<*ei<<"[";
      cout.width(2);
      cout<<bondcontent[*ei]<<"]";
    }
    cout<<endl;
  }
  cout<<"========================================"<<endl;
}

//==============================================================================
//! Prints original graph with content
/*!
  Format suitable for the spin glass
  server:
  
  http://www.informatik.uni-koeln.de/ls_juenger/research/sgs/sgs.html
*/

template<typename Graph, typename vertexcontent_t, typename bondcontent_t> void
print_sgs_graph(const Graph& g, const vertexcontent_t& vertexcontent, const bondcontent_t& bondcontent)
{
  typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
  edge_iterator ei, eend;

  cout<<"type: gauss"<<endl;
  cout<<"size: "<<(int)sqrt((double)num_vertices(g))<<endl;
  cout<<"name: myname\n"<<endl;

  for(tie(ei, eend) = edges(g); ei != eend; ++ei) {
    unsigned int src = source(*ei, g)+1, trg = target(*ei, g)+1;
    if(src > trg) swap(src, trg);
    cout<<src<<" "<<trg<<" "<<bondcontent[*ei]<<endl;
  }
}

//=============================================================================
//! For vertices with parallel edges, selects edge with smallest weight
/*
  Broken bonds at "corners": consider corner of a square lattice with open bc
  of the form
  
 1 a
  *--*--...
 b|  |
  *--*--...
  .  .
  .  .

  (1) in a ground state, bonds a and b can never be broken at the same time,
      since, if they were, spin 1 could be flipped to satisfy bonds a and b
      without any non-local effects

  (2) thus, either a or b or none of them are broken; if it is a or b, it
      will be the bond with the weaker (modulus of) coupling, since otherwise
      1 could be flipped to satisfy the one and break the other bond, lowering
      the energy

  Hence, for the matching problem the bond dual to the stronger of a and b
  can be completely neglected, it is _always_ unbroken in the ground state.

  \warning not necessary any more
 */

template<typename Graph, typename Edgeweight> typename graph_traits<Graph>::edge_descriptor
smallest_weight_edge(const typename graph_traits<Graph>::vertex_descriptor& u,
		     const typename graph_traits<Graph>::vertex_descriptor& v,
		     const Edgeweight& edgeweight, const Graph& g)
{
  typedef typename graph_traits<Graph>::out_edge_iterator out_edge_iterator;
  typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

  // should be 2 for the square lattice, but could be larger for random lattices
  const unsigned int maxedges = 20;
  unsigned int count = 0, minlabel = 0;
  float minweight = numeric_limits<float>::max();
  
  out_edge_iterator ei, eiend;
  edge_descriptor e;

  edge_descriptor lst[maxedges];

  for(tie(ei, eiend) = out_edges(u, g); ei != eiend; ++ei)
    if(target(*ei, g) == v) {
      lst[count] = *ei;
      float weight = edgeweight[*ei];
      if(weight < minweight) {
	minweight = weight;
	minlabel = count;
      }
      ++count;
    }

  return lst[minlabel];
}

//=============================================================================

template<class Dual_Graph, class DualBondContentMap, class FrustratedLabelMap> bool
update_frustration_label(const Dual_Graph& dg, DualBondContentMap& dual_bondcontent, FrustratedLabelMap& frustrated_label, typename graph_traits<Dual_Graph>::vertex_descriptor v)
{
  typedef typename graph_traits<Dual_Graph>::out_edge_iterator out_edge_iterator;
  
  out_edge_iterator ei, eend;
  int count = 0;
  for(tie(ei, eend) = out_edges(v, dg); ei != eend; ++ei) count += (dual_bondcontent[*ei] < 0);
  frustrated_label[v] = (count % 2) ? 1 : 0;
  return frustrated_label[v];
}

//=============================================================================

#endif
