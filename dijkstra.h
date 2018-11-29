// -*- C++ -*-
// $Header: /home/weigel/cvs/progs/simulation/matching/dijkstra.h,v 1.14 2004-12-02 19:20:26 weigel Exp $

/*! \file
   \brief Our implementation of the Dijkstra algorithm for
   shortest paths
*/

#ifndef DIJKSTRA_H
#define DIJKSTRA_H

#include <boost/pending/mutable_queue.hpp>
#include <boost/pending/indirect_cmp.hpp>
#include <boost/graph/properties.hpp>
#include <vector>
#include "match.h"

//==============================================================================
//! Helper function for edge relaxation

template <class Graph, class WeightMap, class PredecessorMap, class DistanceMap>
bool my_relax(typename graph_traits<Graph>::edge_descriptor e, 
	   const Graph& g, const WeightMap& w, 
	   PredecessorMap& p, DistanceMap& d)
{
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  Vertex u = source(e, g), v = target(e, g);
  typedef typename property_traits<DistanceMap>::value_type D;
  typedef typename property_traits<WeightMap>::value_type W;
  D d_u = get(d, u), d_v = get(d, v);
  W w_e = get(w, e);
  
  if(w_e + d_u < d_v) {
    put(d, v, d_u + w_e);
    put(p, v, u);
    return true;
  }
  else return false;
}

//==============================================================================
//! Implementation of Dijsktra's algorithm with distance cut-off and
//! edge predecessor map

template <class Graph, class PredecessorMap, class DistanceMap, class WeightMap, class FrustratedMap>
void my_dijkstra_shortest_paths(const Graph& g, typename graph_traits<Graph>::vertex_descriptor s,
				const PredecessorMap& pm, const DistanceMap& dm,
				const WeightMap& wm, const FrustratedMap& fm, float cutoff, int neighbors, int frustrated)
{
  typedef typename property_map<Graph, vertex_index_t>::type vertexindex_t;
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::edge_descriptor Edge;
  typedef iterator_property_map<default_color_type*, vertexindex_t> ColorMap;
  typedef typename property_traits<ColorMap>::value_type ColorValue;
  typedef color_traits<ColorValue> Color;
  typedef typename property_traits<ColorMap>::value_type ColorValue;  
  typedef typename property_traits<DistanceMap>::value_type D;
  typedef indirect_cmp<DistanceMap, std::less<D> > IndirectCmp;
  typedef mutable_queue<Vertex, std::vector<Vertex>, IndirectCmp, vertexindex_t> MutableQueue;

  vertexindex_t index_map = get(vertex_index, g);
  vector<default_color_type> col(num_vertices(g));
  default_color_type c = white_color;
  ColorMap color = make_iterator_property_map(&col[0], index_map, c);
  IndirectCmp icmp(dm, std::less<D>());
  MutableQueue Q(num_vertices(g), icmp, index_map);

  D m_zero = 0;

  typename graph_traits<Graph>::vertex_iterator ui, ui_end;
  typename graph_traits<Graph>::out_edge_iterator ei, ei_end;

  for (tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui) {
    put(dm, *ui, MAXFLOAT);
    put(pm, *ui, *ui);
    put(color, *ui, white_color);
  }
  put(dm, s, m_zero);

  int seen = 0, seen_frustrated = 0;
  put(color, s, Color::gray());
  Q.push(s);

  while (! Q.empty()) {
    Vertex u = Q.top(); Q.pop();
    for (tie(ei, ei_end) = out_edges(u, g); ei != ei_end; ++ei) {
      Vertex v = target(*ei, g);
      ColorValue v_color = get(color, v);
      if (v_color == Color::white()) { 
	my_relax(*ei, g, wm, pm, dm);
	put(color, v, Color::gray());
	Q.push(v);
      } else {
	if (v_color == Color::gray())
	  if (my_relax(*ei, g, wm, pm, dm)) Q.update(target(*ei, g));
      }
    }
    put(color, u, Color::black());
    ++seen; if(fm[u]) ++seen_frustrated;
    if(seen_frustrated > frustrated and seen > neighbors and dm[u] > cutoff) {
      //cout<<"final distance: "<<dm[u]<<endl;
      break;
    }
  }
  
  while (! Q.empty()) {
    Vertex u = Q.top(); Q.pop();
    dm[u] = MAXFLOAT;
  }
}

//==============================================================================

#endif
