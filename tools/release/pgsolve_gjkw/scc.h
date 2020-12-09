// Author(s): Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file /path/to/file.ext
/// \brief Description comes here

#ifndef SCC_INFO_H
#define SCC_INFO_H

#include <boost/graph/strong_components.hpp>
#include "bfs.h"
#include "degree.h"

template<typename Graph, typename ComponentMap>
inline
Graph quotient_graph(const Graph& g, ComponentMap m, const size_t num_components)
{
  /* Fix quotient graph such that vertex 0 is in component 0 */
  size_t comp_zero = m[0];
  /* Vector containing size for each component */
  std::vector<size_t> scc_size(num_components, 0);
  for(auto i = m.begin(); i != m.end(); ++i)
  {
    if(*i == 0)
      *i = comp_zero;
    else if(*i == comp_zero)
      *i = 0;
    scc_size[*i]++;
  }

  Graph result;
  for(size_t i = 0; i < num_components; ++i)
    boost::add_vertex(result);

  typename boost::graph_traits< Graph >::edge_iterator i, end;
  for(boost::tie(i, end) = boost::edges(g); i != end; ++i)
  {
    typename Graph::vertex_descriptor u = boost::source(*i,g), v = boost::target(*i,g);
    if(m[u] != m[v] || scc_size[m[u]] >= 1) // Add edge if we have a non-trivial component
      boost::add_edge(m[u], m[v], result);
  }
  return result;
}

template<typename Graph>
inline
typename boost::graph_traits<Graph>::vertices_size_type
quotient_height(const Graph& g)
{
  typedef typename boost::graph_traits<Graph>::vertices_size_type vertex_size_t;
  std::vector<vertex_size_t> components (boost::num_vertices(g), 0);
  vertex_size_t ncomponents = boost::strong_components(g, &components[0]);
  return bfs_levels(quotient_graph(g, components, ncomponents));
}

template<typename Graph>
inline
typename boost::graph_traits<Graph>::vertices_size_type
sccs(const Graph& g)
{
  typedef typename boost::graph_traits<Graph>::vertices_size_type vertex_size_t;
  std::vector<vertex_size_t> components (boost::num_vertices(g), 0);
  return boost::strong_components(g, &components[0]);
}

template<typename Graph>
inline
typename boost::graph_traits<Graph>::vertices_size_type
trivial_sccs(const Graph& g)
{
  typedef typename boost::graph_traits<Graph>::vertices_size_type vertex_size_t;
  std::vector<vertex_size_t> components (boost::num_vertices(g), 0);
  boost::strong_components(g, &components[0]);
  return count_elements_occurring_exactly_n_times(components);
}

template<typename Graph>
inline
typename boost::graph_traits<Graph>::vertices_size_type
terminal_sccs(const Graph& g)
{
  typedef typename boost::graph_traits<Graph>::vertices_size_type vertex_size_t;
  std::vector<vertex_size_t> components (boost::num_vertices(g), 0);
  vertex_size_t ncomponents = boost::strong_components(g, &components[0]);
  return vertices_out_degree_n(quotient_graph(g, components, ncomponents));
}

#endif // SCC_INFO_H