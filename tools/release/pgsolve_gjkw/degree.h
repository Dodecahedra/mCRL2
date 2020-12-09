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

#ifndef DEGREE_H
#define DEGREE_H

#include <limits>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

using namespace mcrl2::log;

template <typename Graph>
inline
double avg_degree(const Graph& g)
{
  mCRL2log(verbose) << "Computing average degree" << std::endl;
  double v = boost::num_vertices(g);
  double e = boost::num_edges(g);
  return e/v;
}

template <typename Graph>
inline
typename boost::graph_traits<Graph>::vertices_size_type
max_degree(const Graph& g)
{
  mCRL2log(verbose) << "Computing maximum degree" << std::endl;
  size_t result = 0;

  typename boost::graph_traits<Graph>::vertex_iterator i, end;
  for (boost::tie(i, end) = vertices(g); i != end; ++i)
  {
    result = std::max(result, boost::degree(*i, g));
  }
  return result;
}

template <typename Graph>
inline
typename boost::graph_traits<Graph>::vertices_size_type
min_degree(const Graph& g)
{
  mCRL2log(verbose) << "Computing minimum degree" << std::endl;
  size_t result = std::numeric_limits<size_t>::max();

  typename boost::graph_traits<Graph>::vertex_iterator i, end;
  for (boost::tie(i, end) = vertices(g); i != end; ++i)
  {
    result = std::min(result, boost::degree(*i, g));
  }
  return result;
}

template <typename Graph>
inline
typename boost::graph_traits<Graph>::vertices_size_type
max_in_degree(const Graph& g)
{
  mCRL2log(verbose) << "Computing maximum in-degree" << std::endl;
  size_t result = 0;

  typename boost::graph_traits<Graph>::vertex_iterator i, end;
  for (boost::tie(i, end) = vertices(g); i != end; ++i)
  {
    result = std::max(result, boost::in_degree(*i, g));
  }
  return result;
}

template <typename Graph>
inline
typename boost::graph_traits<Graph>::vertices_size_type
min_in_degree(const Graph& g)
{
  mCRL2log(verbose) << "Computing minimum in-degree" << std::endl;
  size_t result = std::numeric_limits<size_t>::max();

  typename boost::graph_traits<Graph>::vertex_iterator i, end;
  for (boost::tie(i, end) = vertices(g); i != end; ++i)
  {
    result = std::min(result, boost::in_degree(*i, g));
  }
  return result;
}

template <typename Graph>
inline
double avg_in_degree(const Graph& g)
{
  mCRL2log(verbose) << "Computing average in-degree" << std::endl;
  size_t sum = 0;

  typename boost::graph_traits<Graph>::vertex_iterator i, end;
  for (boost::tie(i, end) = vertices(g); i != end; ++i)
  {
    sum += boost::in_degree(*i, g);
  }
  return static_cast<double>(sum)/static_cast<double>(boost::num_vertices(g));
}

template <typename Graph>
inline
typename boost::graph_traits<Graph>::vertices_size_type
max_out_degree(const Graph& g)
{
  mCRL2log(verbose) << "Computing maximum out-degree" << std::endl;
  size_t result = 0;

  typename boost::graph_traits<Graph>::vertex_iterator i, end;
  for (boost::tie(i, end) = vertices(g); i != end; ++i)
  {
    result = std::max(result, boost::out_degree(*i, g));
  }
  return result;
}

template <typename Graph>
inline
typename boost::graph_traits<Graph>::vertices_size_type
min_out_degree(const Graph& g)
{
  mCRL2log(verbose) << "Computing minimum out-degree" << std::endl;
  size_t result = std::numeric_limits<size_t>::max();

  typename boost::graph_traits<Graph>::vertex_iterator i, end;
  for (boost::tie(i, end) = vertices(g); i != end; ++i)
  {
    result = std::min(result, boost::out_degree(*i, g));
  }
  return result;
}

template <typename Graph>
inline
double avg_out_degree(const Graph& g)
{
  mCRL2log(verbose) << "Computing average out-degree" << std::endl;
  size_t sum = 0;

  typename boost::graph_traits<Graph>::vertex_iterator i, end;
  for (boost::tie(i, end) = vertices(g); i != end; ++i)
  {
    sum += boost::out_degree(*i, g);
  }
  return static_cast<double>(sum)/static_cast<double>(boost::num_vertices(g));
}

template <typename Graph>
inline
typename boost::graph_traits<Graph>::vertices_size_type
vertices_out_degree_n(const Graph& g, const size_t n = 0)
{
  mCRL2log(verbose) << "Determining nr of vertices with outdegree " << n << std::endl;
  typename boost::graph_traits<Graph>::vertices_size_type result = 0;
  typename boost::graph_traits< Graph >::vertex_iterator i, end;
  for(boost::tie(i, end) = boost::vertices(g); i != end; ++i)
  {
    if(boost::out_degree(*i, g) == n)
      ++result;
  }
  return result;
}

#endif // DEGREE_H