#ifndef BFS_INFO_H
#define BFS_INFO_H

#include <boost/graph/breadth_first_search.hpp>

#include "utilities.h"

using namespace mcrl2::log; // Logging.

namespace detail
{

template<typename BackEdgeMap, typename DistanceMap>
class back_edge_recorder: public boost::default_bfs_visitor
{
protected:
  size_t& m_back_edges;
  BackEdgeMap& m_back_level_edges_per_length;
  const DistanceMap& m_distance;

public:
  back_edge_recorder(size_t& c, BackEdgeMap& m, const DistanceMap& d)
    : m_back_edges(c), m_back_level_edges_per_length(m), m_distance(d)
  {}

  typedef boost::on_non_tree_edge event_filter;

  template<typename Edge, typename Graph>
  void operator()(Edge e, const Graph& g) {
    typename boost::graph_traits<Graph>::vertex_descriptor
      u = source(e, g), v = target(e, g);
    if(m_distance[u] >= m_distance[v])
    {
      size_t delta = m_distance[u] - m_distance[v];
      ++m_back_edges;
      increment(m_back_level_edges_per_length, delta);
    }
  }
};

template< typename BackEdgeMap, typename DistanceMap >
back_edge_recorder<BackEdgeMap, DistanceMap>
record_back_edges(size_t& n, BackEdgeMap& m, const DistanceMap& d)
{
  return back_edge_recorder<BackEdgeMap, DistanceMap>(n, m, d);
}

// Use the convention that changing the color of a vertex to gray
// means adding to the queue, changing color to black is removing from the queue
// the queue size at position k is the size of the queue once we finish
// processing the k'th vertex.
template<class Tag>
struct queue_append_recorder: public boost::default_bfs_visitor
{
  size_t& cur_queue_size;
  typedef Tag event_filter;

  queue_append_recorder(size_t& cur_queue_size)
    : cur_queue_size(cur_queue_size)
  {}

  template<typename Vertex, typename Graph>
  inline void operator()(Vertex, Graph&)
  {
    ++cur_queue_size;
  }
};

template <typename OutIter>
struct queue_remove_recorder: public boost::default_bfs_visitor
{
  size_t& m_cur_queue_size;
  OutIter m_queue_size_out;
  typedef boost::on_finish_vertex event_filter;

  queue_remove_recorder(size_t& cur_queue_size, OutIter queue_size_out)
    : m_cur_queue_size(cur_queue_size),
      m_queue_size_out(queue_size_out)
  {}

  template<typename Vertex, typename Graph>
  inline void operator()(Vertex, Graph&)
  {
    --m_cur_queue_size;
    *(m_queue_size_out)++ = m_cur_queue_size;
  }
};

struct queue_max_remove_recorder: public boost::default_bfs_visitor
{
  size_t& m_cur_queue_size;
  size_t& m_max_queue_size;
  typedef boost::on_finish_vertex event_filter;

  queue_max_remove_recorder(size_t& cur_queue_size, size_t& max)
    : m_cur_queue_size(cur_queue_size),
      m_max_queue_size(max)
  {}

  template<typename Vertex, typename Graph>
  inline void operator()(Vertex, Graph&)
  {
    --m_cur_queue_size;
    m_max_queue_size = std::max(m_max_queue_size, m_cur_queue_size);
  }
};

template<class Tag>
inline queue_append_recorder<Tag>
record_queue_append(size_t& cur_queue_size, Tag)
{
  return queue_append_recorder<Tag>(cur_queue_size);
}

template<typename QueueSizeMap>
inline queue_remove_recorder<QueueSizeMap>
record_queue_remove(size_t& cur_queue_size, QueueSizeMap queue_sizes)
{
  return queue_remove_recorder<QueueSizeMap>(cur_queue_size, queue_sizes);
}

inline queue_max_remove_recorder
record_queue_max_remove(size_t& cur_queue_size, size_t& max)
{
  return queue_max_remove_recorder(cur_queue_size, max);
}

} // namespace detail

template<typename DistanceMap>
inline
std::map<size_t, size_t> bfs_nodes_per_level(const DistanceMap& d)
{
  mCRL2log(verbose) << "Determining number of BFS nodes per level" << std::endl;
  std::map<size_t, size_t> l;
  for(auto i: d)
    increment(l, i);
  return l;
}

template<typename Graph, typename DistanceMap>
inline
typename boost::graph_traits<Graph>::vertices_size_type
bfs_levels(const Graph& g, typename Graph::vertex_descriptor v, DistanceMap& d)
{
  assert(d.size() == num_vertices(g));
  boost::breadth_first_search(g, v,
          boost::visitor(boost::make_bfs_visitor(
            boost::record_distances(&d[0], boost::on_tree_edge()))
      ));
  return *std::max_element(d.begin(), d.end()) + 1;
}

template<typename Graph>
inline
typename boost::graph_traits<Graph>::vertices_size_type
bfs_levels(const Graph& g, typename Graph::vertex_descriptor v = 0)
{
  mCRL2log(verbose) << "Determining number of bfs levels" << std::endl;
  typedef typename boost::graph_traits<Graph>::vertices_size_type vertex_size_t;
  std::vector<vertex_size_t> d(boost::num_vertices(g), 0);
  return bfs_levels(g, v, d);
}

template<typename Graph, typename BackEdgeMap>
inline
typename boost::graph_traits<Graph>::vertices_size_type
back_level_edges(const Graph& g, typename Graph::vertex_descriptor v, BackEdgeMap& m)
{
  typedef typename boost::graph_traits<Graph>::vertices_size_type vertex_size_t;
  vertex_size_t n = 0;
  std::vector<vertex_size_t> d(boost::num_vertices(g));
  boost::breadth_first_search(g, v,
                              boost::visitor(boost::make_bfs_visitor(std::make_pair(
                                boost::record_distances(&d[0], boost::on_tree_edge()),
                                detail::record_back_edges(n, m, &d[0]))
                          ))
      );
  return n;
}

template<typename Graph>
inline
typename boost::graph_traits<Graph>::vertices_size_type
back_level_edges(const Graph& g, typename Graph::vertex_descriptor v = 0)
{
  mCRL2log(verbose) << "Recording back-level edges" << std::endl;
  typedef typename boost::graph_traits<Graph>::vertices_size_type vertex_size_t;
  std::map<vertex_size_t, vertex_size_t> m;
  return back_level_edges(g, v, m);
}

template<typename Graph>
inline
std::vector<typename boost::graph_traits<Graph>::vertices_size_type>
bfs_queue_sizes(const Graph& g, typename Graph::vertex_descriptor v = 0)
{
  mCRL2log(verbose) << "Recording BFS queue sizes" << std::endl;
  typedef typename boost::graph_traits<Graph>::vertices_size_type vertex_size_t;
  vertex_size_t cur = 0;
  std::vector<vertex_size_t> queue_size_map(boost::num_vertices(g), 0);
  boost::breadth_first_search(g, v,
      boost::visitor(boost::make_bfs_visitor(std::make_pair(
          detail::record_queue_append(cur, boost::on_discover_vertex()),
          detail::record_queue_remove(cur, &queue_size_map[0])
      ))));
  return queue_size_map;
}

template<typename Graph>
inline
typename boost::graph_traits<Graph>::vertices_size_type
bfs_max_queue_size(const Graph& g, typename Graph::vertex_descriptor v = 0)
{
  mCRL2log(verbose) << "Recording BFS queue sizes" << std::endl;
  typedef typename boost::graph_traits<Graph>::vertices_size_type vertex_size_t;
  vertex_size_t cur = 0;
  vertex_size_t max = 0;
  std::vector<vertex_size_t> queue_size_map(boost::num_vertices(g), 0);
  boost::breadth_first_search(g, v,
      boost::visitor(boost::make_bfs_visitor(std::make_pair(
          detail::record_queue_append(cur, boost::on_discover_vertex()),
          detail::record_queue_max_remove(cur, max)
      ))));
  return max;
}

#endif // BFS_INFO_H