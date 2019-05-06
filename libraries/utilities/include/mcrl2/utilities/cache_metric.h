// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_UTILITIES_CACHE_METRIC_H
#define MCRL2_UTILITIES_CACHE_METRIC_H

#include <cstddef>
#include <string>
#include <sstream>

namespace mcrl2
{
namespace utilities
{

/// \brief A helper class to keep track of the number of hits and misses for cache-like data structures.
class cache_metric
{
public:
  /// \brief Should be called when searching the cache was a hit.
  void hit() { ++m_hit_count; }

  /// \brief Should be called when searching the cache was a miss.
  void miss() { ++m_miss_count; }

  /// \brief Resets the cache counters.
  void reset()
  {
    m_hit_count = 0;
    m_miss_count = 0;
  }

  /// \returns A message stating x hits, y misses (z %), where x,y,z indicate the number of hits, misses and percentage respectively.
  std::string message() const
  {
    std::stringstream str;
    std::size_t total_count = m_hit_count + m_miss_count;
    str << m_hit_count << " cache hits out of " << total_count << " calls (" << static_cast<double>(m_hit_count) / static_cast<double>(m_miss_count) * 100 << " %)";
    return str.str();
  }

private:
  std::size_t m_hit_count = 0;
  std::size_t m_miss_count = 0;

};

} // namespace utilities
} // namespace mcrl2


#endif // MCRL2_UTILITIES_CACHE_METRIC_H
