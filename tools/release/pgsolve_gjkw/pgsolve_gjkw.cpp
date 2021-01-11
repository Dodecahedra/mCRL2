// Author(s): Koen Degeling
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file pgsolve_gjkw.cpp
/// Reads a parity game from input and reduces it using the Divergence-preserving 
/// branching bisimulation equivalence using the O(m log m) algorithm 
/// [Groote/Jansen/Keiren/Wijs 2017]

#include "liblts_bisim_gjkw_kripke.h"
#include "mcrl2/lts/lts_fsm.h"
#include "mcrl2/utilities/input_output_tool.h"
#include "pg.h"
#include "pgsolver_io.h"
#include "utilities.h"
#include "mcrl2/lts/lts_algorithm.h"
#include "scc.h"

#include <boost/graph/strong_components.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

using namespace mcrl2::lts; // For reduce() function.
using namespace mcrl2::log; // Logging.

namespace mcrl2 {

class pg_convert
{

  protected:
    /** Storing Parity Game */
    parity_game_t m_pg;
    /** Reduced Kripke structure in LTS */
    lts::lts_fsm_t m_lts;

  public:
  pg_convert(parity_game_t pg)
    : m_pg(pg)
    {}

  parity_game_t get_parity_game()
  {
    return m_pg;
  }

  void convert_ks_to_pg()
  /*  Convert reduced Kripke structure in lts_fsm_t lts to a Parity Game in 
      parity_game_t format. */
  {
    // Clear old Parity Game and construct new one from m_lts
    m_pg.clear();
    std::map<size_t, size_t> state_vertex_map;

    // Loop over states and add vertices to graph
    for(size_t i=0; i!=m_lts.num_states(); ++i )
    {
      std::vector<std::size_t> l = m_lts.state_label(i);
      printf("Label of the state: {%u,%u}.\n", l[0], l[1]);
      size_t v = boost::add_vertex(m_pg);
      printf("New vertex added: %u.\n", v);
      // Set properties of the new state
      m_pg[v].player = l[0] == 0 ? even : odd;
      m_pg[v].prio = l[1];
      // Safe the mapping from state to vertex for later when adding edges.
      state_vertex_map[i] = v;
    }

    for(const transition& t: m_lts.get_transitions())
    {
      size_t u = state_vertex_map[t.from()];
      size_t v = state_vertex_map[t.to()];
      boost::add_edge(u, v, m_pg);
    }
  }

  void reduce_pg_scc(std::string file)
  /* Reduce Parity Game graph to remove Strongly Connected Components. */
  {
    typedef typename boost::graph_traits<parity_game_t>::vertices_size_type vertex_size_t;
    parity_game_t* reduced = new parity_game_t(m_pg);
    boost::graph_traits<parity_game_t>::edge_iterator e, m;
    /* We (deep) copied the graph and go over all edges and remove any edge for which
      the source and target label are not the same (to ensure our SCC have same label)*/
    for(boost::tie(e,m) = edges(m_pg); e != m; ++e)
    {
      vertex_size_t s = source(*e, m_pg);
      vertex_size_t t = target(*e, m_pg);
      if(m_pg[s].prio != m_pg[t].prio || m_pg[s].player != m_pg[t].player)
      {
        // Since we aren't removing edges from m_pg, this is safe.
        boost::remove_edge(s,t,*reduced);
      }
    }
    /* We compute the strong components and go over the edges in our original graph.
      We merge all vertices that are in the same SCC. */
    std::vector<vertex_size_t> components (boost::num_vertices(*reduced), 0);
    size_t num = boost::strong_components(*reduced, &components[0]);
    std::map<size_t, vertex_size_t> component_map;
    mCRL2log(verbose) << "Shrinking the graph" << std::endl;
    boost::graph_traits<parity_game_t >::vertex_iterator i, n;
    for(boost::tie(i, n) = vertices(m_pg); i != n; ++i)
    {
      size_t group = components[*i];
      if(component_map.count(group)==1) 
      // If already present, we shrink the graph by merging i into u.
      {
        vertex_size_t u = component_map[group];
        merge(m_pg, u, *i);
        m_pg[*i].unused = true;
      } 
      else
      // We haven't seen any vertices of this group yet, so we add it to the map.
      {
        component_map[group] = *i;
      }
    }
    // We go over all the vertices and remove those which are unused.
    for(auto i = 0; i < boost::num_vertices(m_pg);)
    {
      if(m_pg[i].unused)
      {
        boost::clear_vertex(i, m_pg);
        boost::remove_vertex(i, m_pg);
      } 
      else 
      {
        i++;
      }
    }
  }

  typedef typename boost::graph_traits<parity_game_t>::
      vertices_size_type vertex_size_t;
  void merge(parity_game_t& pg, vertex_size_t s, vertex_size_t r)
  /* Merge vertex r into vertex s*/
  {
    boost::graph_traits<parity_game_t>::out_edge_iterator oe, m;
    for(boost::tie(oe,m) = out_edges(r, pg); oe != m; ++oe)
    {
      boost::add_edge(s, target(*oe, pg), pg); // Add new edge
    }
    boost::graph_traits<parity_game_t>::in_edge_iterator ie, n;
    for(boost::tie(ie,n) = in_edges(r, pg); ie != n; ++ie)
    {
      boost::add_edge(source(*ie, pg), s, pg); // Add new edge
    }
  }
  
  void convert_pg_to_ks()
  /* Convert Parity Game to a Kripke structure contained in lts_fsm_t. */
  {
    // LTS in which we save the Kripke structure
    std::map<size_t, size_t> vertex_state_map;
    std::map<size_t, size_t> mp;
    std::map<size_t, size_t> mp1;
    m_lts.add_process_parameter("Parity", "Nat");
    m_lts.add_process_parameter("Priority", "Nat");
    size_t zero = m_lts.add_state_element_value(0, "0");
    size_t one = m_lts.add_state_element_value(0, "1");
    mp1[0] = zero;
    mp1[1] = one;
    unsigned int max_priority = 8;
    for(unsigned int i = 0; i != max_priority; ++i)
    {
      size_t l = m_lts.add_state_element_value(1, std::to_string(i));
      mp[i] = l;
    }

    boost::graph_traits<parity_game_t >::vertex_iterator i, n;
    // We loop over vertices and add new states to the lts
    for(boost::tie(i, n) = vertices(m_pg); i != n; ++i)
    {
      size_t player = m_pg[*i].player;
      size_t priority = m_pg[*i].prio;
      state_label_fsm label = 
          state_label_fsm(std::vector<std::size_t> {player,priority});
      size_t state = m_lts.add_state(label);
      // Add the (vertex,state) pair to map
      vertex_state_map[*i] = state;
    }
    m_lts.set_initial_state(0);
    size_t label = m_lts.add_action(action_label_string(" "));
    size_t divergence_label = m_lts.add_action(action_label_string("d"));
    // Loop over all edges in the graph and add transitions to our lts
    boost::graph_traits<parity_game_t>::edge_iterator e, m;
    for(boost::tie(e,m) = edges(m_pg); e != m; ++e)
    {
      size_t s = vertex_state_map[source(*e, m_pg)];
      size_t t = vertex_state_map[target(*e, m_pg)];
      if (s == t) // If we have a self-loop, add divergence-label
      {
        m_lts.add_transition(transition(s,divergence_label,t));
      } 
      else
      {
        m_lts.add_transition(transition(s,label,t));
      }
    }
  }

  void run(std::string file)
  {
    reduce_pg_scc(file);
    convert_pg_to_ks();
    m_lts.save(file);
    // Call algorithm on m_lts.
    mCRL2log(verbose) << "Calling lts, with " << m_lts.num_states() << " states and " << m_lts.num_transitions() << " transitions." << std::endl;
    lts::detail::bisim_partitioner_gjkw_kripke(m_lts, true, true);
    mCRL2log(verbose) << "Replacing LTS!" << std::endl;
     convert_ks_to_pg();
    mCRL2log(verbose) << "Printed file" << std::endl;
  }
};

} // End mcrl2

using mcrl2::utilities::tools::input_output_tool;

template<class ParityGame>
class pgsolve_tool: public input_output_tool
{
  private:
    std::ifstream m_ifstream;
    std::ofstream m_ofstream;

  protected:
    typedef input_output_tool super;

  public:

    pgsolve_tool()
      : super(
          "pgsolve_gjkw",
          "Koen Degeling",
          "Reduce a parity game using the new m log n"
          "algorithm for stuttering bisimulation.",
          "Longer description"
      )
    {}

    bool run()
    {
      std::istream& is = open_input(input_filename(), m_ifstream);
      // Main program flow
      parity_game_t pg;
      mCRL2log(verbose) << "Loading PG from input file..." << std::endl;
      parse_pgsolver(pg, is, timer());
      mcrl2::pg_convert c = mcrl2::pg_convert(pg);
      c.run(output_filename());
      mCRL2log(verbose) << "Terminating" << std::endl;
      return true;
    }
};


int main(int argc, char* argv[])
{
  return pgsolve_tool<parity_game_t>().execute(argc, argv);
}