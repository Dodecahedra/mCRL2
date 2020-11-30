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

#include "mcrl2/utilities/input_output_tool.h"
#include "pg.h"
#include "pgsolver_io.h"
#include "utilities.h"
#include "mcrl2/lts/lts_algorithm.h"
#include "liblts_kripke.cpp"

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

using namespace mcrl2::lts; // For reduce() function.
using namespace mcrl2::log; // Logging.

namespace mcrl2 {

class pg_convert
{

  protected:
    /** Storing bes from Parity Game */
    parity_game_t m_pg;
    /** Reduced Kripke structure in lts */
    lts::lts_fsm_t m_lts;

  public:
  pg_convert(parity_game_t pg)
    : m_pg(pg)
    {}

  parity_game_t get_parity_game()
  {
    return m_pg;
  }

  void convert_pg()
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
  
  void convert_ks()
  /* Convert Parity Game to a Kripke structure contained in lts_fsm_t. */
  {
    // LTS in which we save the Kripke structure
    std::map<size_t, size_t> vertex_state_map;

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
    // Set an initial state.
    size_t label = m_lts.add_action(action_label_string(" "));
    // Loop over all edges in the graph and add transitions to our lts
    boost::graph_traits<parity_game_t>::edge_iterator e, m;
    for(boost::tie(e,m) = edges(m_pg); e != m; ++e)
    {
      size_t s = vertex_state_map[source(*e, m_pg)];
      size_t t = vertex_state_map[target(*e, m_pg)];
      m_lts.add_transition(transition(s,label,t));
    }
  }

  void run(std::string file)
  {
    std::ofstream m_ofstream;
    convert_ks();
    // Call algorithm on m_lts.

    liblts_kripke<lts_fsm_t>(m_lts, true, false);
    convert_pg();
    std::ostream& os = open_output(file, m_ofstream);
    print_pgsolver(m_pg, os);
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
      mCRL2log(verbose) << "Saving output to file..." << std::endl;
      mCRL2log(verbose) << "Terminating..." << std::endl;
      return true;
    }
};


int main(int argc, char* argv[])
{
  return pgsolve_tool<parity_game_t>().execute(argc, argv);
}