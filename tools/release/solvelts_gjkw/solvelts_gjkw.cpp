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

#include "mcrl2/lts/lts_fsm.h"
#include "mcrl2/utilities/input_output_tool.h"
#include "pg.h"
#include "pgsolver_io.h"
#include "utilities.h"
#include "mcrl2/lts/lts_algorithm.h"
#include "scc.h"
#include <unordered_map>

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
  lts::lts_aut_t m_lts;

  public:
  pg_convert(parity_game_t pg)
      : m_pg(pg)
  {}

  parity_game_t get_parity_game()
  {
    return m_pg;
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

  void convert_pg_to_lts()
  /* Convert Parity Game to an LTS of type lts_aut */
  {
    // LTS in which we save the Kripke structure
    std::map<size_t, size_t> vertex_state_map;

    boost::graph_traits<parity_game_t >::vertex_iterator i, n;
    // We loop over vertices and add new states to the lts
    for(boost::tie(i, n) = vertices(m_pg); i != n; ++i)
    {
      size_t state = m_lts.add_state();
      // Add the (vertex,state) pair to map
      vertex_state_map[*i] = state;
    }
    m_lts.set_initial_state(0);
    // Map to store the action label -> size_t
    std::unordered_map<std::string, size_t> label_map;
    // Loop over all edges in the graph and add transitions to our lts
    boost::graph_traits<parity_game_t>::edge_iterator e, m;
    for(boost::tie(e,m) = edges(m_pg); e != m; ++e)
    {
      size_t s = vertex_state_map[source(*e, m_pg)];
      size_t t = vertex_state_map[target(*e, m_pg)];
      size_t player = m_pg[s].player;
      size_t priority = m_pg[s].prio;
      std::stringstream string_label;
      string_label << player << "," << priority;
      size_t label_no;
      if(label_map.count(string_label.str()) == 0)
      { // Label does not yet exist
        label_no = m_lts.add_action(action_label_string(string_label.str()));
        label_map.insert(std::pair<std::string,size_t>(string_label.str(), label_no));
      }
      else
      {
        std::unordered_map<std::string,size_t>::const_iterator label =
                                              label_map.find(string_label.str());
        label_no = label->second;
      }
      m_lts.add_transition(transition(s, label_no, t));
    }
  }

  void convert_lts_to_pg()
  /*  Convert reduced LTS in lts_aut to a Parity Game in
      parity_game_t format. */
  {
    // Clear old Parity Game and construct new one from m_lts
    m_pg.clear();
    std::map<size_t, size_t> state_vertex_map;

    // Loop over states and add vertices to graph
    for(size_t i=0; i!=m_lts.num_states(); ++i )
    {
      size_t v = boost::add_vertex(m_pg);
      // Set properties of the new state
      // Safe the mapping from state to vertex for later when adding edges.
      state_vertex_map[i] = v;
    }

    for(const transition& tr: m_lts.get_transitions())
    {
      size_t s = state_vertex_map[tr.from()];
      size_t t = state_vertex_map[tr.to()];
      action_label_string string_label = m_lts.action_label(tr.label());
      std::vector<size_t> labels = get_labels(string_label); // We get the labels from the string
      if(m_pg[s].player == NULL) // If not yet set
      {
        m_pg[s].player = labels[0] == 0 ? even : odd;
        m_pg[s].prio = labels[1];
      }
      boost::add_edge(s, t, m_pg);
    }
  }

  std::vector<size_t> get_labels(std::string str)
  {
    std::vector<size_t> labels;
    size_t i = str.find(",");
    std::string player = str.substr(0, i);
    std::string priority = str.substr(i+1, str.length()-i-1);
    std::stringstream f(player);
    std::stringstream g(priority);
    size_t pl,pr;
    f >> pl;
    g >> pr;
    labels.push_back(pl);
    labels.push_back(pr);
    return labels;
  }

  void run(std::string file)
  {
    reduce_pg_scc(file); // TODO: Check if this is done in bisim_partitioner
    convert_pg_to_lts(); // We convert the parity game to an lts
    // Call algorithm on m_lts.
    mCRL2log(verbose) << "Calling lts, with " << m_lts.num_states() << " states and " << m_lts.num_transitions() << " transitions." << std::endl;
    lts::detail::bisim_partitioner_gjkw<lts_aut_t> part(m_lts, true, true); // We run the GJKW algorithm
    part.replace_transition_system(true, true);
    mCRL2log(verbose) << "Replacing LTS!" << std::endl;
    convert_lts_to_pg();
    std::ofstream m_ofstream;
    std::ostream& os = open_output(file, m_ofstream);
    print_pgsolver(m_pg, os);
    mCRL2log(verbose) << "Printed file" << std::endl;
  }
};

} // End mcrl2

using mcrl2::utilities::tools::input_output_tool;

class solvelts_tool: public input_output_tool
{
  private:
  std::ifstream m_ifstream;
  std::ofstream m_ofstream;

  protected:
  typedef input_output_tool super;

  public:
    solvelts_tool()
      : super(
      "solvelts",
      "Koen Degeling",
      "Converts a parity game to an lts and uses the GJKW algorithm"
      " to compute the quotient structure and revert back to parity game",
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
  return solvelts_tool().execute(argc, argv);
}