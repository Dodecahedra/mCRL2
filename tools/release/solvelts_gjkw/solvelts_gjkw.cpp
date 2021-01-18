// Author(s): Koen Degeling
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file solvelts_gjkw.cpp
/// Reads a parity game from input and converts it to an LTS and then
/// reduces it using Divergence-preserving branching bisimulation
/// equivalence with the O(m log m) algorithm
/// [Groote/Jansen/Keiren/Wijs 2017]

#include "mcrl2/lts/lts_fsm.h"
#include "mcrl2/utilities/input_output_tool.h"
#include "pg.h"
#include "pgsolver_io.h"
#include "utilities.h"
#include "mcrl2/lts/lts_algorithm.h"
#include <unordered_map>

#include <boost/graph/strong_components.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

using namespace mcrl2::lts;
using namespace mcrl2::log;

namespace mcrl2 {

class pg_convert
{

  protected:
  /** Storing Parity Game */
  parity_game_t m_pg;
  /** Reduced Kripke structure in LTS */
  lts::lts_aut_t m_lts;

  public:
  explicit pg_convert(parity_game_t const& pg)
      : m_pg(pg)
  {}

  void reduce_pg_scc()
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
    std::map<size_t, size_t> svertex_state_map;
    std::unordered_map<std::string, size_t> label_map;

    boost::graph_traits<parity_game_t >::vertex_iterator i, n;
    // We loop over vertices and add new states to the lts
    size_t bot = m_lts.add_action(action_label_string("bot"));
    label_map.insert(std::pair<std::string,size_t>("bot", bot));
    for(boost::tie(i, n) = vertices(m_pg); i != n; ++i)
    {
      size_t state = m_lts.add_state();
      size_t shadow_state = m_lts.add_state();
      // Add the (vertex,state) pair to map
      vertex_state_map[*i] = state;
      svertex_state_map[*i] = shadow_state;
      // Add bottom transitions
      m_lts.add_transition(transition(state, bot, shadow_state));
      // Add label transition
      size_t player = m_pg[*i].player;
      size_t priority = m_pg[*i].prio;
      mCRL2log(verbose) << "player: " << player << " prio " << priority << std::endl;
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
        auto label = label_map.find(string_label.str());
        label_no = label->second;
      }
      m_lts.add_transition(transition(shadow_state, label_no, state));
      mCRL2log(verbose) << "Added state: " << state << ", " << *i << std::endl;
    }
    m_lts.set_initial_state(0);
    size_t tau = m_lts.add_action(action_label_string("tau"));
    m_lts.apply_hidden_label_map(tau);
    // Loop over all edges in the graph and add transitions to our lts
    boost::graph_traits<parity_game_t>::edge_iterator e, m;
    for(boost::tie(e,m) = edges(m_pg); e != m; ++e)
    {
      size_t src = source(*e, m_pg);
      size_t trgt = target(*e, m_pg);
      size_t s = vertex_state_map[src];
      size_t t = vertex_state_map[trgt];
      if(m_pg[src].player == m_pg[trgt].player && m_pg[src].prio == m_pg[trgt].prio)
      { // Add a tau transition
        m_lts.add_transition(transition(s, tau, t));
      }
      else
      { // Add L(t) transition
        size_t player = m_pg[trgt].player;
        size_t priority = m_pg[trgt].prio;
        std::stringstream string_label;
        string_label << player << "," << priority;
        auto label = label_map.find(string_label.str());
        size_t label_no = label->second;
        m_lts.add_transition(transition(s, label_no, t));
      }

    }
  }

  void convert_lts_to_pg()
  /*  Convert reduced LTS in lts_aut to a Parity Game in
      parity_game_t format. */
  {
    // Clear old Parity Game and construct new one from m_lts
    m_pg.clear();
    std::map<size_t, size_t> state_vertex_map;

    for(const transition& tr: m_lts.get_transitions())
    { /* We loop over all transitions and add vertices for the source and target
         states if they did not yet exist. If transition contains a label, we add
         the label information to the target state.*/
      if(state_vertex_map.count(tr.from()) == 0)
      {
        size_t v = boost::add_vertex(m_pg);
        state_vertex_map[tr.from()] = v;
        m_pg[v].unused = true;
      }
      if(state_vertex_map.count(tr.to()) == 0)
      {
        size_t u = boost::add_vertex(m_pg);
        state_vertex_map[tr.to()] = u;
        m_pg[u].unused = true;
      }
      boost::add_edge(state_vertex_map[tr.from()], state_vertex_map[tr.to()], m_pg);
      action_label_string string_label = m_lts.action_label(tr.label());
      mCRL2log(verbose) << "Label: " << string_label << std::endl;
      if(!(string_label == "tau" || string_label == "bot"))
      { // Transition with label L(t)
        size_t u = state_vertex_map[tr.to()];
        std::vector<size_t> labels = get_labels(string_label); // We get the labels from the action string
        if(m_pg[u].unused)
        { // If vertex was not yet initialised
          m_pg[u].player = labels[0] == 0 ? even : odd;
          m_pg[u].prio = labels[1];
          m_pg[u].unused = false; // Mark vertex as initialised
          mCRL2log(verbose) << "Set vertex: " << u << " to: " << labels[0] <<
                                                  " " << labels[1] << std::endl;
        }
      }
    }

    size_t i = 0;
    while(i < boost::num_vertices(m_pg))
    { /* We loop over all the vertices and remove those which have not been
         initialised. */
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

  std::vector<size_t> get_labels(std::string str)
  /* Auxiliary function to get a vector with the {player,priority} from the action
     string. */
  {
    std::vector<size_t> labels;
    size_t i = str.find(',');
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

  void run(std::string file, std::string infile)
  {
    mCRL2log(verbose) << "Filename: " << infile << std::endl;
    mcrl2::utilities::execution_timer timer(infile,
               "/home/koen/Documents/TUe/Master/Quartile2/2IMF00Seminar/mCRL2/timing");
    timer.start("execution");
    reduce_pg_scc(); // Make sure we remove any cycles in the graph.
    convert_pg_to_lts(); // Convert parity game to LTS.
    /* Call GJKW algorithm and start timers for reducing and execution. */
    timer.start("reducing");
    lts::detail::bisim_partitioner_gjkw<lts_aut_t> part(m_lts, true, true);
    part.replace_transition_system(true, true);
    timer.finish("reducing");
    convert_lts_to_pg(); // Convert LTS to Parity Game.
    std::ofstream m_ofstream;
    std::ostream& os = open_output(file, m_ofstream);
    print_pgsolver(m_pg, os);
    timer.finish("execution");
    timer.report();
  }
};

}

using mcrl2::utilities::tools::input_output_tool;

class solvelts_tool: public input_output_tool
{
  private:
  std::ifstream m_ifstream;

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
    std::string infile = input_filename();
    std::istream& is = open_input(infile, m_ifstream);
    parity_game_t pg;
    parse_pgsolver(pg, is, timer());
    mcrl2::pg_convert c = mcrl2::pg_convert(pg);
    c.run(output_filename(), infile);
    mCRL2log(verbose) << "Terminating..." << std::endl;
    return true;
  }
};


int main(int argc, char* argv[])
{
  return solvelts_tool().execute(argc, argv);
}