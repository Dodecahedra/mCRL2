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

using namespace mcrl2::lts; // For reduce() function.
using namespace mcrl2::log;

// namespace mcrl2 {

// template<class LTS_TYPE>
// class pg_convert  : public lts::detail::bisim_gjkw::bisim_partitioner_gjkw_initialise_helper
// {

//   protected:
//     /** Storing bes from Parity Game */
//     boolean_equation_system m_bes;
//     /** Reduced Kripke structure in lts */
//     lts::lts_aut_t m_lts;

//   public:
  
//   bisim_partitioner_gjkw_initialise_helper(LTS_TYPE& l, bool const branching,
//                                                 bool const preserve_divergence)
//   : aut(l),
//     nr_of_states(l.num_states()),
//     orig_nr_of_states(l.num_states()),
//     nr_of_transitions(l.num_transitions()),
//     noninert_out_per_state(l.num_states(), 0),
//     inert_out_per_state(l.num_states(), 0),
//     noninert_in_per_state(l.num_states(), 0),
//     inert_in_per_state(l.num_states(), 0),
//     noninert_out_per_block(1, 0),
//     inert_out_per_block(1, 0),
//     states_per_block(1, l.num_states()),
//     nr_of_nonbottom_states(0)
//   {
    
//   }

//   void run()
//   {
//     mCRL2log(verbose) << "Starting conversion to LTS..." << std::endl;
//     convert_bes();
//     mCRL2log(verbose) << "Reducing LTS using m log n algorithm" << std::endl;
//     reduce(m_lts, lts_eq_divergence_preserving_branching_bisim_gjkw);
//     mCRL2log(verbose) << "Converting to BES" << std::endl;
//     convert_lts();
//   }
// };

// } // End mcrl2

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
      std::ostream& os = open_output(output_filename(), m_ofstream);
      // Main program flow
      parity_game_t pg;
      mCRL2log(verbose) << "Loading PG from input file..." << std::endl;
      parse_pgsolver(pg, is, timer());
      // mcrl2::pg_convert(b).run();
      mCRL2log(verbose) << "Saving output to file..." << std::endl;
      print_pgsolver(pg, os);
      mCRL2log(verbose) << "Terminating..." << std::endl;
      return true;
    }
};


int main(int argc, char* argv[])
{
  return pgsolve_tool<parity_game_t>().execute(argc, argv);
}