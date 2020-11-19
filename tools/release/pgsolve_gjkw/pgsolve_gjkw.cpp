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
#include "mcrl2/bes/pbes_input_output_tool.h"

#include "mcrl2/bes/parse.h"
#include "mcrl2/lts/lts_algorithm.h"

using namespace mcrl2::lts; // For reduce() function.
using namespace mcrl2::bes; // For parsing a PG from input.
using namespace mcrl2::utilities::tools;
using namespace mcrl2::log;
using mcrl2::bes::tools::bes_input_output_tool;

namespace mcrl2 {

class pg_convert
{
  protected:
    /** Storing bes from Parity Game */
    boolean_equation_system m_bes;
    /** Reduced Kripke structure in lts */
    lts::lts_aut_t m_lts;

  public:
  pg_convert(boolean_equation_system& v_bes)
  {
    m_bes = v_bes;
  }

  void convert_bes()
  {
    // Convert reduced lts to BES
  }

  void convert_lts()
  {
    // Convert BES to kripke structure in lts
  }

  void run()
  {
    mCRL2log(verbose) << "Starting conversion to LTS..." << std::endl;
    convert_bes();
    mCRL2log(verbose) << "Reducing LTS using m log n algorithm" << std::endl;
    reduce(m_lts, lts_eq_divergence_preserving_branching_bisim_gjkw);
    mCRL2log(verbose) << "Converting to BES" << std::endl;
    convert_lts();
  }
};

} // End mcrl2


typedef bes_input_output_tool<input_output_tool> super;
class pgsolve_tool: public super
{
  
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
      // Main program flow
      boolean_equation_system b;
      mCRL2log(verbose) << "Loading BES from input file..." << std::endl;
      load_bes(b, input_filename(), bes_input_format());
      // mcrl2::pg_convert(b).run();
      mCRL2log(verbose) << "Saving output to file..." << std::endl;
      save_bes(b, output_filename(), bes_output_format());
      mCRL2log(verbose) << "Terminating..." << std::endl;
      return true;
    }
};


int main(int argc, char* argv[])
{
  return pgsolve_tool().execute(argc, argv);
}