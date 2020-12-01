
#include "mcrl2/lts/lts_algorithm.h"
#include "pg.h"

using namespace mcrl2::lts;
using namespace mcrl2::lts::detail;
using namespace mcrl2::log;

namespace mcrl2 {


// Class containing the new initialise_helper
template<class LTS_TYPE>
class initialise_helper
  : public bisim_gjkw::bisim_partitioner_gjkw_initialise_helper<LTS_TYPE>
{
  public:
    initialise_helper(LTS_TYPE& l, const bool branching,
                                              const bool divergence_preserving)
      : bisim_gjkw::bisim_partitioner_gjkw_initialise_helper<LTS_TYPE>(
        l, branching, divergence_preserving)
    {
      mCRL2log(verbose) << "Using the new initialiser!" << std::endl;
    }

    inline void init_kripke(bool const branching, bool const preserve_divergence)
    {
      mCRL2log(verbose) << "Using the new init_kripke function!" << std::endl;
    }

};

template<class LTS_TYPE>
class liblts_kripke
  : public bisim_partitioner_gjkw<LTS_TYPE>
{
  private:
    initialise_helper<LTS_TYPE> init_helper;
  
  public:
    liblts_kripke(LTS_TYPE& l, const bool branching,
                                      const bool divergence_preserving)
      : bisim_partitioner_gjkw<LTS_TYPE>(l, branching, divergence_preserving),
        init_helper(l, branching, divergence_preserving)
    {
      mCRL2log(verbose) << "This is the new constructor!" << std::endl;
      init_helper.init_kripke(true, false);
      bisim_partitioner_gjkw<LTS_TYPE>::create_initial_partition_gjkw(branching, divergence_preserving);
      bisim_partitioner_gjkw<LTS_TYPE>::refine_partition_until_it_becomes_stable_gjkw();
    }
    ~liblts_kripke()
    {}
};

}
