
#include "mcrl2/lts/lts_algorithm.h"
#include "pg.h"

using namespace mcrl2::lts;
using namespace mcrl2::lts::detail;
using namespace mcrl2::log;

namespace mcrl2 {


// Class containing the new initialise_helper
template<class LTS_TYPE>
class initialise_helper
  : public bisim_partitioner_gjkw_initialise_helper<LTS_TYPE>
{
  private:
    LTS_TYPE& aut;
    state_type nr_of_states;
    const state_type orig_nr_of_states;
    trans_type nr_of_transitions;

    // key and hash function for (action, target state) pair. Required since
    // unordered_map does not directly allow to use pair keys
    class Key 
    {
      public:
        label_type first;
        state_type second;

        Key(const label_type& f, const state_type& s)
          : first(f),
            second(s)
        {} 

        bool operator==(const Key &other) const 
        {
            return first == other.first && second == other.second;
        }
    };

    class KeyHasher 
    {
      public:
        std::size_t operator()(const Key& k) const
        {
            return std::hash<label_type>()(k.first) ^
                                      (std::hash<state_type>()(k.second) << 1);
        }
    };
    // Map used to convert LTS to Kripke structure
    // (also used when converting Kripke structure back to LTS)
    std::unordered_map<Key, state_type, KeyHasher> extra_kripke_states;

    // temporary map to keep track of blocks. maps transition labels (different
    // from tau) to blocks
    std::unordered_map<label_type, state_type> action_block_map;

    std::vector<state_type> noninert_out_per_state, inert_out_per_state;
    std::vector<state_type> noninert_in_per_state, inert_in_per_state;
    std::vector<state_type> noninert_out_per_block, inert_out_per_block;
    std::vector<state_type> states_per_block;
    state_type nr_of_nonbottom_states;
  public:
    initialise_helper(LTS_TYPE& l, bool const branching, 
                                  bool const preserve_divergence)
      : aut(l),
      nr_of_states(l.num_states()),
      orig_nr_of_states(l.num_states()),
      nr_of_transitions(l.num_transitions()),
      noninert_out_per_state(l.num_states(), 0),
      inert_out_per_state(l.num_states(), 0),
      noninert_in_per_state(l.num_states(), 0),
      inert_in_per_state(l.num_states(), 0),
      noninert_out_per_block(1, 0),
      inert_out_per_block(1, 0),
      states_per_block(1, l.num_states()),
      nr_of_nonbottom_states(0)
  {
    extra_kripke_states.clear();
    mCRL2log(verbose) << "Using the new constructor!" << std::endl;
    /* Use empty constructor, as we are already passing a Kripke Structure in 
       lts_fsm_t format. */
  }
};

// Class containig the new implementation of the bisim partitioner
template<class LTS_TYPE>
class liblts_kripke : public bisim_partitioner_gjkw<LTS_TYPE>
{
  protected:
    initialise_helper<LTS_TYPE> init_helper;
    bisim_gjkw::bisim_partitioner_gjkw_initialise_helper<LTS_TYPE>::part_state_t part_st;
    bisim_gjkw::bisim_partitioner_gjkw_initialise_helper<LTS_TYPE>::part_trans_t part_tr;
  public:
    // The constructor constructs the data structures and immediately
    // calculates the bisimulation quotient.  However, it does not change the
    // LTS.
    liblts_kripke(lts_fsm_t& l, bool branching = false,
                                        bool preserve_divergence = false)
      : init_helper(l, branching, preserve_divergence),
        part_st(init_helper.get_nr_of_states()),
        part_tr(init_helper.get_nr_of_transitions())
    {                                                                                assert(branching || !preserve_divergence);
      bisim_partitioner<LTS_TYPE>::create_initial_partition_gjkw(branching, preserve_divergence);
      bisim_partitioner<LTS_TYPE>::refine_partition_until_it_becomes_stable_gjkw();
      mCRL2log(verbose) << "Calling new partitioner." << std::endl;
      mCRL2log(verbose) << "Created initial partition" << std::endl;
      mCRL2log(verbose) << "Refined partition" << std::endl;
    }
};

}
