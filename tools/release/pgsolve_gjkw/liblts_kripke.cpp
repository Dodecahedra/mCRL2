
#include "mcrl2/lts/lts_algorithm.h"
#include "pg.h"

using namespace mcrl2::lts;
using namespace mcrl2::lts::detail;
using namespace mcrl2::log;
using namespace mcrl2::lts::detail::bisim_gjkw;

namespace mcrl2 {

// Class containing the new initialise_helper
template<class LTS_TYPE>
class initialise_helper
  : public bisim_gjkw::bisim_partitioner_gjkw_initialise_helper<LTS_TYPE>
{
  private:
    LTS_TYPE& aut;
    state_type nr_of_states;
    const state_type orig_nr_of_states;
    trans_type nr_of_transitions;

    std::vector<state_type> noninert_out_per_state, inert_out_per_state;
    std::vector<state_type> noninert_in_per_state, inert_in_per_state;
    std::vector<state_type> noninert_out_per_block, inert_out_per_block;
    std::vector<state_type> states_per_block;
    state_type nr_of_nonbottom_states;

    typedef bisim_gjkw::bisim_partitioner_gjkw_initialise_helper<LTS_TYPE> super;

  public:
    initialise_helper(LTS_TYPE& l, const bool branching,
                                              const bool divergence_preserving)
      : bisim_gjkw::bisim_partitioner_gjkw_initialise_helper<LTS_TYPE>(
        l, branching, divergence_preserving),
        aut(l),
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
    {}

    inline void init_kripke(bool const branching, bool const preserve_divergence)
    {
      // Initialisation...
      for (const transition& t: aut.get_transitions())
      {
        ++inert_in_per_state[t.to()];
        if (1 == ++inert_out_per_state[t.from()])
        {
            // this is the first inert outgoing transition of t.from()
            ++nr_of_nonbottom_states;
        }
        ++inert_out_per_block[0];
      }
    }

    inline void init_transitions(part_state_t& part_st, part_trans_t& part_tr,
                          bool const branching, bool const preserve_divergence)
{                                                                               assert(part_st.state_size() == super::get_nr_of_states());
                                                                                assert(part_tr.trans_size() == super::get_nr_of_transitions());
    // initialise blocks and B_to_C slices
    permutation_iter_t begin = part_st.permutation.begin();
    constln_t* const constln = new constln_t(super::get_nr_of_states(), begin,
                              part_st.permutation.end(), part_tr.B_to_C_end());
    if (1 < states_per_block.size())
    {
        constln->make_nontrivial();
    }
    std::vector<block_t*> blocks(states_per_block.size());
    B_to_C_iter_t B_to_C_begin = part_tr.B_to_C.begin();
    for (state_type B = 0; B < states_per_block.size(); ++B)
    {
        permutation_iter_t const end = begin + states_per_block[B];
        blocks[B] = new block_t(constln, begin, end);
        if (0 == noninert_out_per_block[B] && 0 == inert_out_per_block[B])
        {
            blocks[B]->set_inert_begin_and_end(part_tr.B_to_C.begin(),
                                                       part_tr.B_to_C.begin()); assert(blocks[B]->to_constln.empty());
        }
        else
        {
            blocks[B]->set_inert_begin_and_end(B_to_C_begin +
                                                     noninert_out_per_block[B],
                B_to_C_begin + noninert_out_per_block[B] +
                                                       inert_out_per_block[B]);
            blocks[B]->to_constln.emplace_back(B_to_C_begin,
                                                       blocks[B]->inert_end());
            B_to_C_desc_iter_t const slice =
                                        std::prev(blocks[B]->to_constln.end()); assert(B_to_C_begin < slice->end);
            for (; slice->end != B_to_C_begin; ++B_to_C_begin)
            {
                B_to_C_begin->B_to_C_slice = slice;
            }
        }
        begin = end;
    }                                                                           assert(part_st.permutation.end() == begin);
    /* only block 0 has a sequence number and non-bottom states:             */ assert(part_tr.B_to_C.end() == B_to_C_begin);
    blocks[0]->assign_seqnr();
    blocks[0]->set_bottom_begin(blocks[0]->begin() + nr_of_nonbottom_states);
    blocks[0]->set_marked_nonbottom_begin(blocks[0]->bottom_begin());

    // initialise states and succ slices
    part_st.state_info.front().set_pred_begin(part_tr.pred.begin());
    part_st.state_info.front().set_succ_begin(part_tr.succ.begin());
    for (state_type s = 0; super::get_nr_of_states() != s; ++s)
    {
        part_st.state_info[s].set_pred_end(part_st.state_info[s].pred_begin() +
                             noninert_in_per_state[s] + inert_in_per_state[s]);
        part_st.state_info[s].set_inert_pred_begin(part_st.state_info[s].
                                      pred_begin() + noninert_in_per_state[s]);
        // part_st.state_info[s+1].set_pred_begin(part_st.state_info[s].
        //                                                         pred_end());

        succ_iter_t succ_iter = part_st.state_info[s].succ_begin();
        succ_iter_t succ_end = succ_iter +
                            noninert_out_per_state[s] + inert_out_per_state[s];
        part_st.state_info[s].set_succ_end(succ_end);
        part_st.state_info[s].set_current_constln(succ_end);
        part_st.state_info[s].set_inert_succ_begin_and_end(part_st.
            state_info[s].succ_begin() + noninert_out_per_state[s], succ_end);
        if (succ_iter < succ_end)
        {
            --succ_end;
            for (; succ_iter < succ_end; ++succ_iter)
            {
                succ_iter->set_slice_begin_or_before_end(succ_end);
            }                                                                   assert(succ_iter == succ_end);
            succ_end->set_slice_begin_or_before_end(
                                           part_st.state_info[s].succ_begin());
        }                                                                       else  assert(succ_end == succ_iter);
        if (s < aut.num_states())
        {
            // s is not an extra Kripke state.  It is in block 0.
            part_st.state_info[s].block = blocks[0];
            if (0 != inert_out_per_state[s])
            {
                /* non-bottom state:                                         */ assert(0 != nr_of_nonbottom_states);
                --nr_of_nonbottom_states;
                part_st.state_info[s].pos = blocks[0]->begin() +
                                                        nr_of_nonbottom_states;
            }
            else
            {                                                                   // The following assertion is incomplete; only the second
                // bottom state:                                                // assertion (after the assignment) makes sure that not too
                                                                                // many states become part of this slice.
                                                                                assert(0 != states_per_block[0]);
                --states_per_block[0];
                part_st.state_info[s].pos = blocks[0]->begin() +
                                                           states_per_block[0]; assert(part_st.state_info[s].pos >= blocks[0]->bottom_begin());
            }
            *part_st.state_info[s].pos = &part_st.state_info[s];
            // part_st.state_info[s].notblue = 0;
        }
    }

    // initialise transitions (and finalise extra Kripke states)
    for (const transition& t: aut.get_transitions())
    {
      /* inert transition from t.from() to t.to()                      */ assert(0 != inert_in_per_state[t.to()]);
      --inert_in_per_state[t.to()];
      pred_iter_t const t_pred =
                          part_st.state_info[t.to()].inert_pred_begin() +
                                              inert_in_per_state[t.to()]; assert(0 != inert_out_per_state[t.from()]);
      --inert_out_per_state[t.from()];
      succ_iter_t const t_succ =
                        part_st.state_info[t.from()].inert_succ_begin() +
                                            inert_out_per_state[t.from()]; assert(0 != inert_out_per_block[0]);
      --inert_out_per_block[0];
      B_to_C_iter_t const t_B_to_C = blocks[0]->inert_begin() +
                                                  inert_out_per_block[0];

      t_pred->source = &part_st.state_info[t.from()];
      t_pred->succ = t_succ;
      t_succ->target = &part_st.state_info[t.to()];
      t_succ->B_to_C = t_B_to_C;
      // t_B_to_C->B_to_C_slice = (already initialised);
      t_B_to_C->pred = t_pred;
    }
    noninert_out_per_state.clear(); inert_out_per_state.clear();
    noninert_in_per_state.clear();  inert_in_per_state.clear();
    noninert_out_per_block.clear(); inert_out_per_block.clear();
    states_per_block.clear();

    aut.clear_transitions();

    mCRL2log(log::verbose) << "Size of the resulting Kripke structure: "
                               << super::get_nr_of_states() << " states and "
                               << super::get_nr_of_transitions() << " transitions.\n";
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
      init_helper.init_kripke(true, true);
      bisim_partitioner_gjkw<LTS_TYPE>::create_initial_partition_gjkw(branching, divergence_preserving);
      bisim_partitioner_gjkw<LTS_TYPE>::refine_partition_until_it_becomes_stable_gjkw();
    }
    ~liblts_kripke()
    {}
};
}