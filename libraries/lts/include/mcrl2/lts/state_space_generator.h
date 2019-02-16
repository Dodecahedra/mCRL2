// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lts/state_space_generator.h
/// \brief add your file description here.

#ifndef MCRL2_LTS_STATE_SPACE_GENERATOR_H
#define MCRL2_LTS_STATE_SPACE_GENERATOR_H

#include "mcrl2/lps/explorer.h"
#include "mcrl2/lts/lts_builder.h"
#include "mcrl2/trace/trace.h"

namespace mcrl2 {

namespace lts {

inline
std::ostream& operator<<(std::ostream& out, const lps::state& s)
{
  return out << atermpp::pp(s);
}

namespace detail {

inline
bool save_trace(
  trace::Trace& tr,
  const std::string& filename
)
{
  try
  {
    tr.save(filename);
    mCRL2log(log::info) << " and saved trace to '" << filename << "'";
    return true;
  }
  catch(...)
  {
    mCRL2log(log::info) << ", but its trace could not be saved to '" << filename << "'";
  }
  return false;
}

void save_traces(
  trace::Trace& tr,
  const std::string& filename1,
  trace::Trace& tr2,
  const std::string& filename2
)
{
  try
  {
    tr.save(filename1);
    tr2.save(filename2);
    mCRL2log(log::info) << " and saved traces to '" << filename1 << "' and '" << filename2 << "'";
  }
  catch(...)
  {
    mCRL2log(log::info) << ", but its traces could not be saved to '" << filename1 << "' and '" << filename2 << "'";
  }
}

// Facility for constructing a trace to a given state.
class trace_constructor
{
  protected:
    lps::explorer& explorer;
    std::map<lps::state, lps::state> backpointers;

    // Finds a transition s0 --a--> s1, and returns a.
    lps::multi_action find_action(const lps::state& s0, const lps::state& s1)
    {
      for (const std::pair<lps::multi_action, lps::state>& t: explorer.generate_transitions(s0))
      {
        if (t.second == s1)
        {
          return t.first;
        }
      }
      throw mcrl2::runtime_error("no transition found in find_action");
    }

  public:
    explicit trace_constructor(lps::explorer& explorer_)
      : explorer(explorer_)
    {}

    // Constructs a trace ending in s, using the backpointers map.
    trace::Trace construct_trace(const lps::state& s)
    {
      std::deque<lps::state> states{ s };
      std::deque<lps::multi_action> actions;
      while (true)
      {
        const lps::state& s1 = states.front();
        auto i = backpointers.find(s1);
        if (i == backpointers.end())
        {
          break;
        }
        const lps::state& s0 = i->second;
        states.push_front(s0);
        actions.push_front(find_action(s0, s1));
      }

      trace::Trace tr;
      for (std::size_t i = 0; i < actions.size(); i++)
      {
        tr.setState(states[i]);
        tr.addAction(actions[i]);
      }
      tr.setState(states.back());
      return tr;
    }

    // Adds a back pointer for the given edge
    void add_edge(const lps::state& s0, const lps::state& s1)
    {
      backpointers[s1] = s0;
    }

    void clear()
    {
      backpointers.clear();
    }
};

class action_detector
{
  protected:
    const std::set<core::identifier_string>& trace_actions;
    const std::set<lps::multi_action>& trace_multiactions;
    trace_constructor& m_trace_constructor;
    const std::string& filename_prefix;
    std::vector<bool> summand_matches;
    std::size_t m_trace_count = 0;
    std::size_t m_max_trace_count;

    bool match_action(const lps::action_summand& summand) const
    {
      using utilities::detail::contains;
      for (const process::action& a: summand.multi_action().actions())
      {
        if (contains(trace_actions, a.label().name()))
        {
          return true;
        }
      }
      return false;
    }

    bool match_summand(std::size_t i) const
    {
      return summand_matches[i];
    }

    std::string create_filename(const lps::multi_action& a)
    {
      using utilities::detail::contains;
      std::string filename = filename_prefix + "_act_" + std::to_string(m_trace_count++);
      if (utilities::detail::contains(trace_multiactions, a))
      {
        filename += "_" + lps::pp(a);
      }
      for (const process::action& a_i: a.actions())
      {
        if (utilities::detail::contains(trace_actions, a_i.label().name()))
        {
          filename += "_" + core::pp(a_i.label().name());
        }
      }
      filename = filename + ".trc";
      return filename;
    }

  public:
    template <typename Specification>
    action_detector(
      const Specification& lpsspec,
      trace_constructor& trace_constructor_,
      const std::set<core::identifier_string>& trace_actions_,
      const std::set<lps::multi_action>& trace_multiactions_,
      const std::string& filename_prefix_,
      std::size_t max_trace_count
    )
      : trace_actions(trace_actions_),
        trace_multiactions(trace_multiactions_),
        m_trace_constructor(trace_constructor_),
        filename_prefix(filename_prefix_),
        m_max_trace_count(max_trace_count)
    {
      using utilities::detail::contains;
      const auto& summands = lpsspec.process().action_summands();
      summand_matches.reserve(summands.size());
      for (const auto& summand: summands)
      {
        summand_matches.push_back(match_action(summand));
      }
    }

    bool detect_action(const lps::state& s0, std::size_t s0_index, const lps::multi_action& a, const lps::state& s1, std::size_t summand_index)
    {
      using utilities::detail::contains;
      if (!match_summand(summand_index))
      {
        return false;
      }
      bool result = false;

      mCRL2log(log::info) << "Action '" + lps::pp(a) + "' found (state index: " + std::to_string(s0_index) + ")";
      if (m_trace_count < m_max_trace_count)
      {
        trace::Trace tr = m_trace_constructor.construct_trace(s0);
        tr.setState(s1);
        tr.addAction(a);
        std::string filename = create_filename(a);
        save_trace(tr, filename);
        result = true;
      }
      mCRL2log(log::info) << ".\n";
      return result;
    }
};

class deadlock_detector
{
  protected:
    trace_constructor& m_trace_constructor;
    const std::string& filename_prefix;
    std::size_t m_trace_count = 0;
    std::size_t m_max_trace_count;

  public:
    deadlock_detector(
      trace_constructor& trace_constructor_,
      const std::string& filename_prefix_,
      std::size_t max_trace_count
    )
      : m_trace_constructor(trace_constructor_),
        filename_prefix(filename_prefix_),
        m_max_trace_count(max_trace_count)
    {}

    void detect_deadlock(const lps::state& s, std::size_t s_index)
    {
      mCRL2log(log::info) << "Deadlock found (state index: " + std::to_string(s_index) + ")";
      if (m_trace_count < m_max_trace_count)
      {
        trace::Trace tr = m_trace_constructor.construct_trace(s);
        std::string filename = filename_prefix + "_dlk_" + std::to_string(m_trace_count++) + ".trc";
        save_trace(tr, filename);
      }
      mCRL2log(log::info) << ".\n";
    }
};

class nondeterminism_detector
{
  protected:
    trace_constructor& m_trace_constructor;
    const std::string& filename_prefix;
    std::map<lps::multi_action, lps::state> transitions;
    std::size_t m_trace_count = 0;
    std::size_t m_max_trace_count;

  public:
    nondeterminism_detector(
      trace_constructor& trace_constructor_,
      const std::string& filename_prefix_,
      std::size_t max_trace_count = 0
    )
      : m_trace_constructor(trace_constructor_),
        filename_prefix(filename_prefix_),
        m_max_trace_count(max_trace_count)
    {}

    void start_state()
    {
      transitions.clear();
    }

    bool detect_nondeterminism(const lps::state& s0, std::size_t s0_index, const lps::multi_action& a, const lps::state& s1)
    {
      bool result = false;
      auto i = transitions.find(a);
      if (i == transitions.end())
      {
        transitions.insert(std::make_pair(a, s1));
      }
      else if (i->second != s1) // nondeterminism detected
      {
        mCRL2log(log::info) << "Nondeterministic state found (state index: " + std::to_string(s0_index) + ")";
        if (m_trace_count < m_max_trace_count)
        {
          trace::Trace tr = m_trace_constructor.construct_trace(s0);
          tr.setState(s1);
          tr.addAction(a);
          std::string filename = filename_prefix + "_nondeterministic_" + std::to_string(m_trace_count++) + ".trc";
          save_trace(tr, filename);
          result = true;
        }
        mCRL2log(log::info) << ".\n";
      }
      return result;
    }
};

class divergence_detector
{
  protected:
    lps::explorer& explorer;
    const std::string& filename_prefix;
    trace_constructor m_local_trace_constructor;
    std::unordered_map<lps::state, std::size_t> m_divergent_states;
    std::vector<lps::explorer::explorer_summand> m_regular_summands;
    std::vector<lps::explorer::explorer_summand> m_confluent_summands;
    std::size_t m_trace_count = 0;
    std::size_t m_max_trace_count;
    bool m_timed;

  public:
    divergence_detector(
      lps::explorer& explorer_,
      bool timed,
      const std::set<core::identifier_string>& actions,
      const std::string& filename_prefix_,
      std::size_t max_trace_count
    )
      : explorer(explorer_),
        filename_prefix(filename_prefix_),
        m_local_trace_constructor(explorer),
        m_max_trace_count(max_trace_count),
        m_timed(timed)
    {
      using utilities::detail::contains;

      auto is_hidden = [&](const lps::explorer::explorer_summand& summand)
      {
        for (const process::action& a: summand.multi_action.actions())
        {
          if (!contains(actions, a.label().name()))
          {
            return false;
          }
        }
        return true;
      };

      for (const lps::explorer::explorer_summand& summand: explorer.regular_summands())
      {
        if (is_hidden(summand))
        {
          m_regular_summands.push_back(summand);
        }
      }

      for (const lps::explorer::explorer_summand& summand: explorer.confluent_summands())
      {
        if (is_hidden(summand))
        {
          m_confluent_summands.push_back(summand);
        }
      }
    }

    // Returns true if a trace was saved.
    bool detect_divergence(const lps::state& s, std::size_t s_index, trace_constructor& global_trace_constructor)
    {
      bool result = false;
      m_local_trace_constructor.clear();

      auto q = m_divergent_states.find(s);
      if (q != m_divergent_states.end())
      {
        std::string message = "Divergent state found (state index: " + std::to_string(s_index) + "), reachable from divergent state with index " + std::to_string(q->second);
        mCRL2log(log::info) << message << ".\n";
        m_divergent_states.erase(q);
        return false;
      }

      std::unordered_map<lps::state, std::size_t> discovered;
      const lps::state* source = nullptr;
      lps::state last_discovered;

      data::data_expression_list process_parameter_undo = explorer.process_parameter_values();
      explorer.generate_state_space(
        m_timed,
        true,
        s,
        m_regular_summands,
        m_confluent_summands,
        discovered,

        // discover_state
        [&](const lps::state& s, std::size_t /* s_index */)
        {
          last_discovered = s;
          if (source)
          {
            m_local_trace_constructor.add_edge(*source, s);
          }
        },

        // examine_transition
        [&](const lps::state& s0, std::size_t /* s0_index */, const process::timed_multi_action& a, const lps::state& s1, std::size_t /* s1_index */, std::size_t /* summand_index */)
        {
          if (s1 != last_discovered || s1 == s) // found a loop, hence s is a divergent state
          {
            mCRL2log(log::info) << "Divergent state found (state index: " + std::to_string(s_index) + ")";
            if (m_trace_count < m_max_trace_count)
            {
              trace::Trace tr = global_trace_constructor.construct_trace(s);
              trace::Trace tr_loop = m_local_trace_constructor.construct_trace(s0);
              for (const lps::state& u: tr_loop.states())
              {
                m_divergent_states[u] = s_index;
              }
              tr_loop.setState(s1);
              tr_loop.addAction(lps::multi_action(a.actions(), a.time()));
              std::string filename = filename_prefix + "_divergence_" + std::to_string(m_trace_count) + ".trc";
              std::string loop_filename = filename_prefix + "_divergence_loop" + std::to_string(m_trace_count++) + ".trc";
              save_traces(tr, filename, tr_loop, loop_filename);
              result = true;
            }
            mCRL2log(log::info) << ".\n";
            explorer.abort();
          }
        },

        // start_state
        [&](const lps::state& s, std::size_t /* s_index */)
        {
          source = &s;
        }
      );
      explorer.set_process_parameter_values(process_parameter_undo);
      return result;
    }
};

class progress_monitor
{
  protected:
    std::size_t level = 1;    // the current exploration level
    std::size_t level_up = 1; // when count reaches level_up, the level is increased
    std::size_t count = 0;
    std::size_t transition_count = 0;

    std::size_t last_state_count = 0;
    std::size_t last_transition_count = 0;
    time_t last_log_time = time(nullptr) - 1;
    time_t new_log_time;

    lps::exploration_strategy search_strategy;

  public:
    explicit progress_monitor(lps::exploration_strategy search_strategy_)
      : search_strategy(search_strategy_)
    {}

    void examine_transition()
    {
      transition_count++;
    }

    void finish_state(std::size_t state_count, std::size_t todo_list_size)
    {
      if (search_strategy == lps::es_breadth)
      {
        if (++count == level_up)
        {
          mCRL2log(log::debug) << "Number of states at level " << level << " is " << state_count - last_state_count << "\n";
          level++;
          level_up = count + todo_list_size;
          last_state_count = state_count;
          last_transition_count = transition_count;
        }

        if (time(&new_log_time) > last_log_time)
        {
          last_log_time = new_log_time;
          std::size_t lvl_states = state_count - last_state_count;
          std::size_t lvl_transitions = transition_count - last_transition_count;
          mCRL2log(log::status) << std::fixed << std::setprecision(2)
                                << state_count << "st, " << transition_count << "tr"
                                << ", explored " << 100.0 * ((float) count / state_count)
                                << "%. Last level: " << level << ", " << lvl_states << "st, " << lvl_transitions
                                << "tr.\n";
        }
      }
      else
      {
        if (++count % 1000 == 0)
        {
          mCRL2log(log::verbose) << "monitor: currently explored "
                            << count << " state" << ((count==1)?"":"s")
                            << " and " << transition_count << " transition" << ((transition_count==1)?"":"s")
                            << std::endl;
        }
      }
    }

    void finish_exploration(std::size_t state_count)
    {
      if (search_strategy == lps::es_breadth)
      {
        mCRL2log(log::verbose) << "done with state space generation ("
                               << level-1 << " level" << ((level==2)?"":"s") << ", "
                               << state_count << " state" << ((state_count == 1)?"":"s")
                               << " and " << transition_count << " transition" << ((transition_count==1)?"":"s") << ")" << std::endl;
      }
      else
      {
        mCRL2log(log::verbose) << "done with state space generation ("
                          << state_count << " state" << ((state_count == 1)?"":"s")
                          << " and " << transition_count << " transition" << ((transition_count==1)?"":"s") << ")" << std::endl;
      }
    }
};

} // namespace detail

struct state_space_generator
{
  const lps::explorer_options& options;
  lps::explorer explorer;
  detail::trace_constructor m_trace_constructor;

  detail::action_detector m_action_detector;
  detail::deadlock_detector m_deadlock_detector;
  detail::nondeterminism_detector m_nondeterminism_detector;
  std::unique_ptr<detail::divergence_detector> m_divergence_detector;
  detail::progress_monitor m_progress_monitor;
  bool m_timed;

  state_space_generator(const lps::specification& lpsspec, const lps::explorer_options& options_)
    : options(options_), 
      explorer(lpsspec, options_),
      m_trace_constructor(explorer),
      m_action_detector(lpsspec, m_trace_constructor, options.trace_actions, options.trace_multiactions, options.trace_prefix, options.max_traces),
      m_deadlock_detector(m_trace_constructor, options.trace_prefix, options.max_traces),
      m_nondeterminism_detector(m_trace_constructor, options.trace_prefix, options.max_traces),
      m_progress_monitor(options.search_strategy)
  {
    m_timed = lpsspec.process().has_time();
    if (options.detect_divergence)
    {
      m_divergence_detector = std::unique_ptr<detail::divergence_detector>(new detail::divergence_detector(explorer, m_timed, options.actions_internal_for_divergencies, options.trace_prefix, options.max_traces));
    }
  }

  // Explore the specification passed via the constructor, and put the results in builder.
  void explore(lts_builder& builder)
  {
    bool has_outgoing_transitions;
    const lps::state* source = nullptr;

    try
    {
      explorer.generate_state_space(
        m_timed,
        false,

        // discover_state
        [&](const lps::state& s, std::size_t s_index)
        {
          if (options.generate_traces && source)
          {
            m_trace_constructor.add_edge(*source, s);
          }
          if (options.detect_divergence)
          {
            m_divergence_detector->detect_divergence(s, s_index, m_trace_constructor);
          }
        },

        // examine_transition
        [&](const lps::state& s0, std::size_t s0_index, const process::timed_multi_action& a, const lps::state& s1, std::size_t s1_index, std::size_t summand_index)
        {
          builder.add_transition(s0_index, a, s1_index);
          has_outgoing_transitions = true;
          if (options.detect_action)
          {
            m_action_detector.detect_action(s0, s0_index, lps::multi_action(a.actions(), a.time()), s1, summand_index);
          }
          if (options.detect_nondeterminism)
          {
            m_nondeterminism_detector.detect_nondeterminism(s0, s0_index, lps::multi_action(a.actions(), a.time()), s1);
          }
          if (!options.suppress_progress_messages)
          {
            m_progress_monitor.examine_transition();
          }
        },

        // start_state
        [&](const lps::state& s, std::size_t /* s_index */)
        {
          source = &s;
          has_outgoing_transitions = false;
          if (options.detect_nondeterminism)
          {
            m_nondeterminism_detector.start_state();
          }
        },

        // finish_state
        [&](const lps::state& s, std::size_t s_index, std::size_t todo_list_size)
        {
          if (options.detect_deadlock && !has_outgoing_transitions)
          {
            m_deadlock_detector.detect_deadlock(s, s_index);
          }
          if (!options.suppress_progress_messages)
          {
            m_progress_monitor.finish_state(explorer.state_map().size(), todo_list_size);
          }
        }
      );
      m_progress_monitor.finish_exploration(explorer.state_map().size());
      builder.finalize(explorer.state_map());
    }
    catch (const lps::enumerator_error& e)
    {
      mCRL2log(log::error) << "Error while exploring state space: " << e.what() << "\n";
      if (options.save_error_trace)
      {
        const lps::state& s = *source;
        trace::Trace tr = m_trace_constructor.construct_trace(s);
        std::string filename = options.trace_prefix + "_error.trc";
        detail::save_trace(tr, filename);
      }
      mCRL2log(log::info) << ".\n";
    }
  }
};

struct stochastic_state_space_generator
{
  const lps::explorer_options& options;
  lps::explorer explorer;
  detail::trace_constructor m_trace_constructor;

  detail::action_detector m_action_detector;
  detail::deadlock_detector m_deadlock_detector;
  detail::nondeterminism_detector m_nondeterminism_detector;
  std::unique_ptr<detail::divergence_detector> m_divergence_detector;
  detail::progress_monitor m_progress_monitor;
  bool m_timed = false;

  stochastic_state_space_generator(const lps::stochastic_specification& lpsspec, const lps::explorer_options& options_)
    : options(options_),
      explorer(lpsspec, options_),
      m_trace_constructor(explorer),
      m_action_detector(lpsspec, m_trace_constructor, options.trace_actions, options.trace_multiactions, options.trace_prefix, options.max_traces),
      m_deadlock_detector(m_trace_constructor, options.trace_prefix, options.max_traces),
      m_nondeterminism_detector(m_trace_constructor, options.trace_prefix, options.max_traces),
      m_progress_monitor(options.search_strategy)
  {
    if (options.detect_divergence)
    {
      m_divergence_detector = std::unique_ptr<detail::divergence_detector>(new detail::divergence_detector(explorer, m_timed, options.actions_internal_for_divergencies, options.trace_prefix, options.max_traces));
    }
  }

  // Explore the specification passed via the constructor, and put the results in builder.
  void explore(stochastic_lts_builder& builder)
  {
    bool has_outgoing_transitions;
    const lps::state* source = nullptr;

    try
    {
      explorer.generate_stochastic_state_space(
        false,

        // discover_state
        [&](const lps::state& s, std::size_t s_index)
        {
          if (options.generate_traces && source)
          {
            m_trace_constructor.add_edge(*source, s);
          }
          if (options.detect_divergence)
          {
            m_divergence_detector->detect_divergence(s, s_index, m_trace_constructor);
          }
        },

        // examine_transition
        [&](const lps::state& s0, std::size_t s0_index, const process::timed_multi_action& a, const lps::stochastic_state& s1, const std::list<std::size_t>& s1_index, std::size_t summand_index)
        {
          builder.add_transition(s0_index, a, s1_index, s1.probabilities);
          has_outgoing_transitions = true;
          if (options.detect_action)
          {
            m_action_detector.detect_action(s0, s0_index, lps::multi_action(a.actions(), a.time()), s1.states.front(), summand_index);
          }
          if (options.detect_nondeterminism)
          {
            m_nondeterminism_detector.detect_nondeterminism(s0, s0_index, lps::multi_action(a.actions(), a.time()), s1.states.front());
          }
          if (!options.suppress_progress_messages)
          {
            m_progress_monitor.examine_transition();
          }
        },

        // start_state
        [&](const lps::state& s, std::size_t /* s_index */)
        {
          source = &s;
          has_outgoing_transitions = false;
          if (options.detect_nondeterminism)
          {
            m_nondeterminism_detector.start_state();
          }
        },

        // finish_state
        [&](const lps::state& s, std::size_t s_index, std::size_t todo_list_size)
        {
          if (options.detect_deadlock && !has_outgoing_transitions)
          {
            m_deadlock_detector.detect_deadlock(s, s_index);
          }
          if (!options.suppress_progress_messages)
          {
            m_progress_monitor.finish_state(explorer.state_map().size(), todo_list_size);
          }
        },

        // discover_initial_state
        [&](const lps::stochastic_state& s, const std::list<std::size_t>& s_index)
        {
          builder.set_initial_state(s_index, s.probabilities);
        }
      );
      m_progress_monitor.finish_exploration(explorer.state_map().size());
      builder.finalize(explorer.state_map());
    }
    catch (const lps::enumerator_error& e)
    {
      mCRL2log(log::error) << "Error while exploring state space: " << e.what() << "\n";
      if (options.save_error_trace)
      {
        const lps::state& s = *source;
        trace::Trace tr = m_trace_constructor.construct_trace(s);
        std::string filename = options.trace_prefix + "_error.trc";
        detail::save_trace(tr, filename);
      }
      mCRL2log(log::info) << ".\n";
    }
  }
};

} // namespace lts

} // namespace mcrl2

#endif // MCRL2_LTS_STATE_SPACE_GENERATOR_H