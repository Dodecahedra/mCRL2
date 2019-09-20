// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/confluence.h
/// \brief add your file description here.

#ifndef MCRL2_LPS_CONFLUENCE_H
#define MCRL2_LPS_CONFLUENCE_H

#include <memory>
#include "mcrl2/data/rewriter.h"
#include "mcrl2/data/consistency.h"
#include "mcrl2/data/join.h"
#include "mcrl2/data/rewriters/one_point_rule_rewriter.h"
#include "mcrl2/data/rewriters/quantifiers_inside_rewriter.h"
#include "mcrl2/data/standard.h" // equal_to
#include "mcrl2/data/substitution_utility.h"
#include "mcrl2/lps/stochastic_specification.h"
#include "mcrl2/smt/solver.h"

namespace mcrl2 {

namespace lps {

namespace detail {

inline
data::data_expression make_and(const data::data_expression& x1, const data::data_expression& x2, const data::data_expression& x3)
{
  return and_(x1, data::and_(x2, x3));
}

inline
data::data_expression make_and(const data::data_expression& x1, const data::data_expression& x2, const data::data_expression& x3, const data::data_expression& x4)
{
  return and_(data::and_(x1, x2), data::and_(x3, x4));
}

inline
data::data_expression equal_to(const data::data_expression_list& x, const data::data_expression_list& y)
{
  std::vector<data::data_expression> conjuncts;
  auto xi = x.begin();
  auto yi = y.begin();
  for (; xi != x.end(); ++xi, ++yi)
  {
    conjuncts.push_back(data::equal_to(*xi, *yi));
  }
  return data::join_and(conjuncts.begin(), conjuncts.end());
}

inline
std::pair<std::size_t, std::size_t> make_sorted_pair(std::size_t i, std::size_t j)
{
  return (i < j) ? std::make_pair(i, j) : std::make_pair(j, i);
}

inline
data::data_expression make_forall(const data::data_expression& x)
{
  std::set<data::variable> freevars = data::find_free_variables(x);
  return data::make_forall(data::variable_list(freevars.begin(), freevars.end()), x);
}

} // namespace detail

// TODO: reuse this code
template <typename Summand>
const stochastic_distribution& summand_distribution(const Summand& /* summand */)
{
  static stochastic_distribution empty_distribution;
  return empty_distribution;
}

// TODO: reuse this code
template <>
const stochastic_distribution& summand_distribution(const lps::stochastic_action_summand& summand)
{
  return summand.distribution();
}

inline
std::set<data::variable> used_read_variables(const action_summand& summand)
{
  std::set<data::variable> result;
  data::find_all_variables(summand.condition(), std::inserter(result, result.end()));
  data::find_all_variables(summand.multi_action(), std::inserter(result, result.end()));
  for (const auto& a : summand.assignments())
  {
    if (a.lhs() != a.rhs())
    {
      data::find_all_variables(a.rhs(), std::inserter(result, result.end()));
    }
  }
  return result;
}

inline
std::set<data::variable> changed_variables(const action_summand& summand)
{
  std::set<data::variable> result;
  for (const auto& a: summand.assignments())
  {
    result.insert(a.lhs());
  }
  return result;
}

// This is an alternative version of data::make_assignment_list that removes trivial assignments
inline
data::assignment_list make_assignment_list(const data::variable_list& variables, const data::data_expression_list& expressions)
{
  assert(variables.size() == expressions.size());
  std::vector<data::assignment> result;
  auto vi = variables.begin();
  auto ei = expressions.begin();
  for (; vi != variables.end(); ++vi, ++ei)
  {
    if (*vi != *ei)
    {
      result.emplace_back(*vi, *ei);
    }
  }
  return data::assignment_list(result.begin(), result.end());
}

struct confluence_summand
{
  data::variable_list variables;
  data::data_expression condition;
  lps::multi_action multi_action;
  stochastic_distribution distribution;
  data::data_expression_list next_state;
  std::set<data::variable> changed;
  std::set<data::variable> used;

  template <typename ActionSummand>
  confluence_summand(const ActionSummand& summand, const data::variable_list& process_parameters)
    : variables(summand.summation_variables()),
      condition(summand.condition()),
      multi_action(summand.multi_action().actions(), summand.multi_action().time()),
      distribution(summand_distribution(summand)),
      next_state(summand.next_state(process_parameters)),
      changed(changed_variables(summand)),
      used(used_read_variables(summand))
  {}

  bool is_tau() const
  {
    return multi_action.actions().empty();
  }
};

inline
std::string print_confluence_summand(const confluence_summand& summand, const data::variable_list& process_parameters)
{
  std::ostringstream out;
  out << lps::stochastic_action_summand(
    summand.variables,
    summand.condition,
    summand.multi_action,
    make_assignment_list(process_parameters, summand.next_state),
    summand.distribution
  );
  return out.str();
}

/// \brief Indicates whether or not two summands are disjoint.
inline
bool disjoint(const confluence_summand& summand1, const confluence_summand& summand2)
{
  using utilities::detail::has_empty_intersection;
  return has_empty_intersection(summand1.used, summand2.changed)
         && has_empty_intersection(summand2.used, summand1.changed)
         && has_empty_intersection(summand1.changed, summand2.changed);
}

/// \brief Creates an identifier for the ctau action
inline process::action_label make_ctau_act_id()
{
  static atermpp::aterm_appl ctau_act_id = atermpp::aterm_appl(core::detail::function_symbol_ActId(), atermpp::aterm_appl(atermpp::function_symbol("ctau", 0)), atermpp::aterm_list());
  return process::action_label(ctau_act_id);
}

/// \brief Creates the ctau action
inline process::action make_ctau_action()
{
  static atermpp::aterm_appl ctau_action = atermpp::aterm_appl(core::detail::function_symbol_Action(), make_ctau_act_id(), atermpp::aterm_list());
  return process::action(ctau_action);
}

template <typename Specification>
bool has_ctau_action(const Specification& lpsspec)
{
  using utilities::detail::contains;
  return contains(lpsspec.action_labels(), make_ctau_act_id());
}

/// \brief Function object that computes the condition for square confluence
/// \param summand_i A tau summand
/// \param summand_j An arbitrary action summand
struct square_confluence_condition
{
  const data::variable_list& process_parameters;
  data::mutable_indexed_substitution<>& sigma;

  square_confluence_condition(const data::variable_list& process_parameters_, data::mutable_indexed_substitution<>& sigma_)
   : process_parameters(process_parameters_), sigma(sigma_)
  {}

  data::data_expression operator()(const confluence_summand& summand_i, const confluence_summand& summand_j) const
  {
    assert(summand_i.is_tau());

    data::data_expression result;

    const auto& d = process_parameters;

    const auto& ci = summand_i.condition;
    const auto& gi = summand_i.next_state;

    const auto& cj = summand_j.condition;
    const auto& gj = summand_j.next_state;

    data::add_assignments(sigma, d, gi);
    data::data_expression cj_gi = data::replace_variables_capture_avoiding(cj, sigma);
    data::data_expression_list gj_gi = data::replace_variables_capture_avoiding(gj, sigma);

    data::add_assignments(sigma, d, gj);
    data::data_expression ci_gj = data::replace_variables_capture_avoiding(ci, sigma);
    data::data_expression_list gi_gj = data::replace_variables_capture_avoiding(gi, sigma);

    if (summand_j.is_tau())
    {
      data::remove_assignments(sigma, d);
      result = imp(data::and_(ci, cj), or_(detail::equal_to(gi, gj), detail::make_and(ci_gj, cj_gi, detail::equal_to(gj_gi, gi_gj))));
    }
    else
    {
      const auto& fj = summand_j.multi_action.arguments();
      data::data_expression_list fj_gi = data::replace_variables_capture_avoiding(fj, sigma);
      data::remove_assignments(sigma, d);
      result = imp(data::and_(ci, cj), detail::make_and(ci_gj, cj_gi, detail::equal_to(fj, fj_gi), detail::equal_to(gj_gi, gi_gj)));
    }

    return detail::make_forall(result);
  }
};

/// \brief Function object that computes the condition for triangular confluence
/// \param summand_i A tau summand
/// \param summand_j An arbitrary action summand
struct triangular_confluence_condition
{
  const data::variable_list& process_parameters;
  data::mutable_indexed_substitution<>& sigma;

  triangular_confluence_condition(const data::variable_list& process_parameters_, data::mutable_indexed_substitution<>& sigma_)
    : process_parameters(process_parameters_), sigma(sigma_)
  {}

  inline
  data::data_expression operator()(const confluence_summand& summand_i, const confluence_summand& summand_j) const
  {
    assert(summand_i.is_tau());

    data::data_expression result;
    
    const auto& d = process_parameters;

    const auto& ci = summand_i.condition;
    const auto& gi = summand_i.next_state;

    const auto& cj = summand_j.condition;
    const auto& gj = summand_j.next_state;

    data::add_assignments(sigma, d, gi);
    data::data_expression cj_gi = data::replace_variables_capture_avoiding(cj, sigma);
    data::data_expression_list gj_gi = data::replace_variables_capture_avoiding(gj, sigma);

    if (summand_j.is_tau())
    {
      data::remove_assignments(sigma, d);
      result = imp(and_(ci, cj), and_(cj_gi, detail::equal_to(gj_gi, gj)));
    }
    else
    {
      const auto& fj = summand_j.multi_action.arguments();
      data::data_expression_list fj_gi = data::replace_variables_capture_avoiding(fj, sigma);
      data::remove_assignments(sigma, d);
      result = imp(and_(ci, cj), detail::make_and(cj_gi, detail::equal_to(fj, fj_gi), detail::equal_to(gj_gi, gj)));
    }
    
    return detail::make_forall(result);
  }
};

/// \brief Function object that computes the condition for triangular confluence
/// \param summand_i A tau summand
/// \param summand_j An arbitrary action summand
struct trivial_confluence_condition
{
  data::mutable_indexed_substitution<>& sigma;

  explicit trivial_confluence_condition(data::mutable_indexed_substitution<>& sigma_)
    : sigma(sigma_)
  {}

  inline
  data::data_expression operator()(const confluence_summand& summand_i, const confluence_summand& summand_j) const
  {
    assert(summand_i.is_tau());

    data::data_expression result;

    data::data_expression aj_is_tau = summand_j.is_tau() ? data::true_() : data::false_();

    const auto& ci = summand_i.condition;
    const auto& gi = summand_i.next_state;

    const auto& cj = summand_j.condition;
    const auto& gj = summand_j.next_state;

    result = imp(and_(ci, cj), and_(aj_is_tau, detail::equal_to(gi, gj)));

    return detail::make_forall(result);
  }
};

class confluence_checker
{
  protected:
    std::vector<confluence_summand> m_summands;
    data::variable_list m_process_parameters;
    data::rewriter m_rewr;
    mutable data::mutable_indexed_substitution<> m_sigma;
    std::unique_ptr<smt::smt_solver> m_solver;

    // cache for the value of is_confluent for pairs (i, j) with i <= j, and i and j both tau-summands
    mutable std::map<std::pair<std::size_t, std::size_t>, bool> m_cache;

    enum cache_result
    {
      yes, no, indeterminate
    };

    cache_result cache_lookup(std::size_t i, std::size_t j) const
    {
      // check if the combination (i, j) is already in the cache
      auto found = m_cache.find(detail::make_sorted_pair(i, j));
      if (found != m_cache.end())
      {
        return found->second ? yes : no;
      }
      return indeterminate;
    }

    void cache_store(std::size_t i, std::size_t j, bool confluent) const
    {
      m_cache[detail::make_sorted_pair(i, j)] = confluent;
    }

    // check if x is a tautology using the data rewriter
    bool is_true_rewriter(data::data_expression x) const
    {
      data::one_point_rule_rewriter R_one_point;
      data::quantifiers_inside_rewriter R_quantifiers_inside;
      x = m_rewr(R_one_point(R_quantifiers_inside(x)));
      return data::is_true(x);
    }

    // check if x is a tautology using an smt solver
    bool is_true_smt(const data::data_expression& x) const
    {
      if (data::is_forall(x))
      {
        const auto& x_ = atermpp::down_cast<data::forall>(x);
        bool result = true;
        switch (m_solver->solve(x_.variables(), data::sort_bool::not_(x_.body())))
        {
          case smt::answer::SAT: result = true; break;
          case smt::answer::UNSAT: result = false; break;
            // since the formula is negated, we over-approximate unknown results
            // the result of this function will then be an under-approximation
          case smt::answer::UNKNOWN: result = true; break;
        }
        return !result;
      }
      else
      {
        // x has no free variables, so just evaluate the expression
        return data::is_true(m_rewr(x));
      }
    }

    // check if x evaluates true
    bool is_true(const data::data_expression& x) const
    {
      return m_solver ? is_true_smt(x) : is_true_rewriter(x);
    }

    // Returns whether the tau summand with index j is confluent. If not, the second value returned is
    // the index of the summand for which a violation was detected.
    template <typename ConfluenceCondition>
    std::pair<bool, std::size_t> is_confluent(std::size_t j, ConfluenceCondition confluence_condition, bool check_disjointness) const
    {
      const confluence_summand& summand_j = m_summands[j];
      for (std::size_t i = 0; i < m_summands.size(); i++)
      {
        const confluence_summand& summand_i =  m_summands[i];

        // check if the value for (i, j) is already in the cache
        if (summand_j.is_tau())
        {
          auto value = cache_lookup(i, j);
          if (value == yes)
          {
            mCRL2log(log::info) << '.';
            continue;
          }
          else if (value == no)
          {
            return { false, i };
          }
        }

        if (check_disjointness && disjoint(summand_j, summand_i))
        {
          cache_store(i, j, true);
          mCRL2log(log::info) << ':';
          continue;
        }

        data::data_expression condition = confluence_condition(summand_j, summand_i);
        bool confluent = is_true(condition);
        cache_store(i, j, confluent);
        if (confluent)
        {
          mCRL2log(log::info) << '+';
        }
        else
        {
          return { false, i };
        }
      }
      return { true, 0 };
    }

  public:
    std::pair<bool, std::size_t> is_square_confluent(std::size_t j) const
    {
      bool check_disjointness = true;
      return is_confluent(j, square_confluence_condition(m_process_parameters, m_sigma), check_disjointness);
    }

    std::pair<bool, std::size_t> is_triangular_confluent(std::size_t j) const
    {
      bool check_disjointness = false;
      return is_confluent(j, triangular_confluence_condition(m_process_parameters, m_sigma), check_disjointness);
    }

    std::pair<bool, std::size_t> is_trivial_confluent(std::size_t j) const
    {
      bool check_disjointness = false;
      return is_confluent(j, trivial_confluence_condition(m_sigma), check_disjointness);
    }

    template <typename Specification>
    std::vector<std::size_t> compute_confluent_summands(const Specification& lpsspec, char confluence_type)
    {
      std::vector<std::size_t> result;

      m_cache.clear();
      m_summands.clear();
      m_process_parameters = lpsspec.process().process_parameters();
      m_rewr = data::rewriter(lpsspec.data());

      for (const auto& summand: lpsspec.process().action_summands())
      {
        m_summands.emplace_back(summand, m_process_parameters);
      }

      std::size_t n = m_summands.size();
      std::size_t tau_summand_count = 0;
      for (std::size_t j = 0; j < n; j++)
      {
        const auto& summand_j = m_summands[j];
        if (!summand_j.is_tau())
        {
          continue;
        }
        tau_summand_count++;
        mCRL2log(log::info) << "summand " << (j + 1) << " of " << n << " (condition = " << confluence_type << "): ";
        bool confluent;
        std::size_t violating_index;
        switch (confluence_type)
        {
          case 'C': std::tie(confluent, violating_index) = is_square_confluent(j); break;
          case 'T': std::tie(confluent, violating_index) = is_triangular_confluent(j); break;
          case 'Z': std::tie(confluent, violating_index) = is_trivial_confluent(j); break;
          default: throw mcrl2::runtime_error("Unknown confluence type " + std::to_string(confluence_type));
        }
        if (confluent)
        {
          result.push_back(j);
          mCRL2log(log::info) << "Confluent with all summands";
        }
        else
        {
          mCRL2log(log::info) << "Not confluent with summand " << (violating_index + 1);
        }
        mCRL2log(log::info) << std::endl;
      }
      mCRL2log(log::info) << result.size() << " of " << tau_summand_count << " tau summands were found to be confluent";
      return result;
    }

    /// \brief Applies confluent reduction to an LPS
    template <typename Specification>
    void run(Specification& lpsspec, char confluence_type, bool use_smt_solver = false)
    {
      if (has_ctau_action(lpsspec))
      {
        throw mcrl2::runtime_error("An action named \'ctau\' already exists.\n");
      }

      if (use_smt_solver)
      {
        m_solver = std::unique_ptr<smt::smt_solver>(new smt::smt_solver(lpsspec.data()));
      }

      multi_action ctau{make_ctau_action()};
      std::vector<std::size_t> I = compute_confluent_summands(lpsspec, confluence_type);
      auto& summands = lpsspec.process().action_summands();
      for (std::size_t i: I)
      {
        summands[i].multi_action() = ctau;
      }
    }
};

} // namespace lps

} // namespace mcrl2

#endif // MCRL2_LPS_CONFLUENCE_H
