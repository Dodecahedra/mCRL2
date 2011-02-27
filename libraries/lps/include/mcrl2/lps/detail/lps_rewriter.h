// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/detail/lps_rewriter.h
/// \brief Add your file description here.

#ifndef MCRL2_LPS_DETAIL_LPS_REWRITER_H
#define MCRL2_LPS_DETAIL_LPS_REWRITER_H

#include <vector>
// #include "mcrl2/data/assignment.h"
#include "mcrl2/lps/specification.h"

namespace mcrl2
{

namespace lps
{

namespace detail
{

/// \brief Function object for applying a data rewriter to LPS data types.
template <typename DataRewriter>
struct lps_rewriter
{
  const DataRewriter& R;

  lps_rewriter(const DataRewriter& R_)
    : R(R_)
  {}

  /// \brief Applies the rewriter to the elements of a term list
  template <typename TermList>
  TermList rewrite_list_copy(const TermList& l) const
  {
    // TODO: how to make this function efficient?
    typedef typename TermList::value_type value_type;
    atermpp::vector<value_type> v(l.begin(), l.end());
    for (typename std::vector<value_type>::iterator i = v.begin(); i != v.end(); ++i)
    {
      rewrite(*i);
    }
    return TermList(v.begin(), v.end());
  }

  /// \brief Applies the rewriter to the elements of an assignment list, 
  ///           where the element is removed if the rewritten rhs equals the lhs.
  
  data::assignment_list rewrite_list_copy(const data::assignment_list& l) const
  {
    using namespace data;
    // TODO: how to make this function efficient?
    atermpp::vector<assignment> v;
    for (assignment_list::const_iterator i = l.begin(); i != l.end(); ++i)
    {
      assignment a=*i;
      rewrite(a);
      if (a.lhs()!=a.rhs())
      { 
        v.push_back(a);
      }
    }
    return assignment_list(v.begin(), v.end());
  }
  
  /// \brief Applies the rewriter to the elements of a container
  /// \brief Applies the rewriter to the elements of a term list
  template <typename TermList>
  void rewrite_list(TermList& l) const
  {
    l = rewrite_list_copy(l);
  }

  template <typename Container>
  void rewrite_container(Container& c) const
  {
    for (typename Container::iterator i = c.begin(); i != c.end(); ++i)
    {
      rewrite(*i);
    }
  }

  /// \brief Applies the rewriter to a data expression
  /// \param d A data expression
  data::data_expression rewrite_copy(const data::data_expression& d) const
  {
    return R(d);
  }

  /// \brief Applies the rewriter to a data expression
  /// \param d A data expression
  void rewrite(data::data_expression& d) const
  {
    d = rewrite_copy(d);
  }

  /// \brief Applies the rewriter to an assignment
  /// \param a An assignment
  void rewrite(data::assignment& a) const
  {
    a = data::assignment(a.lhs(), rewrite_copy(a.rhs()));
  }

  /// \brief Applies the rewriter to an action
  /// \param a An action
  void rewrite(action& a) const
  {
    a = action(a.label(), rewrite_list_copy(a.arguments()));
  }

  /// \brief Applies the rewriter to a deadlock
  /// \param d A deadlock
  void rewrite(deadlock& d) const
  {
    if (d.has_time())
    {
      rewrite(d.time());
    }
  }

  /// \brief Applies the rewriter to a multi-action
  /// \param a A multi-action
  void rewrite(multi_action& a) const
  {
    if (a.has_time())
    {
      rewrite(a.time());
    }
    rewrite_list(a.actions());
  }

  /// \brief Applies the rewriter to a summand
  /// \param s A summand
  void rewrite(action_summand& s) const
  {
    rewrite(s.condition());
    rewrite(s.multi_action());
    rewrite_list(s.assignments());
  }

  /// \brief Applies the rewriter to a summand
  /// \param s A summand
  void rewrite(deadlock_summand& s) const
  {
    rewrite(s.condition());
    rewrite(s.deadlock());
  }

  /// \brief Applies the rewriter to a process_initializer
  /// \param s A process_initializer
  void rewrite(process_initializer& i) const
  {
    i = process_initializer(rewrite_list_copy(i.assignments()));
  }

  /// \brief Applies the rewriter to a linear_process
  /// \param s A linear_process
  void rewrite(linear_process& p) const
  {
    rewrite_container(p.action_summands());
    rewrite_container(p.deadlock_summands());
  }

  /// \brief Applies the rewriter to a linear process specification
  /// \param spec A linear process specification
  void rewrite(specification& spec) const
  {
    rewrite(spec.process());
    rewrite(spec.initial_process());
  }

  /// \brief Applies the rewriter to the elements of a term list
  template <typename TermList>
  void rewrite(TermList& l) const
  {
    l = rewrite_list_copy(l);
  }

  template <typename Term>
  void operator()(Term& t)
  {
    rewrite(t);
  }
};

/// \brief Utility function to create an lps_rewriter.
template <typename DataRewriter>
lps_rewriter<DataRewriter> make_lps_rewriter(const DataRewriter& R)
{
  return lps_rewriter<DataRewriter>(R);
}

} // namespace detail

} // namespace lps

} // namespace mcrl2

#endif // MCRL2_LPS_DETAIL_LPS_REWRITER_H
