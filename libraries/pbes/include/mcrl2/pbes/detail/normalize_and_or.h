// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/detail/normalize_and_or.h
/// \brief Function to normalize 'and' and 'or' sub expressions.

#ifndef MCRL2_PBES_DETAIL_NORMALIZE_AND_OR_H
#define MCRL2_PBES_DETAIL_NORMALIZE_AND_OR_H

#include <set>
#include <utility>
#include "mcrl2/utilities/optimized_boolean_operators.h"
#include "mcrl2/pbes/builder.h"
#include "mcrl2/pbes/pbes_expression.h"

namespace mcrl2
{

namespace pbes_system
{

namespace detail
{

// Simplifying PBES rewriter.
template <typename Derived>
struct normalize_and_or_builder: public pbes_expression_builder<Derived>
{
  typedef pbes_expression_builder<Derived> super;
  using super::enter;
  using super::leave;
  using super::operator();

  /// \brief Splits a disjunction into a sequence of operands
  /// Given a pbes expression of the form p1 || p2 || .... || pn, this will yield a
  /// set of the form { p1, p2, ..., pn }, assuming that pi does not have a || as main
  /// function symbol.
  /// \param expr A PBES expression
  /// \return A sequence of operands
  atermpp::multiset<pbes_expression> split_or(const pbes_expression& expr)
  {
    using namespace accessors;
    atermpp::multiset<pbes_expression> result;
    core::detail::split(expr, std::insert_iterator<atermpp::multiset<pbes_expression> >(result, result.begin()), is_or, left, right);
    return result;
  }

  /// \brief Splits a conjunction into a sequence of operands
  /// Given a pbes expression of the form p1 && p2 && .... && pn, this will yield a
  /// set of the form { p1, p2, ..., pn }, assuming that pi does not have a && as main
  /// function symbol.
  /// \param expr A PBES expression
  /// \return A sequence of operands
  atermpp::multiset<pbes_expression> split_and(const pbes_expression& expr)
  {
    using namespace accessors;
    atermpp::multiset<pbes_expression> result;
    core::detail::split(expr, std::insert_iterator<atermpp::multiset<pbes_expression> >(result, result.begin()), is_and, left, right);
    return result;
  }

  pbes_expression normalize(const pbes_expression& x)
  {
    if (is_and(x))
    {
      atermpp::multiset<pbes_expression> s = split_and(x);
      return pbes_expr::join_and(s.begin(), s.end());
    }
    else if (is_or(x))
    {
      atermpp::multiset<pbes_expression> s = split_or(x);
      return pbes_expr::join_or(s.begin(), s.end());
    }
    return x;
  }

  // to prevent default operator() being called
  data::data_expression operator()(const data::data_expression& x)
  {
  	return x;
  }

  pbes_expression operator()(const pbes_expression& x)
  {
    return normalize(super::operator()(x));
  }
};

template <typename T>
T normalize_and_or(const T& x,
                   typename boost::enable_if<typename boost::is_base_of<atermpp::aterm_base, T>::type>::type* = 0
                  )
{
  return core::make_apply_builder<normalize_and_or_builder>()(x);
}

template <typename T>
void normalize_and_or(T& x,
                      typename boost::disable_if<typename boost::is_base_of<atermpp::aterm_base, T>::type>::type* = 0
                     )
{
  core::make_apply_builder<normalize_and_or_builder>()(x);
}

} // namespace detail

} // namespace pbes_system

} // namespace mcrl2

#endif // MCRL2_PBES_DETAIL_NORMALIZE_AND_OR_H
