/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


#ifndef CASADI_GENERIC_EXPRESSION_HPP
#define CASADI_GENERIC_EXPRESSION_HPP

#include "calculus.hpp"

namespace casadi {

  /** \brief Empty Base

      This class is extended in SWIG.

      \identifier{ol} */
  struct CASADI_EXPORT GenericExpressionCommon {};

#ifndef SWIG
  /** \brief Expression interface
  *
  This is a common base class for SX, MX and Matrix<>, introducing a uniform syntax and implementing
  common functionality using the curiously recurring template pattern (CRTP) idiom.\n

  \author Joel Andersson
  \date 2012

      \identifier{om} */
template<typename ExType>
class GenericExpression : public GenericExpressionCommon {
  protected:
    // Helper functions
    inline const ExType& self() const { return static_cast<const ExType&>(*this); }
    inline ExType& self() { return static_cast<ExType&>(*this); }
  public:

/**
\addtogroup expression_tools
@{
*/

  ///@{
  /** \brief Addition: (x,y) -> x + y

      Elementwise sum. The shape rules below hold for every binary operation in this group
      (minus, times, rdivide, the comparisons, logic_and, ...):
      - both operands have the same shape, or
      - one operand is a scalar (1-by-1), which is combined with every element, or
      - the operands have the same number of rows and the column count of one is a multiple
        of the other's; the narrower one is repeated horizontally (so a column vector can be
        added to every column of a matrix, but a row vector cannot be added to every row).
      Any other combination raises an error; there is no further broadcasting. Empty
      matrices follow the same rules and give an empty result.

      The result is sparse where both operands are structurally zero; an entry that is a
      structural zero in only one operand becomes an explicit value. Operations that map
      (0, y) or (x, 0) to a nonzero, such as x + 1 or cos, densify the result. Numeric (DM)
      operands are evaluated, symbolic ones (SX, MX) build an expression; DM and SX or MX can
      be mixed, SX and MX cannot.

      \lang{python}
      \doctest
      >>> plus(DM([1, 2]), DM([10, 20]))
      DM([11, 22])
      >>> DM([1, 2]) + 1
      DM([2, 3])
      >>> DM.eye(2) + DM.eye(2)
      DM(
      [[2, 00],
       [00, 2]])
      >>> DM.eye(2) + 1
      DM(
      [[2, 1],
       [1, 2]])
      >>> SX.sym("x", 2, 2) + SX.sym("y", 2)
      SX(
      [[(x_0+y_0), (x_2+y_0)],
       [(x_1+y_1), (x_3+y_1)]])
      >>> (MX.sym("x", 2, 0) + 1).shape
      (2, 0)
      >>> DM.ones(2, 2) + DM.ones(1, 2)
      Traceback (most recent call last):
        ...
      RuntimeError: ...Dimension mismatch...
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> plus(DM([1; 2]), DM([10; 20]))
      [11, 22]
      >> DM([1; 2]) + 1
      [2, 3]
      >> DM.eye(2) + DM.eye(2)
      [[2, 00],
       [00, 2]]
      >> DM.eye(2) + 1
      [[2, 1],
       [1, 2]]
      >> SX.sym('x', 2, 2) + SX.sym('y', 2)
      [[(x_0+y_0), (x_2+y_0)],
       [(x_1+y_1), (x_3+y_1)]]
      >> size(MX.sym('x', 2, 0) + 1)
      2 0
      >> try, DM.ones(2, 2) + DM.ones(1, 2); catch e, disp(e.message); end
      ...Dimension mismatch...
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.plus(ca.DM([1, 2]), ca.DM([10, 20]))
      DM([11, 22])
      julia> ca.DM([1, 2]) + 1
      DM([2, 3])
      julia> ca.DM.eye(2) + ca.DM.eye(2)
      DM(
      [[2, 00],
       [00, 2]])
      julia> ca.DM.eye(2) + 1
      DM(
      [[2, 1],
       [1, 2]])
      julia> ca.SX.sym("x", 2, 2) + ca.SX.sym("y", 2)
      SX(
      [[(x_0+y_0), (x_2+y_0)],
       [(x_1+y_1), (x_3+y_1)]])
      julia> size(ca.MX.sym("x", 2, 0) + 1)
      (2, 0)
      julia> try
                 ca.DM.ones(2, 2) + ca.DM.ones(1, 2)
             catch e
                 println(occursin("Dimension mismatch", sprint(showerror, e)))
             end
      true
      \enddoctest
      \endlang

      \identifier{on} */
  static ExType plus(const ExType &x, const ExType &y) {
    return ExType::binary(OP_ADD, x, y);
  }
  friend inline ExType plus(const ExType &x, const ExType &y) {
    return ExType::plus(x, y);
  }
  friend inline ExType operator+(const ExType &x, const ExType &y) {
    return plus(x, y);
  }
  inline ExType& operator+=(const ExType &y) { return self() = self() + y; }
  ///@}

  ///@{
  /** \brief Subtraction: (x,y) -> x - y

      Elementwise difference; shape and sparsity rules as for plus. For symbolic arguments,
      x - x with the same node simplifies to 0.

      \lang{python}
      \doctest
      >>> minus(DM([5, 6]), DM([1, 2]))
      DM([4, 4])
      >>> DM.eye(2) - 1
      DM(
      [[0, -1],
       [-1, 0]])
      >>> x = SX.sym("x")
      >>> x - x
      SX(0)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> minus(DM([5; 6]), DM([1; 2]))
      [4, 4]
      >> DM.eye(2) - 1
      [[0, -1],
       [-1, 0]]
      >> x = SX.sym('x');
      >> x - x
      0
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.minus(ca.DM([5, 6]), ca.DM([1, 2]))
      DM([4, 4])
      julia> ca.DM.eye(2) - 1
      DM(
      [[0, -1],
       [-1, 0]])
      julia> x = ca.SX.sym("x");
      julia> x - x
      SX(0)
      \enddoctest
      \endlang

      \identifier{oo} */
  static ExType minus(const ExType &x, const ExType &y) {
    return ExType::binary(OP_SUB, x, y);
  }
  friend inline ExType minus(const ExType &x, const ExType &y) {
    return ExType::minus(x, y);
  }
  friend inline ExType operator-(const ExType &x, const ExType &y) {
    return minus(x, y);
  }
  inline ExType& operator-=(const ExType &y) { return self() = self() - y; }
  ///@}

  ///@{
  /** \brief Elementwise multiplication: (x,y) -> x .* y

      Elementwise (Hadamard) product, not the matrix product; use mtimes for that. Shape rules
      as for plus. The result is structurally zero wherever either operand is, so the sparsity
      pattern is the intersection of the patterns.

      \lang{python}
      \doctest
      >>> times(DM([1, 2]), DM([3, 4]))
      DM([3, 8])
      >>> DM([[1, 2], [3, 4]]) * DM([[1, 2], [3, 4]])
      DM(
      [[1, 4],
       [9, 16]])
      >>> DM.eye(2) * DM.ones(2, 2)
      DM(
      [[1, 00],
       [00, 1]])
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> times(DM([1; 2]), DM([3; 4]))
      [3, 8]
      >> DM([1 2; 3 4]) .* DM([1 2; 3 4])
      [[1, 4],
       [9, 16]]
      >> DM([1 2; 3 4]) * DM([1 2; 3 4])
      [[7, 10],
       [15, 22]]
      >> DM.eye(2) .* DM.ones(2, 2)
      [[1, 00],
       [00, 1]]
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.times(ca.DM([1, 2]), ca.DM([3, 4]))
      DM([3, 8])
      julia> ca.DM([1 2; 3 4]) .* ca.DM([1 2; 3 4])
      DM(
      [[1, 4],
       [9, 16]])
      julia> ca.DM.eye(2) .* ca.DM.ones(2, 2)
      DM(
      [[1, 00],
       [00, 1]])
      \enddoctest
      \endlang

      \sa mtimes

      \identifier{op} */
  static ExType times(const ExType &x, const ExType &y) {
    return ExType::binary(OP_MUL, x, y);
  }
  friend inline ExType times(const ExType &x, const ExType &y) {
    return ExType::times(x, y);
  }
  friend inline ExType operator*(const ExType &x, const ExType &y) {
    return times(x, y);
  }
  inline ExType& operator*=(const ExType &y) {return self() = self() * y;}
  ///@}

  ///@{
  /** \brief Elementwise division: (x,y) -> x ./ y

      Elementwise quotient, not a linear solve; use solve or mldivide for that. Shape rules as
      for plus. Division by zero follows IEEE arithmetic (inf or nan), no error is raised.
      Structural zeros of \p x stay structural zeros; dividing by a sparse \p y densifies.

      \lang{python}
      \doctest
      >>> rdivide(DM([1, 2]), 4)
      DM([0.25, 0.5])
      >>> DM([1, 2]) / DM([4, 8])
      DM([0.25, 0.25])
      >>> 1 / DM.eye(2)
      DM(
      [[1, inf],
       [inf, 1]])
      >>> DM.eye(2) / 2
      DM(
      [[0.5, 00],
       [00, 0.5]])
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> rdivide(DM([1; 2]), 4)
      [0.25, 0.5]
      >> DM([1; 2]) ./ DM([4; 8])
      [0.25, 0.25]
      >> 1 ./ DM.eye(2)
      [[1, inf],
       [inf, 1]]
      >> DM.eye(2) ./ 2
      [[0.5, 00],
       [00, 0.5]]
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.rdivide(ca.DM([1, 2]), 4)
      DM([0.25, 0.5])
      julia> ca.DM([1, 2]) ./ ca.DM([4, 8])
      DM([0.25, 0.25])
      julia> 1 ./ ca.DM.eye(2)
      DM(
      [[1, inf],
       [inf, 1]])
      julia> ca.DM.eye(2) / 2
      DM(
      [[0.5, 00],
       [00, 0.5]])
      \enddoctest
      \endlang

      \sa mldivide, mrdivide, solve

      \identifier{oq} */
  static ExType rdivide(const ExType &x, const ExType &y) {
    return ExType::binary(OP_DIV, x, y);
  }
  friend inline ExType rdivide(const ExType &x, const ExType &y) {
    return ExType::rdivide(x, y);
  }
  friend inline ExType operator/(const ExType &x, const ExType &y) {
    return rdivide(x, y);
  }
  inline ExType& operator/=(const ExType &y) {return self() = self() / y;}
  ///@}

  ///@{
  /** \brief Logical less than: (x,y) -> x < y

      Elementwise comparison; shape rules as for plus. Every comparison in this group returns
      a matrix of the same type as the operands with entries 1 (true) and 0 (false): for DM a
      numeric 0/1 matrix, not a boolean, and for SX and MX a symbolic expression that is only
      decided at evaluation time. Comparisons with NaN are false. The result is sparse where
      both operands are structurally zero (0 < 0 is false).

      \lang{python}
      \doctest
      >>> lt(DM([1, 2, 3]), 2)
      DM([1, 0, 0])
      >>> DM([1, 2, 3]) < DM([3, 2, 1])
      DM([1, 0, 0])
      >>> SX.sym("x", 2) < SX.sym("y", 2)
      SX([(x_0<y_0), (x_1<y_1)])
      >>> DM.eye(2) < 1
      DM(
      [[0, 1],
       [1, 0]])
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> lt(DM([1; 2; 3]), 2)
      [1, 0, 0]
      >> DM([1; 2; 3]) < DM([3; 2; 1])
      [1, 0, 0]
      >> SX.sym('x', 2) < SX.sym('y', 2)
      [(x_0<y_0), (x_1<y_1)]
      >> DM.eye(2) < 1
      [[0, 1],
       [1, 0]]
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.lt(ca.DM([1, 2, 3]), 2)
      DM([1, 0, 0])
      julia> ca.DM([1, 2, 3]) < ca.DM([3, 2, 1])
      DM([1, 0, 0])
      julia> ca.SX.sym("x", 2) < ca.SX.sym("y", 2)
      SX([(x_0<y_0), (x_1<y_1)])
      julia> ca.DM.eye(2) < 1
      DM(
      [[0, 1],
       [1, 0]])
      \enddoctest
      \endlang

      \sa le, gt, ge, eq, ne, if_else

      \identifier{or} */
  static ExType lt(const ExType &x, const ExType &y) {
    return ExType::binary(OP_LT, x, y);
  }
  friend inline ExType lt(const ExType &x, const ExType &y) {
    return ExType::lt(x, y);
  }
  friend inline ExType operator<(const ExType &x, const ExType &y) {
    return lt(x, y);
  }
  ///@}

  ///@{
  /** \brief Logical less or equal to: (x,y) -> x <= y

      Elementwise comparison with 0/1 result, see lt. Chained inequalities such as
      lb <= (x <= ub) are recognised by Opti as bound constraints.

      \lang{python}
      \doctest
      >>> le(DM([1, 2, 3]), 2)
      DM([1, 1, 0])
      >>> x = SX.sym("x")
      >>> x <= x
      SX(1)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> le(DM([1; 2; 3]), 2)
      [1, 1, 0]
      >> x = SX.sym('x');
      >> x <= x
      1
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.le(ca.DM([1, 2, 3]), 2)
      DM([1, 1, 0])
      julia> x = ca.SX.sym("x");
      julia> x <= x
      SX(1)
      \enddoctest
      \endlang

      \sa lt

      \identifier{os} */
  static ExType le(const ExType &x, const ExType &y) {
    return ExType::binary(OP_LE, x, y);
  }
  friend inline ExType le(const ExType &x, const ExType &y) {
    return ExType::le(x, y);
  }
  friend inline ExType operator<=(const ExType &x, const ExType &y) {
    return le(x, y);
  }
  ///@}

  ///@{
  /** \brief Logical greater than: (x,y) -> x > y

      Elementwise comparison with 0/1 result, see lt. Implemented as y < x, which is how a
      symbolic expression prints.

      \lang{python}
      \doctest
      >>> gt(DM([1, 2, 3]), 2)
      DM([0, 0, 1])
      >>> gt(SX.sym("x"), SX.sym("y"))
      SX((y<x))
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> gt(DM([1; 2; 3]), 2)
      [0, 0, 1]
      >> SX.sym('x') > SX.sym('y')
      (y<x)
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.gt(ca.DM([1, 2, 3]), 2)
      DM([0, 0, 1])
      julia> ca.gt(ca.SX.sym("x"), ca.SX.sym("y"))
      SX((y<x))
      \enddoctest
      \endlang

      \sa lt

      \identifier{ot} */
  static ExType gt(const ExType &x, const ExType &y) {
    return ExType::lt(y, x);
  }
  friend inline ExType gt(const ExType &x, const ExType &y) {
    return ExType::gt(x, y);
  }
  friend inline ExType operator>(const ExType &x, const ExType &y) {
    return gt(x, y);
  }
  ///@}

  ///@{
  /** \brief Logical greater or equal to: (x,y) -> x >= y

      Elementwise comparison with 0/1 result, see lt. Implemented as y <= x.

      \lang{python}
      \doctest
      >>> ge(DM([1, 2, 3]), 2)
      DM([0, 1, 1])
      >>> ge(MX.sym("x"), MX.sym("y"))
      MX((y<=x))
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> ge(DM([1; 2; 3]), 2)
      [0, 1, 1]
      >> MX.sym('x') >= MX.sym('y')
      (y<=x)
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.ge(ca.DM([1, 2, 3]), 2)
      DM([0, 1, 1])
      julia> ca.ge(ca.MX.sym("x"), ca.MX.sym("y"))
      MX((y<=x))
      \enddoctest
      \endlang

      \sa lt

      \identifier{ou} */
  static ExType ge(const ExType &x, const ExType &y) {
    return ExType::le(y, x);
  }
  friend inline ExType ge(const ExType &x, const ExType &y) {
    return ExType::ge(x, y);
  }
  friend inline ExType operator>=(const ExType &x, const ExType &y) {
    return ge(x, y);
  }
  ///@}

  ///@{
  /** \brief Logical equal to: (x,y) -> x == y

      Elementwise comparison with 0/1 result, see lt: this is value equality of the entries,
      not identity of expressions; use is_equal to test whether two symbolic expressions are
      the same. For symbolic operands x == x with the same node simplifies to 1. NaN is not
      equal to anything, including itself.

      \lang{python}
      \doctest
      >>> eq(DM([1, 2, 3]), 2)
      DM([0, 1, 0])
      >>> DM.eye(2) == DM([[1, 0], [0, 1]])
      DM(
      [[1, 1],
       [1, 1]])
      >>> x = SX.sym("x")
      >>> x == x, x == SX.sym("x")
      (SX(1), SX((x==x)))
      >>> DM(float("nan")) == DM(float("nan"))
      DM(0)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> eq(DM([1; 2; 3]), 2)
      [0, 1, 0]
      >> DM.eye(2) == DM([1 0; 0 1])
      [[1, 1],
       [1, 1]]
      >> x = SX.sym('x');
      >> x == x, x == SX.sym('x')
      1
      (x==x)
      >> DM(nan) == DM(nan)
      0
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.eq(ca.DM([1, 2, 3]), 2)
      DM([0, 1, 0])
      julia> ca.DM.eye(2) == ca.DM([1 0; 0 1])
      DM(
      [[1, 1],
       [1, 1]])
      julia> x = ca.SX.sym("x");
      julia> x == x, x == ca.SX.sym("x")
      (1, (x==x))
      julia> ca.DM(NaN) == ca.DM(NaN)
      DM(0)
      \enddoctest
      \endlang

      \sa ne, is_equal

      \identifier{ov} */
  static ExType eq(const ExType &x, const ExType &y) {
    return ExType::binary(OP_EQ, x, y);
  }
  friend inline ExType eq(const ExType &x, const ExType &y) {
    return ExType::eq(x, y);
  }
  friend inline ExType operator==(const ExType &x, const ExType &y) {
    return eq(x, y);
  }
  ///@}

  ///@{
  /** \brief Logical not equal to: (x,y) -> x != y

      Elementwise comparison with 0/1 result, see lt and eq. For symbolic operands x != x with
      the same node simplifies to 0, so it cannot be used to detect NaN symbolically.

      \lang{python}
      \doctest
      >>> ne(DM([1, 2, 3]), 2)
      DM([1, 0, 1])
      >>> x = MX.sym("x")
      >>> x != x
      MX(0)
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> ne(DM([1; 2; 3]), 2)
      [1, 0, 1]
      >> x = MX.sym('x');
      >> x ~= x
      0
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.ne(ca.DM([1, 2, 3]), 2)
      DM([1, 0, 1])
      julia> x = ca.MX.sym("x");
      julia> x != x
      MX(0)
      \enddoctest
      \endlang

      \sa eq

      \identifier{ow} */
  static ExType ne(const ExType &x, const ExType &y) {
    return ExType::binary(OP_NE, x, y);
  }
  friend inline ExType ne(const ExType &x, const ExType &y) {
    return ExType::ne(x, y);
  }
  friend inline ExType operator!=(const ExType &x, const ExType &y) {
    return ne(x, y);
  }
  ///@}

  ///@{
  /** \brief Logical `and`

   * Returns (an expression evaluating to) 1 where both
   * entries are nonzero and 0 otherwise, elementwise

      Any nonzero value, including NaN, counts as true. Shape rules as for plus; the result is
      structurally zero wherever either operand is. The Python keyword `and` does not work on
      matrices (it needs a single truth value); use this function.

      \lang{python}
      \doctest
      >>> logic_and(DM([1, 2, 0]), DM([3, 0, 0]))
      DM([1, 0, 0])
      >>> logic_and(SX.sym("x"), SX.sym("y"))
      SX((x&&y))
      >>> logic_and(DM.eye(2), DM.ones(2, 2))
      DM(
      [[1, 00],
       [00, 1]])
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> DM([1; 2; 0]) & DM([3; 0; 0])
      [1, 0, 0]
      >> SX.sym('x') & SX.sym('y')
      (x&&y)
      >> and(DM.eye(2), DM.ones(2, 2))
      [[1, 00],
       [00, 1]]
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.logic_and(ca.DM([1, 2, 0]), ca.DM([3, 0, 0]))
      DM([1, 0, 0])
      julia> ca.logic_and(ca.SX.sym("x"), ca.SX.sym("y"))
      SX((x&&y))
      julia> ca.logic_and(ca.DM.eye(2), ca.DM.ones(2, 2))
      DM(
      [[1, 00],
       [00, 1]])
      \enddoctest
      \endlang

      \sa logic_or, logic_not

      \identifier{ox} */
   static ExType logic_and(const ExType &x, const ExType &y) {
     return ExType::binary(OP_AND, x, y);
   }
   friend inline ExType logic_and(const ExType &x, const ExType &y) {
     return ExType::logic_and(x, y);
   }
   friend inline ExType operator&&(const ExType &x, const ExType &y) {
     return logic_and(x, y);
   }
  ///@}

  ///@{
  /** \brief  Logical `or`

   * returns (an expression evaluating to) 1 where at
   * least one entry is nonzero and 0 otherwise, elementwise

      Shape rules as for plus; the result is structurally zero where both operands are.

      \lang{python}
      \doctest
      >>> logic_or(DM([1, 0, 0]), DM([0, 0, 2]))
      DM([1, 0, 1])
      >>> logic_or(MX.sym("x"), MX.sym("y"))
      MX((x||y))
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> DM([1; 0; 0]) | DM([0; 0; 2])
      [1, 0, 1]
      >> or(MX.sym('x'), MX.sym('y'))
      (x||y)
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.logic_or(ca.DM([1, 0, 0]), ca.DM([0, 0, 2]))
      DM([1, 0, 1])
      julia> ca.logic_or(ca.MX.sym("x"), ca.MX.sym("y"))
      MX((x||y))
      \enddoctest
      \endlang

      \sa logic_and, logic_not

      \identifier{oy} */
   static ExType logic_or(const ExType &x, const ExType &y) {
     return ExType::binary(OP_OR, x, y);
   }
   friend inline ExType logic_or(const ExType &x, const ExType &y) {
     return ExType::logic_or(x, y);
   }
   friend inline ExType operator||(const ExType &x, const ExType &y) {
     return logic_or(x, y);
   }
   ///@}

   ///@{
   /** \brief  Logical `not` x -> !x

    * Returns (an expression evaluating to) 1 where the
    * entry is zero and 0 otherwise, elementwise

      Structural zeros map to 1, so a sparse input gives a dense result. The Python keyword
      `not` does not work on matrices; use this function.

      \lang{python}
      \doctest
      >>> logic_not(DM([1, 0, -1]))
      DM([0, 1, 0])
      >>> logic_not(DM.eye(2))
      DM(
      [[0, 1],
       [1, 0]])
      >>> logic_not(SX.sym("x"))
      SX((!x))
      \enddoctest
      \endlang

      \lang{matlab}
      \doctest
      >> ~DM([1; 0; -1])
      [0, 1, 0]
      >> not(DM.eye(2))
      [[0, 1],
       [1, 0]]
      >> ~SX.sym('x')
      (!x)
      \enddoctest
      \endlang

      \lang{julia}
      \doctest
      julia> ca.logic_not(ca.DM([1, 0, -1]))
      DM([0, 1, 0])
      julia> ca.logic_not(ca.DM.eye(2))
      DM(
      [[0, 1],
       [1, 0]])
      julia> ca.logic_not(ca.SX.sym("x"))
      SX((!x))
      \enddoctest
      \endlang

      \sa logic_and, logic_or

       \identifier{oz} */
    static ExType logic_not(const ExType& x) {
      return ExType::unary(OP_NOT, x);
    }
    friend inline ExType logic_not(const ExType& x) {
      return ExType::logic_not(x);
    }
    inline ExType operator!() const {
      return logic_not(self());
    }
    ///@}

    ///@{
    /** \brief Absolute value: x -> abs(x)

        Elementwise; structural zeros are preserved. The unary functions of this group operate
        on every entry, keep the shape and build a symbolic expression for SX and MX. Those
        that map 0 to 0 (abs, sqrt, sin, tan, atan, asin, sinh, tanh, asinh, atanh, floor,
        ceil, sign, erf, erfinv, log1p, expm1) keep the sparsity pattern; the others (cos,
        cosh, exp, log, log10, acos, acosh, logic_not) return a dense result. Named fabs in
        Python, where the built-in abs also applies.

        \lang{python}
        \doctest
        >>> fabs(DM([-1, 0, 2]))
        DM([1, 0, 2])
        >>> fabs(DM.eye(2) * -1)
        DM(
        [[1, 00],
         [00, 1]])
        >>> abs(SX.sym("x"))
        SX(fabs(x))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> abs(DM([-1; 0; 2]))
        [1, 0, 2]
        >> abs(-DM.eye(2))
        [[1, 00],
         [00, 1]]
        >> abs(SX.sym('x'))
        fabs(x)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> abs.(ca.DM([-1, 0, 2]))
        DM([1, 0, 2])
        julia> abs(-ca.DM.eye(2))
        DM(
        [[1, 00],
         [00, 1]])
        julia> abs(ca.SX.sym("x"))
        SX(fabs(x))
        \enddoctest
        \endlang

        \sa sign, copysign

        \identifier{p0} */
    static ExType abs(const ExType& x) {
      return ExType::unary(OP_FABS, x);
    }
    friend inline ExType abs(const ExType& x) {
      return ExType::abs(x);
    }
    friend inline ExType fabs(const ExType& x) {
      return abs(x);
    }
    ///@}

    ///@{
    /** \brief Square root: x -> sqrt(x)

        Elementwise; sparsity preserved. A negative argument gives NaN, no error is raised.
        A symbolic power x**0.5 simplifies to sqrt(x).

        \lang{python}
        \doctest
        >>> sqrt(DM([4, 9]))
        DM([2, 3])
        >>> SX.sym("x")**0.5
        SX(sqrt(x))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sqrt(DM([4; 9]))
        [2, 3]
        >> SX.sym('x').^0.5
        sqrt(x)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sqrt.(ca.DM([4, 9]))
        DM([2, 3])
        julia> ca.SX.sym("x")^0.5
        SX(sqrt(x))
        \enddoctest
        \endlang

        \sa sq, pow

        \identifier{p1} */
    static ExType sqrt(const ExType& x) {
      return ExType::unary(OP_SQRT, x);
    }
    friend inline ExType sqrt(const ExType& x) {
      return ExType::sqrt(x);
    }
    ///@}

    ///@{
    /** \brief Square: x -> x^2

        Elementwise; sparsity preserved. Not exposed as a function in Python; a symbolic
        power with exponent 2 is turned into this operation.

        \lang{python}
        \doctest
        >>> SX.sym("x")**2
        SX(sq(x))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> SX.sym('x').^2
        sq(x)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.SX.sym("x")^2
        SX(sq(x))
        \enddoctest
        \endlang

        \sa sqrt, pow

        \identifier{p2} */
    static ExType sq(const ExType& x) {
      return ExType::unary(OP_SQ, x);
    }
    friend inline ExType sq(const ExType& x) {
      return ExType::sq(x);
    }
    ///@}

    ///@{
    /** \brief Sine: x -> sin(x)

        Elementwise; sparsity preserved (sin(0) == 0).

        \lang{python}
        \doctest
        >>> sin(DM.eye(2))
        DM(
        [[0.841471, 00],
         [00, 0.841471]])
        >>> sin(MX.sym("x", 2))
        MX(sin(x))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sin(DM.eye(2))
        [[0.841471, 00],
         [00, 0.841471]]
        >> sin(MX.sym('x', 2))
        sin(x)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sin.(ca.DM.eye(2))
        DM(
        [[0.841471, 00],
         [00, 0.841471]])
        julia> sin.(ca.MX.sym("x", 2))
        MX(sin(x))
        \enddoctest
        \endlang

        \sa cos, tan, asin

        \identifier{p3} */
    static ExType sin(const ExType& x) {
      return ExType::unary(OP_SIN, x);
    }
    friend inline ExType sin(const ExType& x) {
      return ExType::sin(x);
    }
    ///@}

    ///@{
    /** \brief Cosine: x -> cos(x)

        Elementwise; structural zeros become 1, so the result is dense.

        \lang{python}
        \doctest
        >>> cos(DM.eye(2))
        DM(
        [[0.540302, 1],
         [1, 0.540302]])
        >>> cos(SX.sym("s", Sparsity.diag(2)))
        SX(@1=1,
        [[cos(s_0), @1],
         [@1, cos(s_1)]])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> cos(DM.eye(2))
        [[0.540302, 1],
         [1, 0.540302]]
        >> cos(SX.sym('s', Sparsity.diag(2)))
        @1=1,
        [[cos(s_0), @1],
         [@1, cos(s_1)]]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> cos.(ca.DM.eye(2))
        DM(
        [[0.540302, 1],
         [1, 0.540302]])
        julia> cos.(ca.SX.sym("s", ca.Sparsity.diag(2)))
        SX(@1=1,
        [[cos(s_0), @1],
         [@1, cos(s_1)]])
        \enddoctest
        \endlang

        \sa sin, tan, acos

        \identifier{p4} */
    static ExType cos(const ExType& x) {
      return ExType::unary(OP_COS, x);
    }
    friend inline ExType cos(const ExType& x) {
      return ExType::cos(x);
    }
    ///@}

    ///@{
    /** \brief Tangent: x -> tan(x)

        Elementwise; sparsity preserved.

        \lang{python}
        \doctest
        >>> tan(DM([0, 1]))
        DM([0, 1.55741])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> tan(DM([0; 1]))
        [0, 1.55741]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> tan.(ca.DM([0, 1]))
        DM([0, 1.55741])
        \enddoctest
        \endlang

        \sa atan, sin, cos

        \identifier{p5} */
    static ExType tan(const ExType& x) {
      return ExType::unary(OP_TAN, x);
    }
    friend inline ExType tan(const ExType& x) {
      return ExType::tan(x);
    }
    ///@}

    ///@{
    /** \brief Arc tangent: x -> atan(x)

        Elementwise; sparsity preserved. Result in (-pi/2, pi/2); use atan2 for a
        quadrant-aware angle.

        \lang{python}
        \doctest
        >>> atan(DM([0, 1]))
        DM([0, 0.785398])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> atan(DM([0; 1]))
        [0, 0.785398]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> atan.(ca.DM([0, 1]))
        DM([0, 0.785398])
        \enddoctest
        \endlang

        \sa atan2, tan

        \identifier{p6} */
    static ExType atan(const ExType& x) {
      return ExType::unary(OP_ATAN, x);
    }
    friend inline ExType atan(const ExType& x) {
      return ExType::atan(x);
    }
    ///@}

    ///@{
    /** \brief Arc sine: x -> asin(x)

        Elementwise; sparsity preserved. Arguments outside [-1, 1] give NaN.

        \lang{python}
        \doctest
        >>> asin(DM([0, 1]))
        DM([0, 1.5708])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> asin(DM([0; 1]))
        [0, 1.5708]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> asin.(ca.DM([0, 1]))
        DM([0, 1.5708])
        \enddoctest
        \endlang

        \sa sin, acos

        \identifier{p7} */
    static ExType asin(const ExType& x) {
      return ExType::unary(OP_ASIN, x);
    }
    friend inline ExType asin(const ExType& x) {
      return ExType::asin(x);
    }
    ///@}

    ///@{
    /** \brief Arc cosine: x -> acos(x)

        Elementwise; the result is dense (acos(0) == pi/2). Arguments outside [-1, 1] give NaN.

        \lang{python}
        \doctest
        >>> acos(DM([0, 1]))
        DM([1.5708, 0])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> acos(DM([0; 1]))
        [1.5708, 0]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> acos.(ca.DM([0, 1]))
        DM([1.5708, 0])
        \enddoctest
        \endlang

        \sa cos, asin

        \identifier{p8} */
    static ExType acos(const ExType& x) {
      return ExType::unary(OP_ACOS, x);
    }
    friend inline ExType acos(const ExType& x) {
      return ExType::acos(x);
    }
    ///@}

    ///@{
    /** \brief Hyperbolic tangent: x -> tanh(x)

        Elementwise; sparsity preserved.

        \lang{python}
        \doctest
        >>> tanh(DM([0, 1]))
        DM([0, 0.761594])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> tanh(DM([0; 1]))
        [0, 0.761594]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> tanh.(ca.DM([0, 1]))
        DM([0, 0.761594])
        \enddoctest
        \endlang

        \sa atanh, sinh, cosh

        \identifier{p9} */
    static ExType tanh(const ExType& x) {
      return ExType::unary(OP_TANH, x);
    }
    friend inline ExType tanh(const ExType& x) {
      return ExType::tanh(x);
    }
    ///@}

    ///@{
    /** \brief Hyperbolic sin: x -> sinh(x)

        Elementwise; sparsity preserved.

        \lang{python}
        \doctest
        >>> sinh(DM([0, 1]))
        DM([0, 1.1752])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sinh(DM([0; 1]))
        [0, 1.1752]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sinh.(ca.DM([0, 1]))
        DM([0, 1.1752])
        \enddoctest
        \endlang

        \sa asinh, cosh, tanh

        \identifier{pa} */
    static ExType sinh(const ExType& x) {
      return ExType::unary(OP_SINH, x);
    }
    friend inline ExType sinh(const ExType& x) {
      return ExType::sinh(x);
    }
    ///@}

    ///@{
    /** \brief Hyperbolic cosine: x -> cosh(x)

        Elementwise; the result is dense (cosh(0) == 1).

        \lang{python}
        \doctest
        >>> cosh(DM([0, 1]))
        DM([1, 1.54308])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> cosh(DM([0; 1]))
        [1, 1.54308]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> cosh.(ca.DM([0, 1]))
        DM([1, 1.54308])
        \enddoctest
        \endlang

        \sa acosh, sinh, tanh

        \identifier{pb} */
    static ExType cosh(const ExType& x) {
      return ExType::unary(OP_COSH, x);
    }
    friend inline ExType cosh(const ExType& x) {
      return ExType::cosh(x);
    }
    ///@}

    ///@{
    /** \brief Inverse hyperbolic tangent: x -> atanh(x)

        Elementwise; sparsity preserved. Defined for |x| < 1, infinite at +-1.

        \lang{python}
        \doctest
        >>> atanh(DM([0, 0.5]))
        DM([0, 0.549306])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> atanh(DM([0; 0.5]))
        [0, 0.549306]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> atanh.(ca.DM([0, 0.5]))
        DM([0, 0.549306])
        \enddoctest
        \endlang

        \sa tanh

        \identifier{pc} */
    static ExType atanh(const ExType& x) {
      return ExType::unary(OP_ATANH, x);
    }
    friend inline ExType atanh(const ExType& x) {
      return ExType::atanh(x);
    }
    ///@}

    ///@{
    /** \brief Inverse hyperbolic sin: x -> asinh(x)

        Elementwise; sparsity preserved.

        \lang{python}
        \doctest
        >>> asinh(DM([-1, 1]))
        DM([-0.881374, 0.881374])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> asinh(DM([-1; 1]))
        [-0.881374, 0.881374]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> asinh.(ca.DM([-1, 1]))
        DM([-0.881374, 0.881374])
        \enddoctest
        \endlang

        \sa sinh

        \identifier{pd} */
    static ExType asinh(const ExType& x) {
      return ExType::unary(OP_ASINH, x);
    }
    friend inline ExType asinh(const ExType& x) {
      return ExType::asinh(x);
    }
    ///@}

    ///@{
    /** \brief Inverse hyperbolic cosine: x -> acosh(x)

        Elementwise; the result is dense. Defined for x >= 1, NaN below.

        \lang{python}
        \doctest
        >>> acosh(DM([1, 2]))
        DM([0, 1.31696])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> acosh(DM([1; 2]))
        [0, 1.31696]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> acosh.(ca.DM([1, 2]))
        DM([0, 1.31696])
        \enddoctest
        \endlang

        \sa cosh

        \identifier{pe} */
    static ExType acosh(const ExType& x) {
      return ExType::unary(OP_ACOSH, x);
    }
    friend inline ExType acosh(const ExType& x) {
      return ExType::acosh(x);
    }
    ///@}

    ///@{
    /** \brief Elementwise exponential: x -> exp(x)

        Elementwise; the result is dense (exp(0) == 1).

        \lang{python}
        \doctest
        >>> exp(DM.eye(2))
        DM(
        [[2.71828, 1],
         [1, 2.71828]])
        >>> exp(MX.sym("m", Sparsity.diag(2)))
        MX(exp(dense(m)))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> exp(DM.eye(2))
        [[2.71828, 1],
         [1, 2.71828]]
        >> exp(MX.sym('m', Sparsity.diag(2)))
        exp(dense(m))
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> exp.(ca.DM.eye(2))
        DM(
        [[2.71828, 1],
         [1, 2.71828]])
        julia> exp.(ca.MX.sym("m", ca.Sparsity.diag(2)))
        MX(exp(dense(m)))
        \enddoctest
        \endlang

        \sa log, expm1

        \identifier{pf} */
    static ExType exp(const ExType& x) {
      return ExType::unary(OP_EXP, x);
    }
    friend inline ExType exp(const ExType& x) {
      return ExType::exp(x);
    }
    ///@}

    ///@{
    /** \brief Natural logarithm: x -> log(x)

        Elementwise; the result is dense with -inf at structural zeros. Negative arguments
        give NaN; no error is raised.

        \lang{python}
        \doctest
        >>> log(DM([1, 0]))
        DM([0, -inf])
        >>> log(SX.sym("x"))
        SX(log(x))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> log(DM([1; 0]))
        [0, -inf]
        >> log(SX.sym('x'))
        log(x)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> log.(ca.DM([1, 0]))
        DM([0, -inf])
        julia> log(ca.SX.sym("x"))
        SX(log(x))
        \enddoctest
        \endlang

        \sa exp, log10, log1p

        \identifier{pg} */
    static ExType log(const ExType& x) {
      return ExType::unary(OP_LOG, x);
    }
    friend inline ExType log(const ExType& x) {
      return ExType::log(x);
    }
    ///@}

    ///@{
    /** \brief Base-10 logarithm: x -> log10(x)

        Computed as log(x) * (1 / log(10)); there is no dedicated operation, which is visible
        in symbolic expressions. Dense result like log.

        \lang{python}
        \doctest
        >>> log10(DM([1, 10, 1000]))
        DM([0, 1, 3])
        >>> log10(SX.sym("x"))
        SX((0.434294*log(x)))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> log10(DM([1; 10; 1000]))
        [0, 1, 3]
        >> log10(SX.sym('x'))
        (0.434294*log(x))
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> log10.(ca.DM([1, 10, 1000]))
        DM([0, 1, 3])
        julia> log10(ca.SX.sym("x"))
        SX((0.434294*log(x)))
        \enddoctest
        \endlang

        \sa log

        \identifier{ph} */
    static ExType log10(const ExType& x) {
      return log(x)*(1/std::log(10.));
    }
    friend inline ExType log10(const ExType& x) {
      return ExType::log10(x);
    }
    ///@}

    ///@{
    /** \brief Precision variant for natural logarithm: x -> log(x+1)

        Elementwise; accurate for small x where log(1 + x) would lose digits. Sparsity
        preserved (log1p(0) == 0).

        \lang{python}
        \doctest
        >>> log1p(DM([0, 1e-10]))
        DM([0, 1e-10])
        >>> log1p(SX.sym("x"))
        SX(log1p(x))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> log1p(DM([0; 1e-10]))
        [0, 1e-10]
        >> log1p(SX.sym('x'))
        log1p(x)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> log1p.(ca.DM([0, 1e-10]))
        DM([0, 1e-10])
        julia> log1p(ca.SX.sym("x"))
        SX(log1p(x))
        \enddoctest
        \endlang

        \sa log, expm1

        \identifier{pi} */
    static ExType log1p(const ExType& x) {
      return ExType::unary(OP_LOG1P, x);
    }
    friend inline ExType log1p(const ExType& x) {
      return ExType::log1p(x);
    }
    ///@}

    ///@{
    /** \brief Precision variant for elementwise exponential: x -> exp(x)-1

        Elementwise; accurate for small x where exp(x) - 1 would lose digits. Sparsity
        preserved (expm1(0) == 0).

        \lang{python}
        \doctest
        >>> expm1(DM([0, 1e-10]))
        DM([0, 1e-10])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> expm1(DM([0; 1e-10]))
        [0, 1e-10]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> expm1.(ca.DM([0, 1e-10]))
        DM([0, 1e-10])
        \enddoctest
        \endlang

        \sa exp, log1p

        \identifier{pj} */
    static ExType expm1(const ExType& x) {
      return ExType::unary(OP_EXPM1, x);
    }
    friend inline ExType expm1(const ExType& x) {
      return ExType::expm1(x);
    }
    ///@}

    ///@{
    /** \brief Round down to nearest integer: x -> floor(x)

        Elementwise, towards minus infinity; sparsity preserved. The derivative is zero
        everywhere. Floor division x // y in Python is floor(x / y).

        \lang{python}
        \doctest
        >>> floor(DM([-1.5, 1.5]))
        DM([-2, 1])
        >>> DM([7, -7]) // 2
        DM([3, -4])
        >>> SX.sym("x") // 2
        SX(floor((x/2)))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> floor(DM([-1.5; 1.5]))
        [-2, 1]
        >> floor(DM([7; -7]) ./ 2)
        [3, -4]
        >> floor(SX.sym('x') / 2)
        floor((x/2))
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> floor.(ca.DM([-1.5, 1.5]))
        DM([-2, 1])
        julia> fld.(ca.DM([7, -7]), 2)
        DM([3, -4])
        julia> floor(ca.SX.sym("x") / 2)
        SX(floor((x/2)))
        \enddoctest
        \endlang

        \sa ceil, sign

        \identifier{pk} */
    static ExType floor(const ExType& x) {
      return ExType::unary(OP_FLOOR, x);
    }
    friend inline ExType floor(const ExType& x) {
      return ExType::floor(x);
    }
    ///@}

    ///@{
    /** \brief Round up to nearest integer: x -> ceil(x)

        Elementwise, towards plus infinity; sparsity preserved. The derivative is zero
        everywhere.

        \lang{python}
        \doctest
        >>> ceil(DM([-1.5, 1.5]))
        DM([-1, 2])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> ceil(DM([-1.5; 1.5]))
        [-1, 2]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ceil.(ca.DM([-1.5, 1.5]))
        DM([-1, 2])
        \enddoctest
        \endlang

        \sa floor

        \identifier{pl} */
    static ExType ceil(const ExType& x) {
      return ExType::unary(OP_CEIL, x);
    }
    friend inline ExType ceil(const ExType& x) {
      return ExType::ceil(x);
    }
    ///@}

    ///@{
    /** \brief Error function: x -> erf(x)

        Elementwise; sparsity preserved (erf(0) == 0). Range (-1, 1).

        \lang{python}
        \doctest
        >>> erf(DM([0, 1]))
        DM([0, 0.842701])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> erf(DM([0; 1]))
        [0, 0.842701]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.erf(ca.DM([0, 1]))
        DM([0, 0.842701])
        \enddoctest
        \endlang

        \sa erfinv

        \identifier{pm} */
    static ExType erf(const ExType& x) {
      return ExType::unary(OP_ERF, x);
    }
    friend inline ExType erf(const ExType& x) {
      return ExType::erf(x);
    }
    ///@}

    ///@{
    /** \brief Inverse error function: x -> erfinv(x)

        Elementwise; sparsity preserved. Defined on [-1, 1] with erfinv(+-1) infinite;
        arguments outside give NaN.

        \lang{python}
        \doctest
        >>> erfinv(DM([0, 0.5, 1]))
        DM([0, 0.476936, inf])
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> erfinv(DM([0; 0.5; 1]))
        [0, 0.476936, inf]
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.erfinv(ca.DM([0, 0.5, 1]))
        DM([0, 0.476936, inf])
        \enddoctest
        \endlang

        \sa erf

        \identifier{pn} */
    static ExType erfinv(const ExType& x) {
      return ExType::unary(OP_ERFINV, x);
    }
    friend inline ExType erfinv(const ExType& x) {
      return ExType::erfinv(x);
    }
    ///@}

    ///@{
    /** \brief Sign function:

        sign(x)   := -1 for x<0
        sign(x)   :=  1 for x>0,
        sign(0)   :=  0
        sign(NaN) :=  NaN

        Elementwise; sparsity preserved. The derivative is zero everywhere.

        \lang{python}
        \doctest
        >>> sign(DM([-2, 0, 3]))
        DM([-1, 0, 1])
        >>> sign(SX.sym("x"))
        SX(sign(x))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> sign(DM([-2; 0; 3]))
        [-1, 0, 1]
        >> sign(SX.sym('x'))
        sign(x)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> sign.(ca.DM([-2, 0, 3]))
        DM([-1, 0, 1])
        julia> sign(ca.SX.sym("x"))
        SX(sign(x))
        \enddoctest
        \endlang

        \sa abs, copysign

        \identifier{po} */
    static ExType sign(const ExType& x) {
      return ExType::unary(OP_SIGN, x);
    }
    friend inline ExType sign(const ExType& x) {
      return ExType::sign(x);
    }
    ///@}

    ///@{
    /** \brief Elementwise power: (x,y) -> x.^y

        Elementwise; shape rules as for plus. Named power in Python, where the ** operator
        maps to it as well. A constant exponent is simplified symbolically: integer powers
        become repeated squaring and multiplication, 0.5 becomes sqrt, -1 a reciprocal, other
        constants a constpow operation. Structural zeros of \p x stay zero for a nonzero
        exponent; x**0 is dense (0**0 == 1).

        \lang{python}
        \doctest
        >>> power(DM([1, 2, 3]), 2)
        DM([1, 4, 9])
        >>> DM([1, 2, 3])**DM([3, 2, 1]), 2**DM([1, 2, 3])
        (DM([1, 4, 3]), DM([2, 4, 8]))
        >>> x = SX.sym("x")
        >>> x**2, x**3, x**-1, x**x
        (SX(sq(x)), SX((x*sq(x))), SX((1./x)), SX(pow(x,x)))
        >>> x**2.5
        SX(pow(x,2.5))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> power(DM([1; 2; 3]), 2)
        [1, 4, 9]
        >> DM([1; 2; 3]).^DM([3; 2; 1]), 2.^DM([1; 2; 3])
        [1, 4, 3]
        [2, 4, 8]
        >> x = SX.sym('x');
        >> x.^2, x.^3, x.^-1, x.^x
        sq(x)
        (x*sq(x))
        (1./x)
        pow(x,x)
        >> x.^2.5
        pow(x,2.5)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> ca.power(ca.DM([1, 2, 3]), 2)
        DM([1, 4, 9])
        julia> ca.DM([1, 2, 3]) .^ ca.DM([3, 2, 1]), 2 .^ ca.DM([1, 2, 3])
        ([1, 4, 3], [2, 4, 8])
        julia> x = ca.SX.sym("x");
        julia> x^2, x^3, x^-1, x^x
        (sq(x), (x*sq(x)), (1./x), pow(x,x))
        julia> x^2.5
        SX(pow(x,2.5))
        \enddoctest
        \endlang

        \sa sq, sqrt, constpow, mpower

        \identifier{pp} */
    static ExType pow(const ExType& x, const ExType& y) {
      return ExType::binary(OP_POW, x, y);
    }
    friend inline ExType pow(const ExType& x, const ExType& y) {
      return ExType::pow(x, y);
    }
    ///@}

    ///@{
    /** \brief Remainder after division: (x,y) -> fmod(x,y)

    This Function follows the convention of https://en.cppreference.com/w/c/numeric/math/fmod:
    the result has the sign of \p x and magnitude below |y|.

    This is equivalent to Python's numpy.fmod and Matlab's rem. It differs from the Python
    operator x % y on these objects, which follows the Python convention x - y*floor(x/y)
    (result with the sign of \p y). Elementwise, shape rules as for plus; structural zeros of
    \p x stay zero. A zero divisor gives NaN.

    \lang{python}
    \doctest
    >>> fmod(DM([5, 5, -5, -5]), DM([3, -3, 3, -3]))
    DM([2, 2, -2, -2])
    >>> DM([5, 5, -5, -5]) % DM([3, -3, 3, -3])
    DM([2, -1, 1, -2])
    >>> fmod(SX.sym("x"), 2), SX.sym("x") % 2
    (SX(fmod(x,2)), SX((x-(2.*floor((x/2))))))
    \enddoctest
    \endlang

    \lang{matlab}
    \doctest
    >> x = DM([5; 5; -5; -5]); y = DM([3; -3; 3; -3]);
    >> fmod(x, y)
    [2, 2, -2, -2]
    >> rem(x, y)
    [2, 2, -2, -2]
    >> x - y .* floor(x ./ y)
    [2, -1, 1, -2]
    >> rem(SX.sym('x'), 2)
    fmod(x,2)
    \enddoctest
    \endlang

    \lang{julia}
    \doctest
    julia> ca.fmod(ca.DM([5, 5, -5, -5]), ca.DM([3, -3, 3, -3]))
    DM([2, 2, -2, -2])
    julia> mod.(ca.DM([5, 5, -5, -5]), ca.DM([3, -3, 3, -3]))
    DM([2, -1, 1, -2])
    julia> ca.fmod(ca.SX.sym("x"), 2)
    SX(fmod(x,2))
    \enddoctest
    \endlang

    \sa remainder, floor

        \identifier{pq} */
    static ExType mod(const ExType& x, const ExType& y) {
      return ExType::binary(OP_FMOD, x, y);
    }
    friend inline ExType mod(const ExType& x, const ExType& y) {
      return ExType::mod(x, y);
    }
    friend inline ExType fmod(const ExType& x, const ExType& y) {
      return mod(x, y);
    }
    ///@}

    ///@{
    /** \brief Remainder after division: (x,y) -> remainder(x,y)

    This Function follows the convention of
    https://en.cppreference.com/w/c/numeric/math/remainder: x - n*y with n the integer nearest
    to x/y (ties to even), so the result lies in [-|y|/2, |y|/2].

    Notably:
      - remainder(5,3)   -> -1
      - remainder(5,-3)  -> -1
      - remainder(-5,3)  -> 1
      - remainder(-5,-3) -> 1

    This is equivalent to Python's math.remainder. There is no equivalence in Matlab.
    Elementwise, shape rules as for plus; structural zeros of \p x stay zero.

    \lang{python}
    \doctest
    >>> remainder(DM([5, 5, -5, -5]), DM([3, -3, 3, -3]))
    DM([-1, -1, 1, 1])
    >>> remainder(DM([2.5, 3.5]), 2)
    DM([0.5, -0.5])
    >>> remainder(MX.sym("x"), 2)
    MX(remainder(x,2))
    \enddoctest
    \endlang

    \lang{matlab}
    \doctest
    >> remainder(DM([5; 5; -5; -5]), DM([3; -3; 3; -3]))
    [-1, -1, 1, 1]
    >> remainder(DM([2.5; 3.5]), 2)
    [0.5, -0.5]
    >> remainder(MX.sym('x'), 2)
    remainder(x,2)
    \enddoctest
    \endlang

    \lang{julia}
    \doctest
    julia> ca.remainder(ca.DM([5, 5, -5, -5]), ca.DM([3, -3, 3, -3]))
    DM([-1, -1, 1, 1])
    julia> ca.remainder(ca.DM([2.5, 3.5]), 2)
    DM([0.5, -0.5])
    julia> ca.remainder(ca.MX.sym("x"), 2)
    MX(remainder(x,2))
    \enddoctest
    \endlang

    \sa fmod

        \identifier{24x} */
    static ExType remainder(const ExType& x, const ExType& y) {
      return ExType::binary(OP_REMAINDER, x, y);
    }
    friend inline ExType remainder(const ExType& x, const ExType& y) {
      return ExType::remainder(x, y);
    }
    ///@}

    ///@{
    /** \brief Two argument arc tangent: (y,x) -> atan2(y,x)
     *
     * theta = atan2(y,x) corresponds to x = r cos(theta), y = r sin(theta)

        Angle in (-pi, pi] of the point (x, y), taking the quadrant into account; the first
        argument is the y coordinate, as in C and numpy. Elementwise, shape rules as for plus.
        atan2(0, 0) is 0.

        \lang{python}
        \doctest
        >>> atan2(DM([1, -1, 0]), DM([-1, -1, 1]))
        DM([2.35619, -2.35619, 0])
        >>> atan2(SX.sym("y"), SX.sym("x"))
        SX(atan2(y,x))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> atan2(DM([1; -1; 0]), DM([-1; -1; 1]))
        [2.35619, -2.35619, 0]
        >> atan2(SX.sym('y'), SX.sym('x'))
        atan2(y,x)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> atan.(ca.DM([1, -1, 0]), ca.DM([-1, -1, 1]))
        DM([2.35619, -2.35619, 0])
        julia> atan(ca.SX.sym("y"), ca.SX.sym("x"))
        SX(atan2(y,x))
        \enddoctest
        \endlang

        \sa atan

        \identifier{pr} */
    static ExType atan2(const ExType& y, const ExType& x) {
      return ExType::binary(OP_ATAN2, y, x);
    }
    friend inline ExType atan2(const ExType& y, const ExType& x) {
      return ExType::atan2(y, x);
    }
    ///@}

    ///@{
    /** \brief Conditional assignment: (x,y) -> x ? y : 0

        Elementwise: \p y where \p x is nonzero, 0 elsewhere; the result is structurally zero
        wherever either operand is. The value of \p y in the inactive branch does not
        propagate, so a NaN or inf there yields 0 rather than NaN, and the derivatives are
        treated accordingly. This is the building block of if_else; it is not exposed as a
        function in Python, use if_else(x, y, 0) there.

        \sa if_else

        \identifier{ps} */
    static ExType if_else_zero(const ExType& x, const ExType& y) {
      return ExType::binary(OP_IF_ELSE_ZERO, x, y);
    }
    friend inline ExType if_else_zero(const ExType& x, const ExType& y) {
      return ExType::if_else_zero(x, y);
    }
    ///@}

    ///@{
    /** \brief Smallest of two values: (x,y) -> min(x,y)

        Elementwise minimum; shape rules as for plus. Follows C fmin: when one operand is NaN
        the other is returned, NaN results only when both are NaN. The derivative is 1 with
        respect to the smaller operand (split evenly on ties). The sparsity pattern is kept
        where both operands are structurally zero; fmin with a negative scalar densifies.
        For the minimum over the entries of one matrix use mmin.

        \lang{python}
        \doctest
        >>> fmin(DM([1, 5]), DM([3, 2]))
        DM([1, 2])
        >>> fmin(DM([1, float("nan")]), DM([float("nan"), 2]))
        DM([1, 2])
        >>> fmin(DM.eye(2), 1)
        DM(
        [[1, 00],
         [00, 1]])
        >>> fmin(SX.sym("x"), SX.sym("y"))
        SX(fmin(x,y))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> fmin(DM([1; 5]), DM([3; 2]))
        [1, 2]
        >> min(DM([1; 5]), DM([3; 2]))
        [1, 2]
        >> fmin(DM([1; nan]), DM([nan; 2]))
        [1, 2]
        >> fmin(DM.eye(2), 1)
        [[1, 00],
         [00, 1]]
        >> min(SX.sym('x'), SX.sym('y'))
        fmin(x,y)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> min.(ca.DM([1, 5]), ca.DM([3, 2]))
        DM([1, 2])
        julia> min.(ca.DM([1, NaN]), ca.DM([NaN, 2]))
        DM([1, 2])
        julia> min.(ca.DM.eye(2), 1)
        DM(
        [[1, 00],
         [00, 1]])
        julia> min(ca.SX.sym("x"), ca.SX.sym("y"))
        SX(fmin(x,y))
        \enddoctest
        \endlang

        \sa fmax, mmin

        \identifier{pt} */
    static ExType fmin(const ExType& x, const ExType& y) {
      return ExType::binary(OP_FMIN, x, y);
    }
    friend inline ExType fmin(const ExType& x, const ExType& y) {
      return ExType::fmin(x, y);
    }
    ///@}

    ///@{
    /** \brief Largest of two values: (x,y) -> max(x,y)

        Elementwise maximum; shape rules as for plus. Follows C fmax: when one operand is NaN
        the other is returned. The derivative is 1 with respect to the larger operand (split
        evenly on ties). For the maximum over the entries of one matrix use mmax.

        \lang{python}
        \doctest
        >>> fmax(DM([1, 5]), DM([3, 2]))
        DM([3, 5])
        >>> fmax(DM([1, 5]), 3)
        DM([3, 5])
        >>> fmax(DM([1, float("nan")]), DM([float("nan"), 2]))
        DM([1, 2])
        >>> fmax(MX.sym("x"), MX.sym("y"))
        MX(fmax(x,y))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> fmax(DM([1; 5]), DM([3; 2]))
        [3, 5]
        >> max(DM([1; 5]), 3)
        [3, 5]
        >> fmax(DM([1; nan]), DM([nan; 2]))
        [1, 2]
        >> max(MX.sym('x'), MX.sym('y'))
        fmax(x,y)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> max.(ca.DM([1, 5]), ca.DM([3, 2]))
        DM([3, 5])
        julia> max.(ca.DM([1, 5]), 3)
        DM([3, 5])
        julia> max.(ca.DM([1, NaN]), ca.DM([NaN, 2]))
        DM([1, 2])
        julia> max(ca.MX.sym("x"), ca.MX.sym("y"))
        MX(fmax(x,y))
        \enddoctest
        \endlang

        \sa fmin, mmax

        \identifier{pu} */
    static ExType fmax(const ExType& x, const ExType& y) {
      return ExType::binary(OP_FMAX, x, y);
    }
    friend inline ExType fmax(const ExType& x, const ExType& y) {
      return ExType::fmax(x, y);
    }
    ///@}

    ///@{
    /** \brief Check if two nodes are equivalent up to a given depth.

     * Depth=0 checks if the expressions are identical, i.e. points to the same node.
     *
     * a = x*x
     * b = x*x
     *
     *  is_equal(a,b,0)  will return false, but is_equal(a,b,1) will return true

        Returns a plain boolean, unlike the elementwise comparison eq. With \p depth > 0 two
        different nodes are equivalent when they have the same operation and their
        dependencies are equivalent up to depth-1. For SX and DM matrices the test is applied
        entry by entry; matrices of different shape raise an error, and for DM it is value
        equality. Both arguments must be of the same type.

        \param x First expression
        \param y Second expression
        \param depth Number of levels of the expression graph to compare, 0 for identity

        \lang{python}
        \doctest
        >>> x = SX.sym("x")
        >>> is_equal(x, x), is_equal(x, SX.sym("x"))
        (True, False)
        >>> is_equal(x*x, x*x), is_equal(x*x, x*x, 1)
        (False, True)
        >>> is_equal(DM([1, 2]), DM([1, 2])), is_equal(DM([1, 2]), DM([1, 3]))
        (True, False)
        >>> y = MX.sym("y")
        >>> is_equal(y * y, y * y), is_equal(y * y, y * y, 1)
        (False, True)
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> x = SX.sym('x');
        >> is_equal(x, x), is_equal(x, SX.sym('x'))
        1
        0
        >> is_equal(x*x, x*x), is_equal(x*x, x*x, 1)
        0
        1
        >> is_equal(DM([1; 2]), DM([1; 2])), is_equal(DM([1; 2]), DM([1; 3]))
        1
        0
        >> y = MX.sym('y');
        >> is_equal(y * y, y * y), is_equal(y * y, y * y, 1)
        0
        1
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> x = ca.SX.sym("x");
        julia> ca.is_equal(x, x), ca.is_equal(x, ca.SX.sym("x"))
        (true, false)
        julia> ca.is_equal(x*x, x*x), ca.is_equal(x*x, x*x, 1)
        (false, true)
        julia> ca.is_equal(ca.DM([1, 2]), ca.DM([1, 2])), ca.is_equal(ca.DM([1, 2]), ca.DM([1, 3]))
        (true, false)
        julia> y = ca.MX.sym("y");
        julia> ca.is_equal(y * y, y * y), ca.is_equal(y * y, y * y, 1)
        (false, true)
        \enddoctest
        \endlang

        \sa eq

        \identifier{pv} */
     friend inline bool is_equal(const ExType& x, const ExType& y, casadi_int depth=0) {
       return ExType::is_equal(x, y, depth);
     }
     ///@}

     ///@{
     /** \brief Copy sign: (x,y) -> copysign(x,y)

         Elementwise value with the magnitude of \p x and the sign of \p y, as C copysign;
         a negative zero \p y counts as negative. Shape rules as for plus. The derivative with
         respect to \p y is zero.

         \lang{python}
         \doctest
         >>> copysign(DM([1, -2]), DM([-1, 1]))
         DM([-1, 2])
         >>> copysign(DM(3), DM(-0.0))
         DM(-3)
         >>> copysign(SX.sym("x"), SX.sym("y"))
         SX(copysign(x,y))
         \enddoctest
         \endlang

         \lang{matlab}
         \doctest
         >> copysign(DM([1; -2]), DM([-1; 1]))
         [-1, 2]
         >> copysign(DM(3), DM(-0.0))
         -3
         >> copysign(SX.sym('x'), SX.sym('y'))
         copysign(x,y)
         \enddoctest
         \endlang

         \lang{julia}
         \doctest
         julia> copysign.(ca.DM([1, -2]), ca.DM([-1, 1]))
         DM([-1, 2])
         julia> copysign.(ca.DM(3), ca.DM(-0.0))
         DM(-3)
         julia> ca.copysign(ca.SX.sym("x"), ca.SX.sym("y"))
         SX(copysign(x,y))
         \enddoctest
         \endlang

         \sa sign, abs

         \identifier{2ol} */
     static ExType copysign(const ExType& x, const ExType& y) {
       return ExType::binary(OP_COPYSIGN, x, y);
     }
     friend inline ExType copysign(const ExType& x, const ExType& y) {
       return ExType::copysign(x, y);
     }
     ///@}

     ///@{
     /** \brief Elementwise power with const power: (x,y) -> x^y with y constant

         Same value as pow, but the exponent is assumed constant: the derivative with respect
         to \p y is taken as zero. pow itself produces this operation for a non-integer
         constant exponent, so it is rarely needed directly. Passing a symbolic exponent
         gives a wrong derivative for MX; SX falls back to pow.

         \lang{python}
         \doctest
         >>> constpow(DM(2), 3)
         DM(8)
         >>> x = MX.sym("x")
         >>> constpow(x, 2.5), jacobian(constpow(x, 2.5), x)
         (MX(pow(x,2.5)), MX((2.5*pow(x,1.5))))
         \enddoctest
         \endlang

         \lang{matlab}
         \doctest
         >> constpow(DM(2), 3)
         8
         >> x = MX.sym('x');
         >> constpow(x, 2.5), jacobian(constpow(x, 2.5), x)
         pow(x,2.5)
         (2.5*pow(x,1.5))
         \enddoctest
         \endlang

         \lang{julia}
         \doctest
         julia> ca.constpow(ca.DM(2), 3)
         DM(8)
         julia> x = ca.MX.sym("x");
         julia> ca.constpow(x, 2.5), ca.jacobian(ca.constpow(x, 2.5), x)
         (pow(x,2.5), (2.5*pow(x,1.5)))
         \enddoctest
         \endlang

         \sa pow

         \identifier{2om} */
     static ExType constpow(const ExType& x, const ExType& y) {
       return ExType::binary(OP_CONSTPOW, x, y);
     }
     friend inline ExType constpow(const ExType& x, const ExType& y) {
       return ExType::constpow(x, y);
     }
     ///@}

     ///@{
     /** \brief Debug printing: (x,y) -> x, printing x when evaluated

         Identity operation with a side effect: whenever the expression is evaluated inside a
         Function (and hence inside a solver), every nonzero of \p x is printed on its own
         line as "|> y : value" in storage (column-major) order, with \p y an identifier of
         choice to tell several printme nodes apart. Available as a method x.printme(y) in
         Python; the derivative is that of \p x. Use monitor to print a whole matrix at once.

         \lang{python}
         \doctest
         >>> x = MX.sym("x")
         >>> f = Function("f", [x], [x.printme(7) * 2])
         >>> f(3)
         |> 7 : 3.0000000000000000e+00
         DM(6)
         \enddoctest
         \endlang

         \lang{matlab}
         \doctest
         >> x = MX.sym('x');
         >> f = Function('f', {x}, {x.printme(7) * 2});
         >> f(3)
         |> 7 : 3.0000000000000000e+00
         6
         \enddoctest
         \endlang

         \lang{julia}
         \doctest
         julia> x = ca.MX.sym("x");
         julia> f = ca.Function("f", [x], [ca.printme(x, 7) * 2]);
         julia> f(3)
         |> 7 : 3.0000000000000000e+00
         DM(6)
         \enddoctest
         \endlang

         \identifier{2on} */
     static ExType printme(const ExType& x, const ExType& y) {
       return ExType::binary(OP_PRINTME, x, y);
     }
     friend inline ExType printme(const ExType& x, const ExType& y) {
       return ExType::binary(OP_PRINTME, x, y);
     }
     ///@}

    ///@{
    /** \brief Precision variant for 2 norm: (x,y) -> sqrt(x^2+y^2)

        Elementwise; avoids the overflow and underflow of forming the squares, as C hypot.
        Shape rules as for plus; the result is structurally zero where both operands are.

        \lang{python}
        \doctest
        >>> hypot(DM(3), DM(4))
        DM(5)
        >>> hypot(DM(1e200), DM(1e200))
        DM(1.41421e+200)
        >>> hypot(SX.sym("x"), SX.sym("y"))
        SX(hypot(x,y))
        \enddoctest
        \endlang

        \lang{matlab}
        \doctest
        >> hypot(DM(3), DM(4))
        5
        >> hypot(DM(1e200), DM(1e200))
        1.41421e+200
        >> hypot(SX.sym('x'), SX.sym('y'))
        hypot(x,y)
        \enddoctest
        \endlang

        \lang{julia}
        \doctest
        julia> hypot(ca.DM(3), ca.DM(4))
        DM(5)
        julia> hypot(ca.DM(1e200), ca.DM(1e200))
        DM(1.41421e+200)
        julia> hypot(ca.SX.sym("x"), ca.SX.sym("y"))
        SX(hypot(x,y))
        \enddoctest
        \endlang

        \sa norm_2

        \identifier{pw} */
    static ExType hypot(const ExType& x, const ExType& y) {
      return ExType::binary(OP_HYPOT, x, y);
    }
    friend inline ExType hypot(const ExType& x, const ExType& y) {
      return ExType::hypot(x, y);
    }
    ///@}


/** @} */



};
#endif // SWIG

} // namespace casadi

#endif // CASADI_GENERIC_EXPRESSION_HPP
