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


#ifndef CASADI_SX_ELEMENT_HPP
#define CASADI_SX_ELEMENT_HPP

// exception class
#include "printable.hpp"
#include "exception.hpp"
#include "casadi_limits.hpp"
#include "matrix_fwd.hpp"
#include "generic_expression.hpp"
#include "serializing_stream.hpp"

/** \brief  C/C++

    \identifier{10g} */
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace casadi {

  /** \brief  forward declaration of Node and Matrix

      \identifier{10h} */
  class SXNode; // include will follow in the end

  /** SXElem is exposed only as an empty struct to SWIG */
#ifdef SWIG
  struct SXElem {};
#else // SWIG

  /** \brief The basic scalar symbolic class of CasADi

      \author Joel Andersson
      \date 2010-2014

      \identifier{10i} */
  class CASADI_EXPORT SXElem : public GenericExpression<SXElem>,
                               public Printable<SXElem> {
    friend class SXNode;
    friend class BinarySXNode;
    friend class Matrix<SXElem>;
  public:

    /// \cond CLUTTER
    /** \brief Default constructor (not-a-number)

        Object is initialized as not-a-number.

        \identifier{10j} */
    SXElem();
    /// \endcond

    /** \brief Numerical constant constructor

        \param val Numerical value

        \identifier{10k} */
    SXElem(double val);

    /** \brief Create a symbolic primitive

         \param name Name of the symbolic primitive

        This is the name that will be used by the "operator<<" and "str" methods.
        The name is not used as identifier; you may construct distinct
        SXElem objects with non-unique names.

        \identifier{10l} */
    static SXElem sym(const std::string& name);

    /// \cond INTERNAL
    /// Create an expression from a node: extra dummy argument to avoid ambiguity for 0/NULL
    SXElem(SXNode* node, bool dummy);
    /// \endcond

    /** \brief Copy constructor

        \identifier{10m} */
    SXElem(const SXElem& scalar); // copy constructor

    /// Destructor
    ~SXElem();

    /// \cond INTERNAL
    /// Create an object given a node
    static SXElem create(SXNode* node);
    /// \endcond

    /// Assignment
    SXElem& operator=(const SXElem& scalar);
    SXElem& operator=(double scalar); // needed since otherwise both a = SXElem(double)
                                         // and a = Matrix(double) would be ok

    /// Convert to a 1-by-1 Matrix
    operator Matrix<SXElem>() const;

    /// Type name
    static std::string type_name() {return "SXElem";}

    /// Print a description of the object
    void disp(std::ostream& stream, bool more=false) const;

    /// \cond INTERNAL
    /** \brief  Get a pointer to the node

        \identifier{10n} */
    SXNode* get() const; // note: constant pointer, not pointer to constant object!
                         // (to allow access to the counter)

    /** \brief  Access functions of the node

        \identifier{10o} */
    const SXNode* operator->() const;
    SXNode* operator->();
    /// \endcond

    /** \brief  Perform operations by ID

        \identifier{10p} */
    static SXElem binary(casadi_int op, const SXElem& x, const SXElem& y);
    static SXElem unary(casadi_int op, const SXElem& x);

    /** \brief Check the truth value of this node

     * Introduced to catch bool(x) situations in python

        \identifier{10q} */
    bool __nonzero__() const;

    /** \brief check if this SXElem is a leaf of the SX graph
     *
     * An SXElem qualifies as leaf when it has no dependencies.

        \identifier{10r} */
    bool is_leaf() const;
    bool is_constant() const;
    bool is_integer() const;
    bool is_symbolic() const;
    /** \brief Check whether a binary SXElem is commutative

        \identifier{10s} */
    bool is_commutative() const;
    bool is_zero() const;
    bool is_almost_zero(double tol) const;
    bool is_one() const;
    bool is_minus_one() const;
    bool is_nan() const;
    bool is_inf() const;
    bool is_minus_inf() const;
    const std::string& name() const;
    casadi_int op() const;
    bool is_op(casadi_int op) const;

    /// Checks if expression does not contain NaN or Inf
    bool is_regular() const;

    /** \brief Check if a value is always nonnegative (false negatives are allowed)

        \identifier{10t} */
    bool is_nonnegative() const;
    SXElem dep(casadi_int ch=0) const;

    /// Type conversion to double
    explicit operator double() const;

    /// Type conversion to casadi_int
    explicit operator casadi_int() const;

    /** \brief Check if the node is the sum of two equal expressions

        \identifier{10u} */
    bool is_doubled() const;

    /** \brief Get the number of dependencies of a binary SXElem

        \identifier{10v} */
    casadi_int n_dep() const;

    /** \brief Returns a number that is unique for a given SXNode.

     * If the SXElem does not point to any node, 0 is returned.

        \identifier{10w} */
    casadi_int __hash__() const;

    /** \brief  Negation

        \identifier{10x} */
    SXElem operator-() const;

    /** \brief Element-wise inverse

        \identifier{10y} */
    SXElem inv() const;

    /** \brief Check equality up to a given depth

        \identifier{10z} */
    static bool is_equal(const SXElem& x, const SXElem& y, casadi_int depth=0);

    /// \cond INTERNAL
    /// Get the temporary variable
    int get_temp() const;

    /// Set the temporary variable
    void set_temp(int t) const;

    /// Check if marked (i.e. temporary is negative)
    bool marked() const;

    /// Mark by flipping the sign of the temporary and decreasing by one
    void mark() const;

    /** \brief Assign to another expression, if a duplicate.

     * Check for equality up to a given depth

        \identifier{110} */
    void assignIfDuplicate(const SXElem& scalar, casadi_int depth=1);

    /** \brief Assign the node to something, without invoking the deletion of the node,

     * if the count reaches 0

        \identifier{111} */
    SXNode* assignNoDelete(const SXElem& scalar);
    /// \endcond

    /** \brief SXElem nodes are not allowed to be null

        \identifier{112} */
    inline bool is_null() {return false;}

    /** \brief Ternary if_else: x ? y : z

        \identifier{113} */
    friend inline SXElem if_else(const SXElem& x, const SXElem& y, const SXElem& z) {
      return if_else_zero(x, y) + if_else_zero(!x, z);
    }

    /** \brief Serialize an object

        \identifier{114} */
    void serialize(SerializingStream& s) const;

    static SXElem deserialize(DeserializingStream& s);
  private:
    /// Pointer to node (SXElem is only a reference class)
    SXNode* node;
  };

  template<>
  class CASADI_EXPORT casadi_limits<SXElem>{
  public:
    static bool is_zero(const SXElem& val);
    static bool is_equal(const SXElem& x, const SXElem& y, casadi_int depth);
    static bool is_almost_zero(const SXElem& val, double tol);
    static bool is_one(const SXElem& val);
    static bool is_minus_one(const SXElem& val);
    static bool is_constant(const SXElem& val);
    static bool is_integer(const SXElem& val);
    static bool is_inf(const SXElem& val);
    static bool is_minus_inf(const SXElem& val);
    static bool is_nan(const SXElem& val);

    static const SXElem zero;
    static const SXElem one;
    static const SXElem two;
    static const SXElem minus_one;
    static const SXElem nan;
    static const SXElem inf;
    static const SXElem minus_inf;
  };

#endif // SWIG
/// \endcond

} // namespace casadi

#ifndef SWIG

namespace std {
  template<>
  class CASADI_EXPORT numeric_limits<casadi::SXElem>{
  public:
    static const bool is_specialized = true;
    static casadi::SXElem min() throw();
    static casadi::SXElem max() throw();
    static const int  digits = 0;
    static const int  digits10 = 0;
    static const bool is_signed = false;
    static const bool is_integer = false;
    static const bool is_exact = false;
    static const int radix = 0;
    static casadi::SXElem epsilon() throw();
    static casadi::SXElem round_error() throw();
    static const int  min_exponent = 0;
    static const int  min_exponent10 = 0;
    static const int  max_exponent = 0;
    static const int  max_exponent10 = 0;

    static const bool has_infinity = true;
    static const bool has_quiet_NaN = true;
    static const bool has_signaling_NaN = false;
    //    static const float_denorm_style has_denorm = denorm absent;
    static const bool has_denorm_loss = false;
    static casadi::SXElem infinity() throw();
    static casadi::SXElem quiet_NaN() throw();
    //    static SXElem signaling_NaN() throw();
    //    static SXElem denorm_min() throw();
    static const bool is_iec559 = false;
    static const bool is_bounded = false;
    static const bool is_modulo = false;

    static const bool traps = false;
    static const bool tinyness_before = false;
    static const float_round_style round_style = round_toward_zero;
  };
} //namespace std

#endif // SWIG

#endif // CASADI_SX_ELEMENT_HPP
