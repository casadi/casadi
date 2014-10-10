/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
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
#include "../printable_object.hpp"
#include "../casadi_exception.hpp"
#include "../casadi_limits.hpp"
#include "../matrix/matrix.hpp"
#include "../matrix/generic_expression.hpp"

/** \brief  C/C++ */
#include <iostream>
#include <string>
#include <sstream>
#include <limits>
#include <cmath>
#include <vector>

namespace casadi {

  /** \brief  forward declaration of Node and Matrix */
  class SXNode; // include will follow in the end

#ifdef SWIG
#ifdef WITH_IMPLICITCONV
  %implicitconv SXElement;
#endif // WITH_IMPLICITCONV
#endif // SWIG

  /** \brief The basic scalar symbolic class of CasADi
      \author Joel Andersson
      \date 2010-2014
  */
  class CASADI_CORE_EXPORT SXElement : public GenericExpression<SXElement>,
                                       public PrintableObject<SXElement> {
    friend class SXNode;
    friend class BinarySXNode;
    friend class Matrix<SXElement>;
  public:

    /// \cond CLUTTER
    /** \brief Default constructor (not-a-number)
        Object is initialized as not-a-number.
    */
    SXElement();
    /// \endcond

    /** \brief Numerical constant constructor
        \param val Numerical value
    */
    SXElement(double val);

    /** \brief Create a symbolic primitive
         \param name Name of the symbolic primitive

        This is the name that will be used by the "operator<<" and "toString" methods.
        The name is not used as identifier; you may construct distinct
        SXElement objects with non-unique names.
    */
    static SXElement sym(const std::string& name);

#ifndef SWIG

    /// \cond INTERNAL
    /// Create an expression from a node: extra dummy argument to avoid ambiguity for 0/NULL
    SXElement(SXNode* node, bool dummy);
    /// \endcond

    /** \brief Copy constructor */
    SXElement(const SXElement& scalar); // copy constructor

    /// Destructor
    ~SXElement();

    /// \cond INTERNAL
    /// Create an object given a node
    static SXElement create(SXNode* node);
    /// \endcond

    /// Assignment
    SXElement& operator=(const SXElement& scalar);
    SXElement& operator=(double scalar); // needed since otherwise both a = SXElement(double)
                                         // and a = Matrix(double) would be ok

    /// Convert to a 1-by-1 Matrix
    operator Matrix<SXElement>() const;

    /// Print a representation of the object
    void repr(std::ostream &stream=std::cout, bool trailing_newline=true) const;

    /// Print a description of the object
    void print(std::ostream &stream=std::cout, bool trailing_newline=true) const;

    /** \brief  print to stream, limited */
    void print(std::ostream &stream, long& remaining_calls) const;

    /// \cond INTERNAL
    /** \brief  Get a pointer to the node */
    SXNode* get() const; // note: constant pointer, not pointer to constant object!
                         // (to allow access to the counter)

    /** \brief  Access functions of the node */
    const SXNode* operator->() const;
    SXNode* operator->();
    /// \endcond

#endif // SWIG

    /** \brief  Perform operations by ID */
    static SXElement binary(int op, const SXElement& x, const SXElement& y);
    static SXElement unary(int op, const SXElement& x);

    /** \brief Check the truth value of this node
     * Introduced to catch bool(x) situations in python
     */
    bool __nonzero__() const;

    /** \brief check if this SXElement is a leaf of the SX graph
     *
     * An SXElement qualifies as leaf when it has no dependencies.
     */
    bool isLeaf() const;
    bool isConstant() const;
    bool isInteger() const;
    bool isSymbolic() const;
    bool hasDep() const;
    /** \brief Check whether a binary SXElement is commutative*/
    bool isCommutative() const;
    bool isZero() const;
    bool isAlmostZero(double tol) const;
    bool isOne() const;
    bool isMinusOne() const;
    bool isNan() const;
    bool isInf() const;
    bool isMinusInf() const;
    const std::string& getName() const;
    int getOp() const;
    bool isOp(int op) const;

    /// Checks if expression does not contain NaN or Inf
    bool isRegular() const;

    /** \brief Check if two nodes are equivalent up to a given depth.
     *  Depth=0 checks if the expressions are identical, i.e. points to the same node.
     *
     *  a = x*x
     *  b = x*x
     *
     *  a.isEqual(b, 0)  will return false, but a.isEqual(b, 1) will return true
     */
    bool isEqual(const SXElement& scalar, int depth=0) const;

    /** \brief Check if a value is always nonnegative (false negatives are allowed) */
    bool isNonNegative() const;

    double getValue() const;
    int getIntValue() const;
    SXElement getDep(int ch=0) const;

    /** \brief Check if the node is the sum of two equal expressions */
    bool isDoubled() const;

    /** \brief Get the number of dependencies of a binary SXElement */
    int getNdeps() const;

    /** \brief Returns a number that is unique for a given SXNode.
     * If the SXElement does not point to any node, 0 is returned.
     */
    long __hash__() const;

    /** \brief  Negation */
    SXElement operator-() const;

    //  all binary operations
    SXElement __add__(const SXElement& y) const;
    SXElement __sub__(const SXElement& y) const;
    SXElement __mul__(const SXElement& y) const;
    SXElement __div__(const SXElement& y) const;
    SXElement __lt__(const SXElement& y) const;
    SXElement __le__(const SXElement& y) const;
    SXElement __eq__(const SXElement& y) const;
    SXElement __ne__(const SXElement& y) const;
    using GenericExpression<SXElement>::__gt__;
    using GenericExpression<SXElement>::__ge__;
    using GenericExpression<SXElement>::__mldivide__;
    SXElement __truediv__(const SXElement &y) const {return __div__(y);}
    SXElement __pow__(const SXElement& b) const;
    SXElement __constpow__(const SXElement& b) const;

    SXElement __mrdivide__(const SXElement& b) const {  return *this / b;}
    SXElement __mpower__(const SXElement& b) const {return (*this).__pow__(b);}
    SXElement trans() const { return *this;}

    // The following functions serves two purposes:
    // Numpy compatibility and to allow unambiguous access
    SXElement mul(const SXElement& y) const { return __mul__(y);}
    SXElement exp() const;
    SXElement log() const;
    SXElement sqrt() const;
    SXElement sq() const;
    SXElement sin() const;
    SXElement cos() const;
    SXElement tan() const;
    SXElement arcsin() const;
    SXElement arccos() const;
    SXElement arctan() const;
    SXElement floor() const;
    SXElement ceil() const;
    SXElement fmod(const SXElement &y) const;
    SXElement erf() const;
    SXElement erfinv() const;
    SXElement fabs() const;
    SXElement fmin(const SXElement &y) const;
    SXElement fmax(const SXElement &y) const;
    SXElement inv() const;
    SXElement sinh() const;
    SXElement cosh() const;
    SXElement tanh() const;
    SXElement arcsinh() const;
    SXElement arccosh() const;
    SXElement arctanh() const;
    SXElement arctan2(const SXElement &y) const;
    SXElement log10() const;
    SXElement printme(const SXElement &y) const;
    SXElement sign() const;
    SXElement __copysign__(const SXElement &y) const;
    SXElement constpow(const SXElement& y) const;
    SXElement logic_not() const;
    SXElement logic_and(const SXElement& y) const;
    SXElement logic_or(const SXElement& y) const;
    SXElement if_else_zero(const SXElement& y) const;

    Matrix<SXElement> fmin(const Matrix<SXElement>& b) const;
    Matrix<SXElement> fmax(const Matrix<SXElement>& b) const;
    Matrix<SXElement> constpow(const Matrix<SXElement>& n) const;
    Matrix<SXElement> __copysign__(const Matrix<SXElement>& n) const;
    Matrix<SXElement> arctan2(const Matrix<SXElement>& b) const;

    /// \cond INTERNAL
    /// Get the temporary variable
    int getTemp() const;

    /// Set the temporary variable
    void setTemp(int t);

    /// Check if marked (i.e. temporary is negative)
    bool marked() const;

    /// Mark by flipping the sign of the temporary and decreasing by one
    void mark();

    /** \brief Assign to another expression, if a duplicate.
     * Check for equality up to a given depth */
    void assignIfDuplicate(const SXElement& scalar, int depth=1);

    /** \brief Assign the node to something, without invoking the deletion of the node,
     * if the count reaches 0 */
    SXNode* assignNoDelete(const SXElement& scalar);
    /// \endcond

    /** \brief SXElement nodes are not allowed to be null */
    inline bool isNull() {return false;}

#ifndef SWIG
  private:
    /// Pointer to node (SXElement is only a reference class)
    SXNode* node;

    /** \brief inline if-test */
    /// replaces the ternary conditional operator "?:", which cannot be overloaded
    friend SXElement if_else(const SXElement& cond, const SXElement& if_true,
                             const SXElement& if_false);
#endif // SWIG

  };

/// \cond INTERNAL
#ifndef SWIG
  // Template specializations
  template<>
  CASADI_CORE_EXPORT bool Matrix<SXElement>::__nonzero__() const;

  template<>
  class CASADI_CORE_EXPORT casadi_limits<SXElement>{
  public:
    static bool isZero(const SXElement& val);
    static bool isAlmostZero(const SXElement& val, double tol);
    static bool isOne(const SXElement& val);
    static bool isMinusOne(const SXElement& val);
    static bool isConstant(const SXElement& val);
    static bool isInteger(const SXElement& val);
    static bool isInf(const SXElement& val);
    static bool isMinusInf(const SXElement& val);
    static bool isNaN(const SXElement& val);

    static const SXElement zero;
    static const SXElement one;
    static const SXElement two;
    static const SXElement minus_one;
    static const SXElement nan;
    static const SXElement inf;
    static const SXElement minus_inf;
  };

#endif // SWIG
/// \endcond

  typedef std::vector<SXElement> SXElementVector;
  typedef std::vector<std::vector<SXElement> > SXElementVectorVector;
  typedef std::vector< std::vector<std::vector<SXElement> > > SXElementVectorVectorVector;
  typedef Matrix<SXElement> SX;
  typedef std::vector<Matrix<SXElement> > SXVector;
  typedef std::vector< std::vector<Matrix<SXElement> > > SXVectorVector;

  /// \cond INTERNAL
  typedef SX* SXPtr;
  typedef std::vector<SXPtr> SXPtrV;
  typedef std::vector<SXPtrV> SXPtrVV;
  /// \endcond

  // Specialize functions in GenericMatrix<SX> and SX
  template<> SX GenericMatrix<SX>::sym(const std::string& name, const Sparsity& sp);
  template<> bool SX::isRegular() const;
  template<> bool SX::isSmooth() const;
  template<> bool SX::isSymbolic() const;
  template<> bool SX::isSymbolicSparse() const;
  template<> double SX::getValue() const;
  template<> std::string SX::getName() const;
  template<> void SX::setMaxNumCallsInPrint(long num);
  template<> long SX::getMaxNumCallsInPrint();
  template<> void SX::setEqualityCheckingDepth(int eq_depth);
  template<> int SX::getEqualityCheckingDepth();

} // namespace casadi



#ifndef SWIG

/// \cond INTERNAL
// Template specialization
namespace casadi {
  template<> inline std::string matrixName<SXElement>() { return "SX"; }
} // namespace casadi
/// \endcond

namespace std {
  template<>
  class CASADI_CORE_EXPORT numeric_limits<casadi::SXElement>{
  public:
    static const bool is_specialized = true;
    static casadi::SXElement min() throw();
    static casadi::SXElement max() throw();
    static const int  digits = 0;
    static const int  digits10 = 0;
    static const bool is_signed = false;
    static const bool is_integer = false;
    static const bool is_exact = false;
    static const int radix = 0;
    static casadi::SXElement epsilon() throw();
    static casadi::SXElement round_error() throw();
    static const int  min_exponent = 0;
    static const int  min_exponent10 = 0;
    static const int  max_exponent = 0;
    static const int  max_exponent10 = 0;

    static const bool has_infinity = true;
    static const bool has_quiet_NaN = true;
    static const bool has_signaling_NaN = false;
    //    static const float_denorm_style has_denorm = denorm absent;
    static const bool has_denorm_loss = false;
    static casadi::SXElement infinity() throw();
    static casadi::SXElement quiet_NaN() throw();
    //    static SXElement signaling_NaN() throw();
    //    static SXElement denorm_min() throw();
    static const bool is_iec559 = false;
    static const bool is_bounded = false;
    static const bool is_modulo = false;

    static const bool traps = false;
    static const bool tinyness_before = false;
    static const float_round_style round_style = round_toward_zero;
  };
} //namespace std

/** \brief  The following functions needs the class so they cannot be included
 * in the beginning of the header */
#include "sx_node.hpp"

#endif // SWIG

#endif // CASADI_SX_ELEMENT_HPP
