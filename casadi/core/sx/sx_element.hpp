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

  /** SXElement is exposed only as an empty struct to SWIG */
#ifdef SWIG
  struct SXElement {};
#else // SWIG

  /** \brief The basic scalar symbolic class of CasADi
      \author Joel Andersson
      \date 2010-2014
  */
  class CASADI_EXPORT SXElement : public GenericExpression<SXElement>,
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
    void repr(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const;

    /// Print a description of the object
    void print(std::ostream &stream=CASADI_COUT, bool trailing_newline=true) const;

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
    SXElement zz_plus(const SXElement& y) const;
    SXElement zz_minus(const SXElement& y) const;
    SXElement zz_times(const SXElement& y) const;
    SXElement zz_rdivide(const SXElement& y) const;
    SXElement zz_lt(const SXElement& y) const;
    SXElement zz_le(const SXElement& y) const;
    SXElement zz_eq(const SXElement& y) const;
    SXElement zz_ne(const SXElement& y) const;
    SXElement __truediv__(const SXElement &y) const {return zz_rdivide(y);}
    SXElement zz_power(const SXElement& b) const;
    SXElement __constpow__(const SXElement& b) const;

    SXElement __mrdivide__(const SXElement& b) const {  return *this / b;}
    SXElement zz_mpower(const SXElement& b) const {return pow(*this, b);}

    // The following functions serves two purposes:
    // Numpy compatibility and to allow unambiguous access
    SXElement zz_mul(const SXElement& y) const { return zz_times(y);}
    SXElement zz_exp() const;
    SXElement zz_log() const;
    SXElement zz_sqrt() const;
    SXElement sq() const;
    SXElement zz_sin() const;
    SXElement zz_cos() const;
    SXElement zz_tan() const;
    SXElement zz_asin() const;
    SXElement zz_acos() const;
    SXElement zz_atan() const;
    SXElement zz_floor() const;
    SXElement zz_ceil() const;
    SXElement zz_mod(const SXElement &y) const;
    SXElement zz_erf() const;
    SXElement zz_erfinv() const;
    SXElement zz_abs() const;
    SXElement zz_min(const SXElement &y) const;
    SXElement zz_max(const SXElement &y) const;
    SXElement inv() const;
    SXElement zz_sinh() const;
    SXElement zz_cosh() const;
    SXElement zz_tanh() const;
    SXElement zz_asinh() const;
    SXElement zz_acosh() const;
    SXElement zz_atanh() const;
    SXElement zz_atan2(const SXElement &y) const;
    SXElement zz_log10() const;
    SXElement printme(const SXElement &y) const;
    SXElement zz_sign() const;
    SXElement __copysign__(const SXElement &y) const;
    SXElement constpow(const SXElement& y) const;
    SXElement zz_not() const;
    SXElement zz_and(const SXElement& y) const;
    SXElement zz_or(const SXElement& y) const;
    SXElement zz_if_else_zero(const SXElement& y) const;

    Matrix<SXElement> zz_min(const Matrix<SXElement>& b) const;
    Matrix<SXElement> zz_max(const Matrix<SXElement>& b) const;
    Matrix<SXElement> constpow(const Matrix<SXElement>& n) const;
    Matrix<SXElement> __copysign__(const Matrix<SXElement>& n) const;
    Matrix<SXElement> zz_atan2(const Matrix<SXElement>& b) const;
    bool zz_isEqual(const SXElement& scalar, int depth=0) const;
    SXElement zz_simplify() const;

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

  private:
    /// Pointer to node (SXElement is only a reference class)
    SXNode* node;

    /**
    \ingroup expression_tools
    @{
    */
    /** \brief inline if-test */
    /// replaces the ternary conditional operator "?:", which cannot be overloaded
    friend SXElement if_else(const SXElement& cond, const SXElement& if_true,
                             const SXElement& if_false);
    /** @} */
  };

  /// \cond INTERNAL
  // Template specializations
  template<>
  CASADI_EXPORT bool Matrix<SXElement>::__nonzero__() const;

  template<>
  class CASADI_EXPORT casadi_limits<SXElement>{
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

  // Specialize functions in GenericMatrix<SX> and SX
  template<> SX GenericMatrix<SX>::sym(const std::string& name, const Sparsity& sp);
  template<> bool SX::isRegular() const;
  template<> bool SX::isSmooth() const;
  template<> bool SX::isLeaf() const;
  template<> bool SX::isCommutative() const;
  template<> bool SX::isSymbolic() const;
  template<> bool SX::isSymbolicSparse() const;
  template<> double SX::getValue(int k) const;
  template<> int SX::getIntValue() const;
  template<> std::vector<double> SX::nonzeros() const;
  template<> std::vector<int> SX::nonzeros_int() const;
  template<> SX SX::getDep(int ch) const;
  template<> int SX::getNdeps() const;
  template<> std::string SX::getName() const;
  template<> void SX::setMaxNumCallsInPrint(long num);
  template<> long SX::getMaxNumCallsInPrint();
  template<> void SX::setEqualityCheckingDepth(int eq_depth);
  template<> int SX::getEqualityCheckingDepth();
  template<> long SX::getElementHash() const;
  template<> void SX::zz_expand(SX &weights, SX& terms) const;
  template<> SX SX::zz_pw_const(const SX &tval, const SX &val) const;
  template<> SX SX::zz_pw_lin(const SX &tval, const SX &val) const;
  template<> SX SX::zz_if_else(const SX &if_true,
                               const SX &if_false) const;
  template<> SX SX::zz_heaviside() const;
  template<> SX SX::zz_rectangle() const;
  template<> SX SX::zz_triangle() const;
  template<> SX SX::zz_ramp() const;
  template<> SX SX::zz_gauss_quadrature(const SX &x, const SX &a,
                                        const SX &b, int order,
                                        const SX& w) const;
  template<> SX SX::zz_simplify() const;
  template<> SX SX::zz_substitute(const SX& v, const SX& vdef) const;
  template<> std::vector<SX > SX::zz_substitute(const std::vector<SX >& ex,
                                                const std::vector<SX >& v,
                                                const std::vector<SX >& vdef);
  template<> void SX::zz_substituteInPlace(const std::vector<SX >& v,
                                           std::vector<SX >& vdef,
                                           std::vector<SX >& ex,
                                           bool reverse);
  template<> SX SX::zz_spy() const;
  template<> bool SX::zz_dependsOn(const SX &arg) const;
  template<> std::vector<SX > SX::zz_getSymbols() const;
  template<> std::vector<SX > SX::zz_getSymbols(const std::vector<SX >& e);
  template<> SX SX::zz_jacobian(const SX &arg) const;
  template<> SX SX::zz_gradient(const SX &arg) const;
  template<> SX SX::zz_tangent(const SX &arg) const;
  template<> SX SX::zz_hessian(const SX &arg) const;
  template<> void SX::zz_hessian(const SX &arg, SX &H, SX &g) const;
  template<> SX SX::zz_jacobianTimesVector(const SX &arg, const SX &v,
                                           bool transpose_jacobian) const;
  template<> SX SX::zz_taylor(const SX& x, const SX& a, int order) const;
  template<> SX SX::zz_mtaylor(const SX& x, const SX& a, int order) const;
  template<> SX SX::zz_mtaylor(const SX& x, const SX& a, int order,
                               const std::vector<int>& order_contributions) const;
  template<> int SX::zz_countNodes() const;
  template<> std::string
  SX::zz_getOperatorRepresentation(const std::vector<std::string>& args) const;
  template<> void SX::zz_extractShared(std::vector<SX >& ex,
                                       std::vector<SX >& v,
                                       std::vector<SX >& vdef,
                                       const std::string& v_prefix,
                                       const std::string& v_suffix);
  template<> void SX::zz_printCompact(std::ostream &stream) const;
  template<> SX SX::zz_poly_coeff(const SX&x) const;
  template<> SX SX::zz_poly_roots() const;
  template<> SX SX::zz_eig_symbolic() const;
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
  class CASADI_EXPORT numeric_limits<casadi::SXElement>{
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
