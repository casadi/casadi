/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#ifndef SX_HPP
#define SX_HPP

// exception class
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

namespace CasADi{

  /** \brief  forward declaration of Node and Matrix */
  class SXNode; // include will follow in the end

  /** \brief The basic scalar symbolic class of CasADi
      \author Joel Andersson 
      \date 2010
  */ 

#ifdef SWIG
#ifdef WITH_IMPLICITCONV
  %implicitconv SX;
#endif // WITH_IMPLICITCONV
#endif // SWIG


  class SX : public GenericExpression<SX>{
    friend class SXNode;
    friend class BinarySXNode;

  public:
    
    /// Constructors
    /** \brief Default constructor (not-a-number)

        Object is initialised as not-a-number.
    */
    SX();
    /** \brief Numerical constant constructor
        \param val Numerical value
    */
    SX(double val);
    
    /** \brief Symbolic constructor
         \param name Name of the symbol

        This is the name that wil be used by the "operator<<" and "toSTring" methods.
        The name is not used as identifier; you may construct distinct SX objects with non-unique names.
    */
    explicit SX(const std::string& name); // variable (must be explicit, otherwise 0/NULL would be ambigous)
    /** \brief Symbolic constructor
         \param Name of the symbol

        This is the name that wil be used by the "operator<<" and "toSTring" methods.
        The name is not used as identifier; you may construct distinct SX objects with non-unique names.
    */
#ifndef SWIG

    /// Create an expression from a node: extra dummy argument to avoid ambigousity for 0/NULL
    SX(SXNode* node, bool dummy);
    
    /** \brief Copy constructor */
    SX(const SX& scalar); // copy constructor

    /// Destructor
    ~SX();

    /// Create an object given a node
    static SX create(SXNode* node);
    
    // Assignment
    SX& operator=(const SX& scalar);
    SX& operator=(double scalar); // needed since otherwise both a = SX(double) and a = Matrix(double) would be ok

    // Convert to a 1-by-1 Matrix
    operator Matrix<SX>() const;
    
    /** \brief  print to stream */
    friend std::ostream& operator<<(std::ostream &stream, const SX &scalar);

    /** \brief  print to stream, limited */
#ifndef SWIG
    void print(std::ostream &stream, long& remaining_calls) const;
#endif // SWIG
    
    /** \brief  string representation (SWIG workaround) */
    std::string toString() const;
    
    /** \brief  Get a pointer to the node */
    SXNode* get() const; // note: constant pointer, not pointer to constant object! (to allow access to the counter)

    /** \brief  Access functions of the node */
    const SXNode* operator->() const;
    SXNode* operator->();
#endif // SWIG
    
    /** \brief  Perform operations by ID */
    static SX binary(int op, const SX& x, const SX& y);
    static SX unary(int op, const SX& x);
    
    /** \brief Check the truth value of this node
     * Introduced to catch bool(x) situations in python
     */
    bool __nonzero__() const;
    
    /** \brief check if this SX is a leaf of the SX graph
     *
     * An SX qualifies as leaf when it has no dependencies.
     */
    bool isLeaf() const;
    bool isConstant() const;
    bool isInteger() const;
    bool isSymbolic() const;
    bool hasDep() const;
    /** \brief Check wether a binary SX is commutative*/
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
    
    /** \brief Check if two nodes are equivalent up to a given depth. 
     *  Depth=0 checks if the expressions are identical, i.e. points to the same node.
     * 
     *  a = x*x
     *  b = x*x
     *
     *  a.isEqual(b,0)  will return false, but a.isEqual(b,1) will return true
     */
    bool isEqual(const SX& scalar, int depth=0) const;
    
    /** \brief Check if a value is always nonnegative (false negatives are allowed) */
    bool isNonNegative() const;
    
    double getValue() const;
    int getIntValue() const;
    SX getDep(int ch=0) const;
    
    /** \brief Check if the node is the sum of two equal expressions */
    bool isDoubled() const;
    
    /** \brief Get the number of dependencies of a binary SX */
    int getNdeps() const;
    
    /** \brief Returns a number that is unique for a given SXNode. 
     * If the SX does not point to any node, 0 is returned.
     */
    long __hash__() const;

    /** \brief  Negation */
    SX operator-() const;

    //  all binary operations
    SX __add__(const SX& y) const;
    SX __sub__(const SX& y) const;
    SX __mul__(const SX& y) const;
    SX __div__(const SX& y) const;
    SX __lt__(const SX& y) const;
    SX __le__(const SX& y) const;
    SX __eq__(const SX& y) const;
    SX __ne__(const SX& y) const;
    using GenericExpression<SX>::__gt__;
    using GenericExpression<SX>::__ge__;
    using GenericExpression<SX>::__mldivide__;
    SX __truediv__(const SX &y) const {return __div__(y);};
    SX __pow__(const SX& b) const;
    SX __constpow__(const SX& b) const;
    
    SX __mrdivide__(const SX& b) const{  return *this / b;}
    SX __mpower__(const SX& b) const {return (*this).__pow__(b);}
    SX trans() const{ return *this;}
    
    /// The following functions serves two purposes: Numpy compatibility and to allow unambigous access
    SX mul(const SX& y) const{ return __mul__(y);}
    SX exp() const;
    SX log() const;
    SX sqrt() const;
    SX sq() const;
    SX sin() const;
    SX cos() const;
    SX tan() const;
    SX arcsin() const;
    SX arccos() const;
    SX arctan() const;
    SX floor() const;
    SX ceil() const;
    SX erf() const;
    SX erfinv() const;
    SX fabs() const;
    SX fmin(const SX &y) const;
    SX fmax(const SX &y) const;
    SX inv() const;
    SX sinh() const;
    SX cosh() const;
    SX tanh() const;
    SX arcsinh() const;
    SX arccosh() const;
    SX arctanh() const;
    SX arctan2(const SX &y) const;
    SX log10() const;
    SX printme(const SX &y) const;
    SX sign() const;
    SX __copysign__(const SX &y) const;
    SX constpow(const SX& y) const;
    SX logic_not() const;
    SX logic_and(const SX& y) const;
    SX logic_or(const SX& y) const;
    SX if_else_zero(const SX& y) const;

    Matrix<SX> fmin(const Matrix<SX>& b) const;
    Matrix<SX> fmax(const Matrix<SX>& b) const;
    Matrix<SX> constpow(const Matrix<SX>& n) const;
    Matrix<SX> __copysign__(const Matrix<SX>& n) const;
    Matrix<SX> arctan2(const Matrix<SX>& b) const;
        
    // Get the temporary variable
    int getTemp() const;
    
    // Set the temporary variable
    void setTemp(int t);
    
    // Check if marked (i.e. temporary is negative)
    bool marked() const;
    
    // Mark by flipping the sign of the temporary and decreasing by one
    void mark();
    
    /** \brief Assign to another expression, if a duplicate. Check for equality up to a given depth */
    void assignIfDuplicate(const SX& scalar, int depth=1);
    
    /** \brief Set or reset the maximum number of calls to the printing function when printing an expression */
    static void setMaxNumCallsInPrint(long num=10000);

    /** \brief Get the maximum number of calls to the printing function when printing an expression */
    static long getMaxNumCallsInPrint();
    
    /** \brief Set or reset the depth to which equalities are being checked for simplifications */
    static void setEqualityCheckingDepth(int eq_depth=1);

    /** \brief Get the depth to which equalities are being checked for simplifications */
    static int getEqualityCheckingDepth();
    
    /** \brief Assign the node to something, without invoking the deletion of the node, if the count reaches 0 */
    SXNode* assignNoDelete(const SX& scalar);
    
    /** \brief SX nodes are not allowed to be null */
    inline bool isNull(){return false;}

#ifndef SWIG
  private:
    // Maximum number of calls
    static long max_num_calls_in_print_;
    
    // Depth when checking equalities
    static int eq_depth_;
    
    // Pointer to node (SX is only a reference class)
    SXNode* node;
    
    /** \brief inline if-test */
    friend SX if_else(const SX& cond, const SX& if_true, const SX& if_false); // replaces the ternary conditional operator "?:", which cannot be overloaded
#endif // SWIG

  };


#ifdef SWIG
  %extend SX {
    std::string __str__()  { return $self->toString(); }
    std::string __repr__() { return $self->toString(); }
    double __float__() { return $self->getValue();}
    int __int__() { return $self->getIntValue();}
  }
#endif // SWIG

#ifndef SWIG
  // Template specializations
  template<>
  bool __nonzero__(const SX& val);

  template<>
  class casadi_limits<SX>{
  public:
    static bool isZero(const SX& val);
    static bool isAlmostZero(const SX& val, double tol);
    static bool isOne(const SX& val);
    static bool isMinusOne(const SX& val);
    static bool isConstant(const SX& val);
    static bool isInteger(const SX& val);
    static bool isInf(const SX& val);
    static bool isMinusInf(const SX& val);
    static bool isNaN(const SX& val);

    static const SX zero;
    static const SX one;
    static const SX two;
    static const SX minus_one;
    static const SX nan;
    static const SX inf; 
    static const SX minus_inf;
  };

#endif // SWIG

  typedef std::vector<SX> SXVector;
  typedef std::vector<std::vector<SX> > SXVectorVector;
  typedef std::vector< std::vector<std::vector<SX> > > SXVectorVectorVector;
  typedef Matrix<SX> SXMatrix;
  typedef std::vector<Matrix<SX> > SXMatrixVector;
  typedef std::vector< std::vector<Matrix<SX> > > SXMatrixVectorVector;

  typedef SXMatrix* SXMatrixPtr;
  typedef std::vector<SXMatrixPtr> SXMatrixPtrV;
  typedef std::vector<SXMatrixPtrV> SXMatrixPtrVV;

} // namespace CasADi



#ifndef SWIG

// Template specialization
namespace CasADi{
  template<> inline const char* typeName<SX>() { return "SX"; }
} // namespace CasADi

namespace std{
  template<>
  class numeric_limits<CasADi::SX>{
  public:
    static const bool is_specialized = true;
    static CasADi::SX min() throw();
    static CasADi::SX max() throw();
    static const int  digits = 0;
    static const int  digits10 = 0;
    static const bool is_signed = false;
    static const bool is_integer = false;
    static const bool is_exact = false;
    static const int radix = 0;
    static CasADi::SX epsilon() throw();
    static CasADi::SX round_error() throw();
    static const int  min_exponent = 0;
    static const int  min_exponent10 = 0;
    static const int  max_exponent = 0;
    static const int  max_exponent10 = 0;
    
    static const bool has_infinity = true;
    static const bool has_quiet_NaN = true;
    static const bool has_signaling_NaN = false;
    //    static const float_denorm_style has_denorm = denorm absent;
    static const bool has_denorm_loss = false;
    static CasADi::SX infinity() throw();
    static CasADi::SX quiet_NaN() throw();
    //    static SX signaling_NaN() throw();
    //    static SX denorm_min() throw();
    static const bool is_iec559 = false;
    static const bool is_bounded = false;
    static const bool is_modulo = false;

    static const bool traps = false;
    static const bool tinyness_before = false;
    static const float_round_style round_style = round_toward_zero;
  };
} //namespace std

/** \brief  The following functions needs the class so they cannot be included in the beginning of the header */
#include "sx_node.hpp"

#endif // SWIG


#endif // SX_HPP
