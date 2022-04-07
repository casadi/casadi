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


#ifndef CASADI_CALCULUS_HPP
#define CASADI_CALCULUS_HPP

#include <iostream>
#include <string>
#include <cmath>
#include <limits>
#include <algorithm>
#include "casadi_common.hpp"

// Define pi if the compiler fails to do so

/// \cond INTERNAL

namespace casadi {
#ifndef SWIG
  /// Define pi
#ifdef M_PI
  const double pi = M_PI;
#else
  const double pi = 3.14159265358979323846;
#endif

  /// infinity
  const double inf = std::numeric_limits<double>::infinity();

  /// Not a number
  const double nan = std::numeric_limits<double>::quiet_NaN();

  /// Machine epsilon
  const double eps = std::numeric_limits<double>::epsilon();
#endif // SWIG

  /// Enum for quick access to any node
  enum class Operation : unsigned char {
    // Simple assignment
    OP_ASSIGN,

    // Standard unary and binary functions
    OP_ADD,  OP_SUB,  OP_MUL,  OP_DIV,
    OP_NEG,  OP_EXP,  OP_LOG,  OP_POW, OP_CONSTPOW,
    OP_SQRT,  OP_SQ,  OP_TWICE,
    OP_SIN,  OP_COS,  OP_TAN,
    OP_ASIN,  OP_ACOS,  OP_ATAN,
    OP_LT, OP_LE, OP_EQ, OP_NE, OP_NOT, OP_AND, OP_OR,
    OP_FLOOR,  OP_CEIL,  OP_FMOD, OP_FABS, OP_SIGN, OP_COPYSIGN, OP_IF_ELSE_ZERO,
    OP_ERF,  OP_FMIN,  OP_FMAX,
    OP_INV,
    OP_SINH,  OP_COSH,  OP_TANH,
    OP_ASINH, OP_ACOSH, OP_ATANH,
    OP_ATAN2,

    // Double constant
    OP_CONST,

    // Function input and output
    OP_INPUT, OP_OUTPUT,

    // Free parameter
    OP_PARAMETER,

    // Embedded function call
    OP_CALL,

    // Find first nonzero in a vector
    OP_FIND,

    // Find first nonzero in a vector
    OP_LOW,

    // Embedded function call in parallel
    OP_MAP,

    // Matrix multiplication
    OP_MTIMES,

    // Solve linear system of equations
    OP_SOLVE,

    // Matrix transpose
    OP_TRANSPOSE,

    // Matrix determinant
    OP_DETERMINANT,

    // Matrix inverse
    OP_INVERSE,

    // Inner product
    OP_DOT,

    // Bilinear form
    OP_BILIN,

    // Rank-1 update
    OP_RANK1,

    // Horizontal concatenation
    OP_HORZCAT,

    // Vertical concatenation of vectors
    OP_VERTCAT,

    // Diagonal concatenation
    OP_DIAGCAT,

    // Horizontal split
    OP_HORZSPLIT,

    // Vertical split of vectors
    OP_VERTSPLIT,

    // Diagonal split
    OP_DIAGSPLIT,

    // Reshape an expression
    OP_RESHAPE,

    // Submatrix reference
    OP_SUBREF,

    // Submatrix assignment
    OP_SUBASSIGN,

    // Nonzero reference
    OP_GETNONZEROS,

    // Parametric nonzero reference
    OP_GETNONZEROS_PARAM,

    // Nonzero addition
    OP_ADDNONZEROS,

    // parametric nonzero addition
    OP_ADDNONZEROS_PARAM,

    // Nonzero assignment
    OP_SETNONZEROS,

    // Parametric nonzero assignment
    OP_SETNONZEROS_PARAM,

    // Set sparse
    OP_PROJECT,

    // Assertion
    OP_ASSERTION,

    // Monitor
    OP_MONITOR,

    // Norms
    OP_NORM2, OP_NORM1, OP_NORMINF, OP_NORMF,

    // min/max
    OP_MMIN, OP_MMAX,

    // Horizontal repeat
    OP_HORZREPMAT,

    // Horizontal repeat sum
    OP_HORZREPSUM,

    OP_ERFINV,
    OP_PRINTME,
    OP_LIFT,

    OP_EINSTEIN,

    OP_BSPLINE,

    OP_CONVEXIFY,

    OP_LOG1P,

    OP_EXPM1,

    OP_HYPOT,
  
    OP_LOGSUMEXP,

    OP_INVALID
  };
  #define NUM_BUILT_IN_OPS (OP_INVALID)

  #define OP_

#ifndef SWIG

  ///@{
  /** \brief Enable using elementary numerical operations without std:: prefix */
  using std::isfinite;
  using std::sqrt;
  using std::sin;
  using std::cos;
  using std::tan;
  using std::atan;
  using std::asin;
  using std::acos;
  using std::sinh;
  using std::cosh;
  using std::tanh;
  using std::exp;
  using std::log;
  using std::log10;
  using std::abs;
  using std::fabs;
  using std::floor;
  using std::ceil;
  using std::pow;
  using std::fmod;
  using std::atan2;
  using std::erf;
  using std::fmin;
  using std::fmax;
  using std::fabs;
  using std::atanh;
  using std::asinh;
  using std::acosh;
  using std::isnan;
  using std::isinf;
  using std::log1p;
  using std::expm1;
  using std::hypot;
  ///@}

  ///@{
  // Implement "missing" operations

  /// Sign function, note that sign(nan) == nan
  inline double sign(double x) { return x<0 ? -1 : x>0 ? 1 : x;}
  ///@}

  ///@}

  ///@{
  /** \brief  CasADi additions */
  inline double simplify(double x) { return x;}
  inline double constpow(double x, double y) { return pow(x, y);}
  inline double printme(double x, double y) {
    std::ios::fmtflags f(uout().flags());
    uout() << "|> " << y << " : ";
    uout() << std::setprecision(std::numeric_limits<double>::digits10 + 1) << std::scientific;
    uout() << x << std::endl;
    uout().flags(f);
    return x;
  }
  inline bool is_equal(double x, double y, casadi_int depth=0) { return x==y;}

  #ifdef HAS_COPYSIGN
  using std::copysign;
  #else
  /// copysign function
  inline double copysign(double x, double y) { return y>=0 ? fabs(x) : -fabs(x);}
  #endif //HAS_COPYSIGN

  // Integer maximum and minimum
  template <typename T>
  inline T casadi_max(T x, T y) { return std::max(x, y);}
  template <typename T>
  inline T casadi_min(T x, T y) { return std::min(x, y);}

  /// Conditional assignment
  inline double if_else_zero(double x, double y) { return x==0 ? 0 : y;}
  inline double if_else(double x, double y, double z) { return x==0 ? z : y;}
#ifdef HAS_ERFINV
  using ::erfinv;
#else // HAS ERFINV
  inline double erfinv(double x) throw() {
    // Approximation found in Sourceforge and modified: Not very efficient
    if (x>=1) {
      return x==1 ? inf : nan;
    } else if (x<=-1) {
      return x==-1 ? -inf : nan;
    } else if (x<-0.7) {
      double z = sqrt(-log((1.0+x)/2.0));
      return -(((1.641345311*z+3.429567803)*z-1.624906493)*z-1.970840454)/
          ((1.637067800*z+3.543889200)*z+1.0);
    } else {
      double y;
      if (x<0.7) {
        double z = x*x;
        y = x*(((-0.140543331*z+0.914624893)*z-1.645349621)*z+0.886226899)/
            ((((-0.329097515*z+0.012229801)*z+1.442710462)*z-2.118377725)*z+1.0);
      } else {
        double z = sqrt(-log((1.0-x)/2.0));
        y = (((1.641345311*z+3.429567803)*z-1.624906493)*z-1.970840454)/
            ((1.637067800*z+3.543889200)*z+1.0);
      }

      //polish x to full accuracy
      y = y - (erf(y) - x) / (2.0/sqrt(pi) * exp(-y*y));
      y = y - (erf(y) - x) / (2.0/sqrt(pi) * exp(-y*y));
      return y;
    }
  }
#endif // HAS_ERFINV
  ///@}

  template<typename T>
  T twice(const T& x) {
    return x+x;
  }

  template<typename T>
  T sq(const T& x) {
    return x*x;
  }

  template<Operation I>
  struct UnaryOperation {
    /// Function evaluation
    template<typename T> static inline void fcn(const T& x, T& f);

    /// Partial derivatives
    template<typename T> static inline void der(const T& x, const T& f, T* d);
  };

  template<Operation I>
  struct BinaryOperation {
    /// Function evaluation
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) {
        UnaryOperation<I>::fcn(x, f);}

    /// Partial derivatives - binary function
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        UnaryOperation<I>::der(x, f, d); d[1]=0; }
  };

  template<Operation I>
  struct BinaryOperationE {
    /// Function evaluation
    template<typename T> static inline T fcn(const T& x, const T& y) {
      T ret;
      BinaryOperation<I>::fcn(x, y, ret);
      return ret;
    }
  };

  /// Calculate function and derivative
  template<Operation I>
  struct DerBinaryOperation {
    /// Perform the operation
    template<typename T> static inline void derf(const T& x, const T& y, T& f, T* d) {

      /** First save to temp since f might have the same address as x or y,
      * in which case it will be incorrect in the second call
      */
      T tmp;

      /// Evaluate the function
      BinaryOperation<I>::fcn(x, y, tmp);

      /// Evaluate the partial derivatives
      BinaryOperation<I>::der(x, y, tmp, d);

      /// Now save f
      f = tmp;
    }
  };

  /// Perform a binary operation on two scalars
  template<Operation I>
  struct BinaryOperationSS {
    /// Function evaluation
    template<typename T> static inline void fcn(const T& x, const T& y, T& f, casadi_int n) {
      BinaryOperation<I>::fcn(x, y, f);
    }

    /// Partial derivatives - binary function
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d,
        casadi_int n) {
      BinaryOperation<I>::der(x, y, f, d);
    }
  };


  /// Perform a binary operation on two vectors
  template<Operation I>
  struct BinaryOperationVV {
    /// Function evaluation
    template<typename T> static inline void fcn(const T* x, const T* y, T* f, casadi_int n) {
      for (casadi_int i=0; i<n; ++i) {
        BinaryOperation<I>::fcn(*x++, *y++, *f++);
      }
    }

    /// Partial derivatives - binary function
    template<typename T> static inline void der(const T* x, const T* y,
        const T* f, T* d, casadi_int n) {
      for (casadi_int i=0; i<n; ++i, d+=2) {
        BinaryOperation<I>::der(*x++, *y++, *f++, d);
      }
    }
  };

  /// Perform a binary operation on a vector and a scalar
  template<Operation I>
  struct BinaryOperationVS {
    /// Function evaluation
    template<typename T> static inline void fcn(const T* x, const T& y, T* f, casadi_int n) {
      for (casadi_int i=0; i<n; ++i) {
        BinaryOperation<I>::fcn(*x++, y, *f++);
      }
    }

    /// Partial derivatives - binary function
    template<typename T> static inline void der(const T* x, const T& y,
        const T* f, T* d, casadi_int n) {
      for (casadi_int i=0; i<n; ++i, d+=2) {
        BinaryOperation<I>::der(*x++, y, *f++, d);
      }
    }
  };

  /// Perform a binary operation on a scalar and a vector
  template<Operation I>
  struct BinaryOperationSV {
    /// Function evaluation
    template<typename T> static inline void fcn(const T& x, const T* y, T* f, casadi_int n) {
      for (casadi_int i=0; i<n; ++i) {
        BinaryOperation<I>::fcn(x, *y++, *f++);
      }
    }

    /// Partial derivatives - binary function
    template<typename T> static inline void der(const T& x, const T* y,
        const T* f, T* d, casadi_int n) {
      for (casadi_int i=0; i<n; ++i, d+=2) {
        BinaryOperation<I>::der(x, *y++, *f++, d);
      }
    }
  };

  ///@{
  /// Smoothness (by default true)
  template<Operation I> struct SmoothChecker { static const bool check=true;};
  template<>      struct SmoothChecker<Operation::OP_LT>{ static const bool check=false;};
  template<>      struct SmoothChecker<Operation::OP_LE>{ static const bool check=false;};
  template<>      struct SmoothChecker<Operation::OP_FLOOR>{ static const bool check=false;};
  template<>      struct SmoothChecker<Operation::OP_CEIL>{ static const bool check=false;};
  template<>      struct SmoothChecker<Operation::OP_FMOD>{ static const bool check=false;};
  template<>      struct SmoothChecker<Operation::OP_EQ>{ static const bool check=false;};
  template<>      struct SmoothChecker<Operation::OP_NE>{ static const bool check=false;};
  template<>      struct SmoothChecker<Operation::OP_SIGN>{ static const bool check=false;};
  template<>      struct SmoothChecker<Operation::OP_COPYSIGN>{ static const bool check=false;};
  template<>      struct SmoothChecker<Operation::OP_NOT>{ static const bool check=false;};
  template<>      struct SmoothChecker<Operation::OP_AND>{ static const bool check=false;};
  template<>      struct SmoothChecker<Operation::OP_OR>{ static const bool check=false;};
  template<>      struct SmoothChecker<Operation::OP_IF_ELSE_ZERO>{ static const bool check=false;};
  ///@}

  ///@{
  /// If evaluated with the first argument zero, is the result zero?
  template<Operation I> struct F0XChecker { static const bool check=false;};
  template<>      struct F0XChecker<Operation::OP_ASSIGN>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_MUL>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_DIV>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_NEG>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_POW>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_CONSTPOW>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_SQRT>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_SQ>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_TWICE>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_SIN>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_TAN>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_ATAN>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_ASIN>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_FLOOR>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_CEIL>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_FMOD>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_FABS>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_SIGN>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_COPYSIGN>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_ERF>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_SINH>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_TANH>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_ASINH>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_ATANH>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_ERFINV>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_AND>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_IF_ELSE_ZERO>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_LOG1P>{ static const bool check=true;};
  template<>      struct F0XChecker<Operation::OP_EXPM1>{ static const bool check=true;};
  ///@}

  ///@{
  /// If evaluated with the second argument zero, is the result zero?
  template<Operation I> struct FX0Checker { static const bool check=false;};
  template<>      struct FX0Checker<Operation::OP_MUL>{ static const bool check=true;};
  template<>      struct FX0Checker<Operation::OP_AND>{ static const bool check=true;};
  template<>      struct FX0Checker<Operation::OP_IF_ELSE_ZERO>{ static const bool check=true;};
  ///@}

  ///@{
  /// If evaluated with both arguments zero, is the result zero?
  template<Operation I> struct F00Checker {
    static const bool check=F0XChecker<I>::check || FX0Checker<I>::check;
  };
  template<>      struct F00Checker<Operation::OP_ADD>{ static const bool check=true;};
  template<>      struct F00Checker<Operation::OP_SUB>{ static const bool check=true;};
  template<>      struct F00Checker<Operation::OP_FMIN>{ static const bool check=true;};
  template<>      struct F00Checker<Operation::OP_FMAX>{ static const bool check=true;};
  template<>      struct F00Checker<Operation::OP_AND>{ static const bool check=true;};
  template<>      struct F00Checker<Operation::OP_OR>{ static const bool check=true;};
  template<>      struct F00Checker<Operation::OP_COPYSIGN>{ static const bool check=true;};
  template<>      struct F00Checker<Operation::OP_LT>{ static const bool check=true;};
  template<>      struct F00Checker<Operation::OP_HYPOT>{ static const bool check=true;};
  ///@}

  ///@{
  /// Is commutative
  template<Operation I> struct CommChecker { static const bool check=false;};
  template<>      struct CommChecker<Operation::OP_ADD>{ static const bool check=true;};
  template<>      struct CommChecker<Operation::OP_MUL>{ static const bool check=true;};
  template<>      struct CommChecker<Operation::OP_EQ>{ static const bool check=true;};
  template<>      struct CommChecker<Operation::OP_NE>{ static const bool check=true;};
  template<>      struct CommChecker<Operation::OP_AND>{ static const bool check=true;};
  template<>      struct CommChecker<Operation::OP_OR>{ static const bool check=true;};
  template<>      struct CommChecker<Operation::OP_HYPOT>{ static const bool check=true;};
  ///@}

  ///@{
  /// Always non-negative (false by default)
  template<Operation I> struct NonnegativeChecker { static const bool check=false;};
  template<>      struct NonnegativeChecker<Operation::OP_SQRT>{ static const bool check=true;};
  template<>      struct NonnegativeChecker<Operation::OP_SQ>{ static const bool check=true;};
  template<>      struct NonnegativeChecker<Operation::OP_EXP>{ static const bool check=true;};
  template<>      struct NonnegativeChecker<Operation::OP_LT>{ static const bool check=true;};
  template<>      struct NonnegativeChecker<Operation::OP_LE>{ static const bool check=true;};
  template<>      struct NonnegativeChecker<Operation::OP_EQ>{ static const bool check=true;};
  template<>      struct NonnegativeChecker<Operation::OP_NE>{ static const bool check=true;};
  template<>      struct NonnegativeChecker<Operation::OP_NOT>{ static const bool check=true;};
  template<>      struct NonnegativeChecker<Operation::OP_AND>{ static const bool check=true;};
  template<>      struct NonnegativeChecker<Operation::OP_OR>{ static const bool check=true;};
  template<>      struct NonnegativeChecker<Operation::OP_HYPOT>{ static const bool check=true;};
  ///@}

  ///@{
  /// Is the operation binary as opposed to unary
  template<Operation I> struct NargChecker { static const casadi_int check=1;};
  template<>      struct NargChecker<Operation::OP_ADD>{ static const casadi_int check=2;};
  template<>      struct NargChecker<Operation::OP_SUB>{ static const casadi_int check=2;};
  template<>      struct NargChecker<Operation::OP_MUL>{ static const casadi_int check=2;};
  template<>      struct NargChecker<Operation::OP_DIV>{ static const casadi_int check=2;};
  template<>      struct NargChecker<Operation::OP_POW>{ static const casadi_int check=2;};
  template<>      struct NargChecker<Operation::OP_CONSTPOW>{ static const casadi_int check=2;};
  template<>      struct NargChecker<Operation::OP_EQ>{ static const casadi_int check=2;};
  template<>      struct NargChecker<Operation::OP_NE>{ static const casadi_int check=2;};
  template<>      struct NargChecker<Operation::OP_AND>{ static const casadi_int check=2;};
  template<>      struct NargChecker<Operation::OP_OR>{ static const casadi_int check=2;};
  template<>      struct NargChecker<Operation::OP_FMIN>{ static const casadi_int check=2;};
  template<>      struct NargChecker<Operation::OP_FMAX>{ static const casadi_int check=2;};
  template<>      struct NargChecker<Operation::OP_PRINTME>{ static const casadi_int check=2;};
  template<>      struct NargChecker<Operation::OP_ATAN2>{ static const casadi_int check=2;};
  template<>      struct NargChecker<Operation::OP_IF_ELSE_ZERO>{ static const casadi_int check=2;};
  template<>      struct NargChecker<Operation::OP_FMOD>{ static const casadi_int check=2;};
  template<>      struct NargChecker<Operation::OP_COPYSIGN>{ static const casadi_int check=2;};
  template<>      struct NargChecker<Operation::OP_CONST>{ static const casadi_int check=0;};
  template<>      struct NargChecker<Operation::OP_PARAMETER>{ static const casadi_int check=0;};
  template<>      struct NargChecker<Operation::OP_INPUT>{ static const casadi_int check=0;};
  template<>      struct NargChecker<Operation::OP_HYPOT>{ static const casadi_int check=2;};
  ///@}

  /// Simple assignment
  template<>
  struct UnaryOperation<Operation::OP_ASSIGN>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = x;}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0] = 1; }
  };

  /// Addition
  template<>
  struct BinaryOperation<Operation::OP_ADD>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = x+y;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=d[1]=1;}
  };

  /// Subtraction
  template<>
  struct BinaryOperation<Operation::OP_SUB>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = x-y;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=1; d[1]=-1;}
  };

  /// Multiplication
  template<>
  struct BinaryOperation<Operation::OP_MUL>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = x*y;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=y; d[1]=x;}
  };

  /// Division
  template<>
  struct BinaryOperation<Operation::OP_DIV>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = x/y;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=1/y; d[1]=-f/y;}
  };

  /// Negation
  template<>
  struct UnaryOperation<Operation::OP_NEG>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = -x;}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0]=-1;}
  };

  /// Natural exponent
  template<>
  struct UnaryOperation<Operation::OP_EXP>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = exp(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0]=f;}
  };

  /// Natural logarithm
  template<>
  struct UnaryOperation<Operation::OP_LOG>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = log(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0]=1/x;}
  };

  /// Power, defined only for x>=0
  template<>
  struct BinaryOperation<Operation::OP_POW>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = pow(x, y);}
    // See issue #104 why d[0] is no longer y*f/x
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=y*pow(x, y-1); d[1]=log(x)*f;}
  };

  /// Power, defined only for y constant
  template<>
  struct BinaryOperation<Operation::OP_CONSTPOW>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = pow(x, y);}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=y*pow(x, y-1); d[1]=0;}
  };

  /// Square root
  template<>
  struct UnaryOperation<Operation::OP_SQRT>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = sqrt(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0]=1/(twice(f));}
  };

  /// Square
  template<>
  struct UnaryOperation<Operation::OP_SQ>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = sq(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0]=twice(x);}
  };

  /// Times two
  template<>
  struct UnaryOperation<Operation::OP_TWICE>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = 2.*x;}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0] = 2; }
  };

  /// Sine
  template<>
  struct UnaryOperation<Operation::OP_SIN>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = sin(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0]=cos(x);}
  };

  /// Cosine
  template<>
  struct UnaryOperation<Operation::OP_COS>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = cos(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0]=-sin(x);}
  };

  /// Tangent
  template<>
  struct UnaryOperation<Operation::OP_TAN>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = tan(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d)
    { d[0] = 1/sq(cos(x));}
  };

  /// Arcus sine
  template<>
  struct UnaryOperation<Operation::OP_ASIN>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = asin(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0]=1/sqrt(1-x*x);}
  };

  /// Arcus cosine
  template<>
  struct UnaryOperation<Operation::OP_ACOS>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = acos(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d)
    { d[0]=-1/sqrt(1-x*x);}
  };

  /// Arcus tangent
  template<>
  struct UnaryOperation<Operation::OP_ATAN>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = atan(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0] = 1/(1+x*x);}
  };

  /// Less than
  template<>
  struct BinaryOperation<Operation::OP_LT>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = x < y;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=d[1]=0;}
  };

  /// Less or equal to
  template<>
  struct BinaryOperation<Operation::OP_LE>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = x <= y;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=d[1]=0;}
  };

  /// Floor function
  template<>
  struct UnaryOperation<Operation::OP_FLOOR>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = floor(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0] = 0;}
  };

  /// Ceil function
  template<>
  struct UnaryOperation<Operation::OP_CEIL>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = ceil(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0] = 0;}
  };

  /// Remainder of division
  template<>
  struct BinaryOperation<Operation::OP_FMOD>{
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = fmod(x, y);}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
      d[0]=1; d[1]=(f-x)/y;}
  };

  /// Equal to
  template<>
  struct BinaryOperation<Operation::OP_EQ>{
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = x==y;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=d[1]=0;}
  };

  /// Not equal to
  template<>
  struct BinaryOperation<Operation::OP_NE>{
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = x!=y;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=d[1]=0;}
  };

  /// Logical not
  template<>
  struct UnaryOperation<Operation::OP_NOT>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = !x;}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0] = 0;}
  };

  /// Logical and
  template<>
  struct BinaryOperation<Operation::OP_AND>{
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = x && y;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=d[1]=0;}
  };

  /// Logical or
  template<>
  struct BinaryOperation<Operation::OP_OR>{
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = x || y;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=d[1]=0;}
  };

  /// Error function
  template<>
  struct UnaryOperation<Operation::OP_ERF>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = erf(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) {
        d[0] = (2/sqrt(pi))*exp(-x*x);}
  };

  /// Absolute value
  template<>
  struct UnaryOperation<Operation::OP_FABS>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = fabs(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) {
        d[0]=sign(x);}
  };

  /// Sign
  template<>
  struct UnaryOperation<Operation::OP_SIGN>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = sign(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0]=0;}
  };

  /// Copysign
  template<>
  struct BinaryOperation<Operation::OP_COPYSIGN>{
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = copysign(x, y);}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        T e = 1; d[0]=copysign(e, y); d[1]=0;}
  };

  /// Minimum
  template<>
  struct BinaryOperation<Operation::OP_FMIN>{
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = fmin(x, y);}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
      T a = x<=y;
      T b = y<=x;
      T c = a+b;
      d[0]=a/c; d[1]=b/c;}
  };

  /// Maximum
  template<>
  struct BinaryOperation<Operation::OP_FMAX>{
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = fmax(x, y);}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
      T a = y<=x;
      T b = x<=y;
      T c = a+b;
      d[0]=a/c; d[1]=b/c;}
  };

  /// Element-wise inverse
  template<>
  struct UnaryOperation<Operation::OP_INV>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = 1./x;}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0] = -f*f; }
  };

  /// Hyperbolic sine
  template<>
  struct UnaryOperation<Operation::OP_SINH>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = sinh(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0] = cosh(x); }
  };

  /// Hyperbolic cosine
  template<>
  struct UnaryOperation<Operation::OP_COSH>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = cosh(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0] = sinh(x); }
  };

  /// Hyperbolic tangent
  template<>
  struct UnaryOperation<Operation::OP_TANH>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = tanh(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0] = 1-f*f; }
  };

  /// Inverse hyperbolic sine
  template<>
  struct UnaryOperation<Operation::OP_ASINH>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = asinh(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) {
        d[0] = 1/sqrt(1+x*x); }
  };

  /// Inverse hyperbolic cosine
  template<>
  struct UnaryOperation<Operation::OP_ACOSH>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = acosh(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) {
        d[0] = 1/sqrt(x-1)/sqrt(x+1); }
  };

  /// Inverse hyperbolic tangent
  template<>
  struct UnaryOperation<Operation::OP_ATANH>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = atanh(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0] = 1/(1-x*x); }
  };

  /// Inverse of error function
  template<>
  struct UnaryOperation<Operation::OP_ERFINV>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = erfinv(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) {
        d[0] = (sqrt(pi)/2)*exp(f*f); }
  };

  /// log1p(x) = log(1+x)
  template<>
  struct UnaryOperation<Operation::OP_LOG1P>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = log1p(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) {
        d[0] = 1/(1+x);}
  };

  /// expm1(x) = exp(x)-1
  template<>
  struct UnaryOperation<Operation::OP_EXPM1>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = expm1(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) {
        d[0] = exp(x); }
  };

  /// Identity operator with the side effect of printing
  template<>
  struct BinaryOperation<Operation::OP_PRINTME>{
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) {f = printme(x, y); }
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=1; d[1]=0;}
  };

  /// Arctan2
  template<>
  struct BinaryOperation<Operation::OP_ATAN2>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = atan2(x, y);}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        T t = x*x+y*y; d[0]=y/t; d[1]=-x/t;}
  };

  /// Conditional assignment
  template<>
  struct BinaryOperation<Operation::OP_IF_ELSE_ZERO>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) {
        f = if_else_zero(x, y);}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=0; d[1]=x;}
  };

  /// Inverse of error function
  template<>
  struct BinaryOperation<Operation::OP_LIFT>{
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = x;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0] = 1; d[1] = 0; }
  };

  /// hypot(x,y) = sqrt(x^2+y^2)
  template<>
  struct BinaryOperation<Operation::OP_HYPOT>{
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = hypot(x, y);}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0] = x/f; d[1] = y/f; }
  };

  template<template<Operation> class F, typename T>
  T operation_getter(Operation op) {
    switch (static_cast<Operation>(op)) {
    case Operation::OP_ASSIGN:        return F<Operation::OP_ASSIGN>::check;
    case Operation::OP_ADD:           return F<Operation::OP_ADD>::check;
    case Operation::OP_SUB:           return F<Operation::OP_SUB>::check;
    case Operation::OP_MUL:           return F<Operation::OP_MUL>::check;
    case Operation::OP_DIV:           return F<Operation::OP_DIV>::check;
    case Operation::OP_NEG:           return F<Operation::OP_NEG>::check;
    case Operation::OP_EXP:           return F<Operation::OP_EXP>::check;
    case Operation::OP_LOG:           return F<Operation::OP_LOG>::check;
    case Operation::OP_POW:           return F<Operation::OP_POW>::check;
    case Operation::OP_CONSTPOW:      return F<Operation::OP_CONSTPOW>::check;
    case Operation::OP_SQRT:          return F<Operation::OP_SQRT>::check;
    case Operation::OP_SQ:            return F<Operation::OP_SQ>::check;
    case Operation::OP_TWICE:         return F<Operation::OP_TWICE>::check;
    case Operation::OP_SIN:           return F<Operation::OP_SIN>::check;
    case Operation::OP_COS:           return F<Operation::OP_COS>::check;
    case Operation::OP_TAN:           return F<Operation::OP_TAN>::check;
    case Operation::OP_ASIN:          return F<Operation::OP_ASIN>::check;
    case Operation::OP_ACOS:          return F<Operation::OP_ACOS>::check;
    case Operation::OP_ATAN:          return F<Operation::OP_ATAN>::check;
    case Operation::OP_LT:            return F<Operation::OP_LT>::check;
    case Operation::OP_LE:            return F<Operation::OP_LE>::check;
    case Operation::OP_EQ:            return F<Operation::OP_EQ>::check;
    case Operation::OP_NE:            return F<Operation::OP_NE>::check;
    case Operation::OP_NOT:           return F<Operation::OP_NOT>::check;
    case Operation::OP_AND:           return F<Operation::OP_AND>::check;
    case Operation::OP_OR:            return F<Operation::OP_OR>::check;
    case Operation::OP_FLOOR:         return F<Operation::OP_FLOOR>::check;
    case Operation::OP_CEIL:          return F<Operation::OP_CEIL>::check;
    case Operation::OP_FMOD:          return F<Operation::OP_FMOD>::check;
    case Operation::OP_FABS:          return F<Operation::OP_FABS>::check;
    case Operation::OP_SIGN:          return F<Operation::OP_SIGN>::check;
    case Operation::OP_COPYSIGN:      return F<Operation::OP_COPYSIGN>::check;
    case Operation::OP_IF_ELSE_ZERO:  return F<Operation::OP_IF_ELSE_ZERO>::check;
    case Operation::OP_ERF:           return F<Operation::OP_ERF>::check;
    case Operation::OP_FMIN:          return F<Operation::OP_FMIN>::check;
    case Operation::OP_FMAX:          return F<Operation::OP_FMAX>::check;
    case Operation::OP_INV:           return F<Operation::OP_INV>::check;
    case Operation::OP_SINH:          return F<Operation::OP_SINH>::check;
    case Operation::OP_COSH:          return F<Operation::OP_COSH>::check;
    case Operation::OP_TANH:          return F<Operation::OP_TANH>::check;
    case Operation::OP_ASINH:         return F<Operation::OP_ASINH>::check;
    case Operation::OP_ACOSH:         return F<Operation::OP_ACOSH>::check;
    case Operation::OP_ATANH:         return F<Operation::OP_ATANH>::check;
    case Operation::OP_ATAN2:         return F<Operation::OP_ATAN2>::check;
    case Operation::OP_CONST:         return F<Operation::OP_CONST>::check;
    case Operation::OP_INPUT:         return F<Operation::OP_INPUT>::check;
    case Operation::OP_OUTPUT:        return F<Operation::OP_OUTPUT>::check;
    case Operation::OP_PARAMETER:     return F<Operation::OP_PARAMETER>::check;
    case Operation::OP_CALL:          return F<Operation::OP_CALL>::check;
    case Operation::OP_FIND:          return F<Operation::OP_FIND>::check;
    case Operation::OP_LOW:           return F<Operation::OP_LOW>::check;
    case Operation::OP_MAP:           return F<Operation::OP_MAP>::check;
    case Operation::OP_MTIMES:        return F<Operation::OP_MTIMES>::check;
    case Operation::OP_SOLVE:         return F<Operation::OP_SOLVE>::check;
    case Operation::OP_TRANSPOSE:     return F<Operation::OP_TRANSPOSE>::check;
    case Operation::OP_DETERMINANT:   return F<Operation::OP_DETERMINANT>::check;
    case Operation::OP_INVERSE:       return F<Operation::OP_INVERSE>::check;
    case Operation::OP_DOT:           return F<Operation::OP_DOT>::check;
    case Operation::OP_BILIN:         return F<Operation::OP_BILIN>::check;
    case Operation::OP_RANK1:         return F<Operation::OP_RANK1>::check;
    case Operation::OP_HORZCAT:       return F<Operation::OP_HORZCAT>::check;
    case Operation::OP_VERTCAT:       return F<Operation::OP_VERTCAT>::check;
    case Operation::OP_DIAGCAT:       return F<Operation::OP_DIAGCAT>::check;
    case Operation::OP_HORZSPLIT:     return F<Operation::OP_HORZSPLIT>::check;
    case Operation::OP_VERTSPLIT:     return F<Operation::OP_VERTSPLIT>::check;
    case Operation::OP_DIAGSPLIT:     return F<Operation::OP_DIAGSPLIT>::check;
    case Operation::OP_RESHAPE:       return F<Operation::OP_RESHAPE>::check;
    case Operation::OP_SUBREF:        return F<Operation::OP_SUBREF>::check;
    case Operation::OP_SUBASSIGN:     return F<Operation::OP_SUBASSIGN>::check;
    case Operation::OP_GETNONZEROS:   return F<Operation::OP_GETNONZEROS>::check;
    case Operation::OP_GETNONZEROS_PARAM:   return F<Operation::OP_GETNONZEROS_PARAM>::check;
    case Operation::OP_ADDNONZEROS:   return F<Operation::OP_ADDNONZEROS>::check;
    case Operation::OP_ADDNONZEROS_PARAM:   return F<Operation::OP_ADDNONZEROS>::check;
    case Operation::OP_SETNONZEROS:   return F<Operation::OP_SETNONZEROS>::check;
    case Operation::OP_SETNONZEROS_PARAM:   return F<Operation::OP_SETNONZEROS>::check;
    case Operation::OP_PROJECT:       return F<Operation::OP_PROJECT>::check;
    case Operation::OP_ASSERTION:     return F<Operation::OP_ASSERTION>::check;
    case Operation::OP_MONITOR:       return F<Operation::OP_MONITOR>::check;
    case Operation::OP_NORM2:         return F<Operation::OP_NORM2>::check;
    case Operation::OP_NORM1:         return F<Operation::OP_NORM1>::check;
    case Operation::OP_NORMINF:       return F<Operation::OP_NORMINF>::check;
    case Operation::OP_NORMF:         return F<Operation::OP_NORMF>::check;
    case Operation::OP_MMIN:          return F<Operation::OP_MMIN>::check;
    case Operation::OP_MMAX:          return F<Operation::OP_MMAX>::check;
    case Operation::OP_HORZREPMAT:    return F<Operation::OP_HORZREPMAT>::check;
    case Operation::OP_HORZREPSUM:    return F<Operation::OP_HORZREPSUM>::check;
    case Operation::OP_ERFINV:        return F<Operation::OP_ERFINV>::check;
    case Operation::OP_PRINTME:       return F<Operation::OP_PRINTME>::check;
    case Operation::OP_LIFT:          return F<Operation::OP_LIFT>::check;
    case Operation::OP_EINSTEIN:      return F<Operation::OP_EINSTEIN>::check;
    case Operation::OP_BSPLINE:       return F<Operation::OP_BSPLINE>::check;
    case Operation::OP_CONVEXIFY:     return F<Operation::OP_CONVEXIFY>::check;
    case Operation::OP_LOG1P:         return F<Operation::OP_LOG1P>::check;
    case Operation::OP_EXPM1:         return F<Operation::OP_EXPM1>::check;
    case Operation::OP_HYPOT:         return F<Operation::OP_HYPOT>::check;
    case Operation::OP_LOGSUMEXP:     return F<Operation::OP_LOGSUMEXP>::check;
    }
    return T();
  }

  template<template<Operation> class F>
  bool operation_checker(Operation op) {
    return operation_getter<F, bool>(op);
  }

  /// Easy access to all the functions for a particular type
  template<typename T>
  struct casadi_math {

    /** \brief Evaluate a built in function (scalar-scalar) */
    static inline void fun(Operation op, const T& x, const T& y, T& f);

    /** \brief Evaluate a built in function (vector-vector) */
    static inline void fun(Operation op, const T* x, const T* y, T* f, casadi_int n);

    /** \brief Evaluate a built in function (vector-scalar) */
    static inline void fun(Operation op, const T* x, const T& y, T* f, casadi_int n);

    /** \brief Evaluate a built in function (scalar-vector) */
    static inline void fun(Operation op, const T& x, const T* y, T* f, casadi_int n);

    /** \brief Evaluate a built in derivative function */
    static inline void der(Operation op, const T& x, const T& y, const T& f, T* d);

    /** \brief Evaluate the function and the derivative function */
    static inline void derF(Operation op, const T& x, const T& y, T& f, T* d);

    /** \brief Is binary operation? */
    static inline bool is_binary(Operation op);

    /** \brief Is unary operation? */
    static inline bool is_unary(Operation op);

    /** \brief Number of dependencies */
    static inline casadi_int ndeps(Operation op);

    /** \brief Print */
    static inline std::string print(Operation op, const std::string& x,
                             const std::string& y);
    static inline std::string print(Operation op, const std::string& x);
    static inline std::string name(Operation op);
    static inline std::string pre(Operation op);
    static inline std::string sep(Operation op);
    static inline std::string post(Operation op);
  };

  /// Specialize the class so that it can be used with integer type
  template<>
  struct casadi_math<casadi_int>{

    /** \brief Evaluate a built in function */
    static inline void fun(Operation op, const casadi_int& x,
        const casadi_int& y, casadi_int& f) {
      double ff(0);
      casadi_math<double>::fun(op, static_cast<double>(x), static_cast<double>(y), ff);
      f = static_cast<casadi_int>(ff);
    }

    static inline void fun(Operation op, const casadi_int* x, const casadi_int* y,
        casadi_int* f, casadi_int n) {
      for (casadi_int i=0; i<n; ++i) {
        double ff(0);
        casadi_math<double>::fun(op, static_cast<double>(*x++), static_cast<double>(*y++), ff);
        *f++ = static_cast<casadi_int>(ff);
      }
    }

    static inline void fun(Operation op, const casadi_int* x, const casadi_int& y,
        casadi_int* f, casadi_int n) {
      for (casadi_int i=0; i<n; ++i) {
        double ff;
        casadi_math<double>::fun(op, static_cast<double>(*x++), static_cast<double>(y), ff);
        *f++ = static_cast<casadi_int>(ff);
      }
    }

    static inline void fun(Operation op, const casadi_int& x, const casadi_int* y,
        casadi_int* f, casadi_int n) {
      for (casadi_int i=0; i<n; ++i) {
        double ff;
        casadi_math<double>::fun(op, static_cast<double>(x), static_cast<double>(*y++), ff);
        *f++ = static_cast<casadi_int>(ff);
      }
    }

    /** \brief Evaluate a built in derivative function */
    static inline void der(Operation op, const casadi_int& x, const casadi_int& y,
        const casadi_int& f, casadi_int* d) {
      double d_real[2] = {static_cast<double>(d[0]), static_cast<double>(d[1])};
      casadi_math<double>::der(op, static_cast<double>(x), static_cast<double>(y),
                               static_cast<double>(f), d_real);
      d[0] = static_cast<casadi_int>(d_real[0]);
      d[1] = static_cast<casadi_int>(d_real[1]);
    }

    /** \brief Evaluate the function and the derivative function */
    static inline void derF(Operation op, const casadi_int& x, const casadi_int& y,
        casadi_int& f, casadi_int* d) {
      double d_real[2] = {static_cast<double>(d[0]), static_cast<double>(d[1])};
      double f_real = static_cast<double>(f);
      casadi_math<double>::derF(op, static_cast<double>(x), static_cast<double>(y), f_real, d_real);
      f = static_cast<casadi_int>(f_real);
      d[0] = static_cast<casadi_int>(d_real[0]);
      d[1] = static_cast<casadi_int>(d_real[1]);
    }

    /** \brief Number of dependencies */
    static inline casadi_int ndeps(Operation op) {
      return casadi_math<double>::ndeps(op);
    }

    /** \brief Print */
    static inline std::string print(Operation op, const std::string& x,
                                    const std::string& y) {
      return casadi_math<double>::print(op, x, y);
    }
    static inline std::string print(Operation op, const std::string& x) {
      return casadi_math<double>::print(op, x);
    }
    static inline std::string pre(Operation op) {
      return casadi_math<double>::pre(op);
    }
    static inline std::string name(Operation op) {
      return casadi_math<double>::name(op);
    }
    static inline std::string sep(Operation op) {
      return casadi_math<double>::sep(op);
    }
    static inline std::string post(Operation op) {
      return casadi_math<double>::post(op);
    }
  };

  // Template implementations

  template<typename T>
  inline void casadi_math<T>::fun(Operation op, const T& x, const T& y, T& f) {
    // NOTE: We define the implementation in a preprocessor macro to be able to force inlining,
    //  and to allow extensions in the VM
#define CASADI_MATH_FUN_BUILTIN_GEN(CNAME, X, Y, F, N)                  \
    case Operation::OP_ASSIGN:    CNAME<Operation::OP_ASSIGN>::fcn(X, Y, F, N);        break; \
  case Operation::OP_ADD:       CNAME<Operation::OP_ADD>::fcn(X, Y, F, N);           break;   \
  case Operation::OP_SUB:       CNAME<Operation::OP_SUB>::fcn(X, Y, F, N);           break;   \
  case Operation::OP_MUL:       CNAME<Operation::OP_MUL>::fcn(X, Y, F, N);           break;   \
  case Operation::OP_DIV:       CNAME<Operation::OP_DIV>::fcn(X, Y, F, N);           break;   \
  case Operation::OP_NEG:       CNAME<Operation::OP_NEG>::fcn(X, Y, F, N);           break;   \
  case Operation::OP_EXP:       CNAME<Operation::OP_EXP>::fcn(X, Y, F, N);           break;   \
  case Operation::OP_LOG:       CNAME<Operation::OP_LOG>::fcn(X, Y, F, N);           break;   \
  case Operation::OP_POW:       CNAME<Operation::OP_POW>::fcn(X, Y, F, N);           break;   \
  case Operation::OP_CONSTPOW:  CNAME<Operation::OP_CONSTPOW>::fcn(X, Y, F, N);      break;   \
  case Operation::OP_SQRT:      CNAME<Operation::OP_SQRT>::fcn(X, Y, F, N);          break;   \
  case Operation::OP_SQ:        CNAME<Operation::OP_SQ>::fcn(X, Y, F, N);            break;   \
  case Operation::OP_TWICE:     CNAME<Operation::OP_TWICE>::fcn(X, Y, F, N);         break;   \
  case Operation::OP_SIN:       CNAME<Operation::OP_SIN>::fcn(X, Y, F, N);           break;   \
  case Operation::OP_COS:       CNAME<Operation::OP_COS>::fcn(X, Y, F, N);           break;   \
  case Operation::OP_TAN:       CNAME<Operation::OP_TAN>::fcn(X, Y, F, N);           break;   \
  case Operation::OP_ASIN:      CNAME<Operation::OP_ASIN>::fcn(X, Y, F, N);          break;   \
  case Operation::OP_ACOS:      CNAME<Operation::OP_ACOS>::fcn(X, Y, F, N);          break;   \
  case Operation::OP_ATAN:      CNAME<Operation::OP_ATAN>::fcn(X, Y, F, N);          break;   \
  case Operation::OP_LT:        CNAME<Operation::OP_LT>::fcn(X, Y, F, N);            break;   \
  case Operation::OP_LE:        CNAME<Operation::OP_LE>::fcn(X, Y, F, N);            break;   \
  case Operation::OP_EQ:        CNAME<Operation::OP_EQ>::fcn(X, Y, F, N);            break;   \
  case Operation::OP_NE:        CNAME<Operation::OP_NE>::fcn(X, Y, F, N);            break;   \
  case Operation::OP_NOT:       CNAME<Operation::OP_NOT>::fcn(X, Y, F, N);           break;   \
  case Operation::OP_AND:       CNAME<Operation::OP_AND>::fcn(X, Y, F, N);           break;   \
  case Operation::OP_OR:        CNAME<Operation::OP_OR>::fcn(X, Y, F, N);            break;   \
  case Operation::OP_IF_ELSE_ZERO: CNAME<Operation::OP_IF_ELSE_ZERO>::fcn(X, Y, F, N); break; \
  case Operation::OP_FLOOR:     CNAME<Operation::OP_FLOOR>::fcn(X, Y, F, N);         break;   \
  case Operation::OP_CEIL:      CNAME<Operation::OP_CEIL>::fcn(X, Y, F, N);          break;   \
  case Operation::OP_FMOD:      CNAME<Operation::OP_FMOD>::fcn(X, Y, F, N);          break;   \
  case Operation::OP_FABS:      CNAME<Operation::OP_FABS>::fcn(X, Y, F, N);          break;   \
  case Operation::OP_SIGN:      CNAME<Operation::OP_SIGN>::fcn(X, Y, F, N);          break;   \
  case Operation::OP_COPYSIGN:  CNAME<Operation::OP_COPYSIGN>::fcn(X, Y, F, N);      break;   \
  case Operation::OP_ERF:       CNAME<Operation::OP_ERF>::fcn(X, Y, F, N);           break;   \
  case Operation::OP_FMIN:      CNAME<Operation::OP_FMIN>::fcn(X, Y, F, N);          break;   \
  case Operation::OP_FMAX:      CNAME<Operation::OP_FMAX>::fcn(X, Y, F, N);          break;   \
  case Operation::OP_INV:       CNAME<Operation::OP_INV>::fcn(X, Y, F, N);           break;   \
  case Operation::OP_SINH:      CNAME<Operation::OP_SINH>::fcn(X, Y, F, N);          break;   \
  case Operation::OP_COSH:      CNAME<Operation::OP_COSH>::fcn(X, Y, F, N);          break;   \
  case Operation::OP_TANH:      CNAME<Operation::OP_TANH>::fcn(X, Y, F, N);          break;   \
  case Operation::OP_ASINH:     CNAME<Operation::OP_ASINH>::fcn(X, Y, F, N);         break;   \
  case Operation::OP_ACOSH:     CNAME<Operation::OP_ACOSH>::fcn(X, Y, F, N);         break;   \
  case Operation::OP_ATANH:     CNAME<Operation::OP_ATANH>::fcn(X, Y, F, N);         break;   \
  case Operation::OP_ATAN2:     CNAME<Operation::OP_ATAN2>::fcn(X, Y, F, N);         break;   \
  case Operation::OP_ERFINV:    CNAME<Operation::OP_ERFINV>::fcn(X, Y, F, N);        break;   \
  case Operation::OP_LIFT:      CNAME<Operation::OP_LIFT>::fcn(X, Y, F, N);          break;   \
  case Operation::OP_PRINTME:   CNAME<Operation::OP_PRINTME>::fcn(X, Y, F, N);       break;   \
  case Operation::OP_LOG1P:     CNAME<Operation::OP_LOG1P>::fcn(X, Y, F, N);       break;   \
  case Operation::OP_EXPM1:     CNAME<Operation::OP_EXPM1>::fcn(X, Y, F, N);       break;   \
  case Operation::OP_HYPOT:     CNAME<Operation::OP_HYPOT>::fcn(X, Y, F, N);       break;

#define CASADI_MATH_FUN_BUILTIN(X, Y, F) CASADI_MATH_FUN_BUILTIN_GEN(BinaryOperationSS, X, Y, F, 1)

    switch (op) {
    CASADI_MATH_FUN_BUILTIN(x, y, f)
      }
  }

  template<typename T>
  inline void casadi_math<T>::fun(Operation op, const T* x, const T* y, T* f, casadi_int n) {
    switch (op) {
      CASADI_MATH_FUN_BUILTIN_GEN(BinaryOperationVV, x, y, f, n)
        }
  }

  template<typename T>
  inline void casadi_math<T>::fun(Operation op, const T* x, const T& y, T* f, casadi_int n) {
    switch (op) {
      CASADI_MATH_FUN_BUILTIN_GEN(BinaryOperationVS, x, y, f, n)
        }
  }

  template<typename T>
  inline void casadi_math<T>::fun(Operation op, const T& x, const T* y, T* f, casadi_int n) {
    switch (op) {
      CASADI_MATH_FUN_BUILTIN_GEN(BinaryOperationSV, x, y, f, n)
        }
  }


  template<typename T>
  inline void casadi_math<T>::der(Operation op, const T& x, const T& y, const T& f, T* d) {
    // NOTE: We define the implementation in a preprocessor macro to be able to force inlining,
    // and to allow extensions in the VM
#define CASADI_MATH_DER_BUILTIN(X, Y, F, D)                             \
    case Operation::OP_ASSIGN:    BinaryOperation<Operation::OP_ASSIGN>::der(X, Y, F, D);     break; \
  case Operation::OP_ADD:       BinaryOperation<Operation::OP_ADD>::der(X, Y, F, D);        break; \
  case Operation::OP_SUB:       BinaryOperation<Operation::OP_SUB>::der(X, Y, F, D);        break; \
  case Operation::OP_MUL:       BinaryOperation<Operation::OP_MUL>::der(X, Y, F, D);        break; \
  case Operation::OP_DIV:       BinaryOperation<Operation::OP_DIV>::der(X, Y, F, D);        break; \
  case Operation::OP_NEG:       BinaryOperation<Operation::OP_NEG>::der(X, Y, F, D);        break; \
  case Operation::OP_EXP:       BinaryOperation<Operation::OP_EXP>::der(X, Y, F, D);        break; \
  case Operation::OP_LOG:       BinaryOperation<Operation::OP_LOG>::der(X, Y, F, D);        break; \
  case Operation::OP_POW:       BinaryOperation<Operation::OP_POW>::der(X, Y, F, D);        break; \
  case Operation::OP_CONSTPOW:  BinaryOperation<Operation::OP_CONSTPOW>::der(X, Y, F, D);   break; \
  case Operation::OP_SQRT:      BinaryOperation<Operation::OP_SQRT>::der(X, Y, F, D);       break; \
  case Operation::OP_SQ:        BinaryOperation<Operation::OP_SQ>::der(X, Y, F, D);         break; \
  case Operation::OP_TWICE:     BinaryOperation<Operation::OP_TWICE>::der(X, Y, F, D);      break; \
  case Operation::OP_SIN:       BinaryOperation<Operation::OP_SIN>::der(X, Y, F, D);        break; \
  case Operation::OP_COS:       BinaryOperation<Operation::OP_COS>::der(X, Y, F, D);        break; \
  case Operation::OP_TAN:       BinaryOperation<Operation::OP_TAN>::der(X, Y, F, D);        break; \
  case Operation::OP_ASIN:      BinaryOperation<Operation::OP_ASIN>::der(X, Y, F, D);       break; \
  case Operation::OP_ACOS:      BinaryOperation<Operation::OP_ACOS>::der(X, Y, F, D);       break; \
  case Operation::OP_ATAN:      BinaryOperation<Operation::OP_ATAN>::der(X, Y, F, D);       break; \
  case Operation::OP_LT:        BinaryOperation<Operation::OP_LT>::der(X, Y, F, D);         break; \
  case Operation::OP_LE:        BinaryOperation<Operation::OP_LE>::der(X, Y, F, D);         break; \
  case Operation::OP_EQ:        BinaryOperation<Operation::OP_EQ>::der(X, Y, F, D);         break; \
  case Operation::OP_NE:        BinaryOperation<Operation::OP_NE>::der(X, Y, F, D);         break; \
  case Operation::OP_NOT:       BinaryOperation<Operation::OP_NOT>::der(X, Y, F, D);        break; \
  case Operation::OP_AND:       BinaryOperation<Operation::OP_AND>::der(X, Y, F, D);        break; \
  case Operation::OP_OR:        BinaryOperation<Operation::OP_OR>::der(X, Y, F, D);         break; \
  case Operation::OP_IF_ELSE_ZERO: BinaryOperation<Operation::OP_IF_ELSE_ZERO>::der(X, Y, F, D);         break; \
  case Operation::OP_FLOOR:     BinaryOperation<Operation::OP_FLOOR>::der(X, Y, F, D);      break; \
  case Operation::OP_CEIL:      BinaryOperation<Operation::OP_CEIL>::der(X, Y, F, D);       break; \
  case Operation::OP_FMOD:      BinaryOperation<Operation::OP_FMOD>::der(X, Y, F, D);       break; \
  case Operation::OP_FABS:      BinaryOperation<Operation::OP_FABS>::der(X, Y, F, D);       break; \
  case Operation::OP_SIGN:      BinaryOperation<Operation::OP_SIGN>::der(X, Y, F, D);       break; \
  case Operation::OP_COPYSIGN:  BinaryOperation<Operation::OP_COPYSIGN>::der(X, Y, F, D);   break; \
  case Operation::OP_ERF:       BinaryOperation<Operation::OP_ERF>::der(X, Y, F, D);        break; \
  case Operation::OP_FMIN:      BinaryOperation<Operation::OP_FMIN>::der(X, Y, F, D);       break; \
  case Operation::OP_FMAX:      BinaryOperation<Operation::OP_FMAX>::der(X, Y, F, D);       break; \
  case Operation::OP_INV:       BinaryOperation<Operation::OP_INV>::der(X, Y, F, D);        break; \
  case Operation::OP_SINH:      BinaryOperation<Operation::OP_SINH>::der(X, Y, F, D);       break; \
  case Operation::OP_COSH:      BinaryOperation<Operation::OP_COSH>::der(X, Y, F, D);       break; \
  case Operation::OP_TANH:      BinaryOperation<Operation::OP_TANH>::der(X, Y, F, D);       break; \
  case Operation::OP_ASINH:     BinaryOperation<Operation::OP_ASINH>::der(X, Y, F, D);      break; \
  case Operation::OP_ACOSH:     BinaryOperation<Operation::OP_ACOSH>::der(X, Y, F, D);      break; \
  case Operation::OP_ATANH:     BinaryOperation<Operation::OP_ATANH>::der(X, Y, F, D);      break; \
  case Operation::OP_ATAN2:     BinaryOperation<Operation::OP_ATAN2>::der(X, Y, F, D);      break; \
  case Operation::OP_ERFINV:    BinaryOperation<Operation::OP_ERFINV>::der(X, Y, F, D);     break; \
  case Operation::OP_LIFT:      BinaryOperation<Operation::OP_LIFT>::der(X, Y, F, D);       break; \
  case Operation::OP_PRINTME:   BinaryOperation<Operation::OP_PRINTME>::der(X, Y, F, D);    break; \
  case Operation::OP_LOG1P:     BinaryOperation<Operation::OP_LOG1P>::der(X, Y, F, D);      break; \
  case Operation::OP_EXPM1:     BinaryOperation<Operation::OP_EXPM1>::der(X, Y, F, D);      break; \
  case Operation::OP_HYPOT:     BinaryOperation<Operation::OP_HYPOT>::der(X, Y, F, D);      break;
    switch (op) {
    CASADI_MATH_DER_BUILTIN(x, y, f, d)
      }
  }


    template<typename T>
      inline void casadi_math<T>::derF(Operation op, const T& x, const T& y, T& f, T* d) {
    // NOTE: We define the implementation in a preprocessor macro to be able to force inlining,
    // and to allow extensions in the VM
#define CASADI_MATH_DERF_BUILTIN(X, Y, F, D)                            \
case Operation::OP_ASSIGN:    DerBinaryOperation<Operation::OP_ASSIGN>::derf(X, Y, F, D);        break; \
case Operation::OP_ADD:       DerBinaryOperation<Operation::OP_ADD>::derf(X, Y, F, D);        break; \
case Operation::OP_SUB:       DerBinaryOperation<Operation::OP_SUB>::derf(X, Y, F, D);        break; \
case Operation::OP_MUL:       DerBinaryOperation<Operation::OP_MUL>::derf(X, Y, F, D);        break; \
case Operation::OP_DIV:       DerBinaryOperation<Operation::OP_DIV>::derf(X, Y, F, D);        break; \
case Operation::OP_NEG:       DerBinaryOperation<Operation::OP_NEG>::derf(X, Y, F, D);        break; \
case Operation::OP_EXP:       DerBinaryOperation<Operation::OP_EXP>::derf(X, Y, F, D);        break; \
case Operation::OP_LOG:       DerBinaryOperation<Operation::OP_LOG>::derf(X, Y, F, D);        break; \
case Operation::OP_POW:       DerBinaryOperation<Operation::OP_POW>::derf(X, Y, F, D);        break; \
case Operation::OP_CONSTPOW:  DerBinaryOperation<Operation::OP_CONSTPOW>::derf(X, Y, F, D);   break; \
case Operation::OP_SQRT:      DerBinaryOperation<Operation::OP_SQRT>::derf(X, Y, F, D);       break; \
case Operation::OP_SQ:        DerBinaryOperation<Operation::OP_SQ>::derf(X, Y, F, D);         break; \
case Operation::OP_TWICE:     DerBinaryOperation<Operation::OP_TWICE>::derf(X, Y, F, D);      break; \
case Operation::OP_SIN:       DerBinaryOperation<Operation::OP_SIN>::derf(X, Y, F, D);        break; \
case Operation::OP_COS:       DerBinaryOperation<Operation::OP_COS>::derf(X, Y, F, D);        break; \
case Operation::OP_TAN:       DerBinaryOperation<Operation::OP_TAN>::derf(X, Y, F, D);        break; \
case Operation::OP_ASIN:      DerBinaryOperation<Operation::OP_ASIN>::derf(X, Y, F, D);       break; \
case Operation::OP_ACOS:      DerBinaryOperation<Operation::OP_ACOS>::derf(X, Y, F, D);       break; \
case Operation::OP_ATAN:      DerBinaryOperation<Operation::OP_ATAN>::derf(X, Y, F, D);       break; \
case Operation::OP_LT:        DerBinaryOperation<Operation::OP_LT>::derf(X, Y, F, D);         break; \
case Operation::OP_LE:        DerBinaryOperation<Operation::OP_LE>::derf(X, Y, F, D);         break; \
case Operation::OP_EQ:        DerBinaryOperation<Operation::OP_EQ>::derf(X, Y, F, D);         break; \
case Operation::OP_NE:        DerBinaryOperation<Operation::OP_NE>::derf(X, Y, F, D);         break; \
case Operation::OP_NOT:       DerBinaryOperation<Operation::OP_NOT>::derf(X, Y, F, D);        break; \
case Operation::OP_AND:       DerBinaryOperation<Operation::OP_AND>::derf(X, Y, F, D);        break; \
case Operation::OP_OR:        DerBinaryOperation<Operation::OP_OR>::derf(X, Y, F, D);         break; \
case Operation::OP_IF_ELSE_ZERO: DerBinaryOperation<Operation::OP_IF_ELSE_ZERO>::derf(X, Y, F, D);         break; \
case Operation::OP_FLOOR:     DerBinaryOperation<Operation::OP_FLOOR>::derf(X, Y, F, D);      break; \
case Operation::OP_CEIL:      DerBinaryOperation<Operation::OP_CEIL>::derf(X, Y, F, D);       break; \
case Operation::OP_FMOD:      DerBinaryOperation<Operation::OP_FMOD>::derf(X, Y, F, D);       break; \
case Operation::OP_FABS:      DerBinaryOperation<Operation::OP_FABS>::derf(X, Y, F, D);       break; \
case Operation::OP_SIGN:      DerBinaryOperation<Operation::OP_SIGN>::derf(X, Y, F, D);       break; \
case Operation::OP_COPYSIGN:  DerBinaryOperation<Operation::OP_COPYSIGN>::derf(X, Y, F, D);   break; \
case Operation::OP_ERF:       DerBinaryOperation<Operation::OP_ERF>::derf(X, Y, F, D);        break; \
case Operation::OP_FMIN:      DerBinaryOperation<Operation::OP_FMIN>::derf(X, Y, F, D);       break; \
case Operation::OP_FMAX:      DerBinaryOperation<Operation::OP_FMAX>::derf(X, Y, F, D);       break; \
case Operation::OP_INV:       DerBinaryOperation<Operation::OP_INV>::derf(X, Y, F, D);        break; \
case Operation::OP_SINH:      DerBinaryOperation<Operation::OP_SINH>::derf(X, Y, F, D);       break; \
case Operation::OP_COSH:      DerBinaryOperation<Operation::OP_COSH>::derf(X, Y, F, D);       break; \
case Operation::OP_TANH:      DerBinaryOperation<Operation::OP_TANH>::derf(X, Y, F, D);       break; \
case Operation::OP_ASINH:     DerBinaryOperation<Operation::OP_ASINH>::derf(X, Y, F, D);      break; \
case Operation::OP_ACOSH:     DerBinaryOperation<Operation::OP_ACOSH>::derf(X, Y, F, D);      break; \
case Operation::OP_ATANH:     DerBinaryOperation<Operation::OP_ATANH>::derf(X, Y, F, D);      break; \
case Operation::OP_ATAN2:     DerBinaryOperation<Operation::OP_ATAN2>::derf(X, Y, F, D);      break; \
case Operation::OP_ERFINV:    DerBinaryOperation<Operation::OP_ERFINV>::derf(X, Y, F, D);     break; \
case Operation::OP_LIFT:      DerBinaryOperation<Operation::OP_LIFT>::derf(X, Y, F, D);       break; \
case Operation::OP_PRINTME:   DerBinaryOperation<Operation::OP_PRINTME>::derf(X, Y, F, D);    break; \
case Operation::OP_LOG1P:     DerBinaryOperation<Operation::OP_LOG1P>::derf(X, Y, F, D);      break; \
case Operation::OP_EXPM1:     DerBinaryOperation<Operation::OP_EXPM1>::derf(X, Y, F, D);      break; \
case Operation::OP_HYPOT:     DerBinaryOperation<Operation::OP_HYPOT>::derf(X, Y, F, D);      break;
    switch (op) {
      CASADI_MATH_DERF_BUILTIN(x, y, f, d)
        }
  }

  #define CASADI_MATH_BINARY_BUILTIN              \
    case Operation::OP_ADD:                                  \
    case Operation::OP_SUB:                                  \
    case Operation::OP_MUL:                                  \
    case Operation::OP_DIV:                                  \
    case Operation::OP_POW:                                  \
    case Operation::OP_CONSTPOW:                             \
    case Operation::OP_LT:                                   \
    case Operation::OP_LE:                                   \
    case Operation::OP_EQ:                                   \
    case Operation::OP_NE:                                   \
    case Operation::OP_AND:                                  \
    case Operation::OP_OR:                                   \
    case Operation::OP_COPYSIGN:                             \
    case Operation::OP_FMOD:                                 \
    case Operation::OP_FMIN:                                 \
    case Operation::OP_FMAX:                                 \
    case Operation::OP_ATAN2:                                \
    case Operation::OP_PRINTME:                              \
    case Operation::OP_LIFT:                                 \
    case Operation::OP_HYPOT:

  #define CASADI_MATH_UNARY_BUILTIN              \
    case Operation::OP_ASSIGN:                              \
    case Operation::OP_NEG:                                 \
    case Operation::OP_EXP:                                 \
    case Operation::OP_LOG:                                 \
    case Operation::OP_SQRT:                                \
    case Operation::OP_SQ:                                  \
    case Operation::OP_TWICE:                               \
    case Operation::OP_SIN:                                 \
    case Operation::OP_COS:                                 \
    case Operation::OP_TAN:                                 \
    case Operation::OP_ASIN:                                \
    case Operation::OP_ACOS:                                \
    case Operation::OP_ATAN:                                \
    case Operation::OP_FLOOR:                               \
    case Operation::OP_CEIL:                                \
    case Operation::OP_NOT:                                 \
    case Operation::OP_ERF:                                 \
    case Operation::OP_FABS:                                \
    case Operation::OP_SIGN:                                \
    case Operation::OP_INV:                                 \
    case Operation::OP_SINH:                                \
    case Operation::OP_COSH:                                \
    case Operation::OP_TANH:                                \
    case Operation::OP_ASINH:                               \
    case Operation::OP_ACOSH:                               \
    case Operation::OP_ATANH:                               \
    case Operation::OP_ERFINV:                              \
    case Operation::OP_LOG1P:                               \
    case Operation::OP_EXPM1:

  template<typename T>
  bool casadi_math<T>::is_binary(Operation op) {
    switch (op) {
      CASADI_MATH_BINARY_BUILTIN
      case Operation::OP_IF_ELSE_ZERO:
        return true;
      default:
        return false;
    }
  }

  template<typename T>
  bool casadi_math<T>::is_unary(Operation op) {
    switch (op) {
      CASADI_MATH_UNARY_BUILTIN
      return true;
    default:
      return false;
    }
  }

  template<typename T>
  inline casadi_int casadi_math<T>::ndeps(Operation op) {
    switch (op) {
      case Operation::OP_CONST:
      case Operation::OP_PARAMETER:
      case Operation::OP_INPUT:
        return 0;
      CASADI_MATH_BINARY_BUILTIN
      case Operation::OP_IF_ELSE_ZERO:
        return 2;
      default:
        return 1;
    }
  }

  template<typename T>
  inline std::string
  casadi_math<T>::print(Operation op,
                        const std::string& x, const std::string& y) {
    casadi_assert_dev(ndeps(op)==2);
    return pre(op) + x + sep(op) + y + post(op);
  }

  template<typename T>
  inline std::string
  casadi_math<T>::print(Operation op, const std::string& x) {
    casadi_assert_dev(ndeps(op)==1);
    return pre(op) + x + post(op);
  }

  static inline const char* get_operation_name_c(Operation op) {
    switch (op) {
    case Operation::OP_ASSIGN:         return "assign";
    case Operation::OP_ADD:            return "add";
    case Operation::OP_SUB:            return "sub";
    case Operation::OP_MUL:            return "mul";
    case Operation::OP_DIV:            return "div";
    case Operation::OP_NEG:            return "neg";
    case Operation::OP_EXP:            return "exp";
    case Operation::OP_LOG:            return "log";
    case Operation::OP_CONSTPOW:
    case Operation::OP_POW:            return "pow";
    case Operation::OP_SQRT:           return "sqrt";
    case Operation::OP_SQ:             return "sq";
    case Operation::OP_TWICE:          return "twice";
    case Operation::OP_SIN:            return "sin";
    case Operation::OP_COS:            return "cos";
    case Operation::OP_TAN:            return "tan";
    case Operation::OP_ASIN:           return "asin";
    case Operation::OP_ACOS:           return "acos";
    case Operation::OP_ATAN:           return "atan";
    case Operation::OP_LT:             return "lt";
    case Operation::OP_LE:             return "le";
    case Operation::OP_EQ:             return "eq";
    case Operation::OP_NE:             return "ne";
    case Operation::OP_NOT:            return "not";
    case Operation::OP_AND:            return "and";
    case Operation::OP_OR:             return "or";
    case Operation::OP_FLOOR:          return "floor";
    case Operation::OP_CEIL:           return "ceil";
    case Operation::OP_FMOD:           return "fmod";
    case Operation::OP_FABS:           return "fabs";
    case Operation::OP_SIGN:           return "sign";
    case Operation::OP_COPYSIGN:       return "copysign";
    case Operation::OP_IF_ELSE_ZERO:   return "if_else_zero";
    case Operation::OP_ERF:            return "erf";
    case Operation::OP_FMIN:           return "fmin";
    case Operation::OP_FMAX:           return "fmax";
    case Operation::OP_INV:            return "inv";
    case Operation::OP_SINH:           return "sinh";
    case Operation::OP_COSH:           return "cosh";
    case Operation::OP_TANH:           return "tanh";
    case Operation::OP_ASINH:          return "asinh";
    case Operation::OP_ACOSH:          return "acosh";
    case Operation::OP_ATANH:          return "atanh";
    case Operation::OP_ATAN2:          return "atan2";
    case Operation::OP_CONST:          return "const";
    case Operation::OP_INPUT:          return "input";
    case Operation::OP_OUTPUT:         return "output";
    case Operation::OP_PARAMETER:      return "parameter";
    case Operation::OP_CALL:           return "call";
    case Operation::OP_MTIMES:         return "mtimes";
    case Operation::OP_SOLVE:          return "solve";
    case Operation::OP_TRANSPOSE:      return "transpose";
    case Operation::OP_DETERMINANT:    return "determinant";
    case Operation::OP_INVERSE:        return "inverse";
    case Operation::OP_DOT:            return "dot";
    case Operation::OP_HORZCAT:        return "horzcat";
    case Operation::OP_VERTCAT:        return "vertcat";
    case Operation::OP_DIAGCAT:        return "diagcat";
    case Operation::OP_HORZSPLIT:      return "horzsplit";
    case Operation::OP_VERTSPLIT:      return "vertsplit";
    case Operation::OP_DIAGSPLIT:      return "diagsplit";
    case Operation::OP_RESHAPE:        return "reshape";
    case Operation::OP_SUBREF:         return "subref";
    case Operation::OP_SUBASSIGN:      return "subassign";
    case Operation::OP_GETNONZEROS:    return "getnonzeros";
    case Operation::OP_GETNONZEROS_PARAM:    return "getnonzeros_param";
    case Operation::OP_ADDNONZEROS:    return "addnonzeros";
    case Operation::OP_ADDNONZEROS_PARAM:    return "addnonzeros_param";
    case Operation::OP_SETNONZEROS:    return "setnonzeros";
    case Operation::OP_SETNONZEROS_PARAM:    return "setnonzeros_param";
    case Operation::OP_PROJECT:        return "project";
    case Operation::OP_ASSERTION:      return "assertion";
    case Operation::OP_NORM2:          return "norm2";
    case Operation::OP_NORM1:          return "norm1";
    case Operation::OP_NORMINF:        return "norminf";
    case Operation::OP_NORMF:          return "normf";
    case Operation::OP_ERFINV:         return "erfinv";
    case Operation::OP_PRINTME:        return "printme";
    case Operation::OP_LIFT:           return "lift";
    case Operation::OP_EINSTEIN:       return "einstein";
    case Operation::OP_BSPLINE:        return "bspline";
    case Operation::OP_CONVEXIFY:      return "convexify";
    case Operation::OP_LOG1P:          return "log1p";
    case Operation::OP_EXPM1:          return "expm1";
    case Operation::OP_HYPOT:          return "hypot";
    case Operation::OP_LOGSUMEXP:      return "logsumexp";
    }
    return nullptr;
  }

  static inline std::string get_operation_name(Operation op) {
    return std::move(std::string(get_operation_name_c(op)));
  }

  template<typename T>
  inline std::string casadi_math<T>::name(Operation op) {
    return std::move(get_operation_name(op));
  }

  template<typename T>
  inline std::string casadi_math<T>::pre(Operation op) {
    switch (op) {
    case Operation::OP_ASSIGN:    return "";
    case Operation::OP_ADD:       return "(";
    case Operation::OP_SUB:       return "(";
    case Operation::OP_MUL:       return "(";
    case Operation::OP_DIV:       return "(";
    case Operation::OP_NEG:       return "(-";
    case Operation::OP_TWICE:     return "(2.*";
    case Operation::OP_LT:        return "(";
    case Operation::OP_LE:        return "(";
    case Operation::OP_EQ:        return "(";
    case Operation::OP_NE:        return "(";
    case Operation::OP_NOT:       return "(!";
    case Operation::OP_AND:       return "(";
    case Operation::OP_OR:        return "(";
    case Operation::OP_IF_ELSE_ZERO: return "(";
    case Operation::OP_INV:       return "(1./";
    default: return name(op) + "(";
    }
  }

  template<typename T>
  inline std::string casadi_math<T>::sep(Operation op) {
    switch (op) {
    case Operation::OP_ADD:       return "+";
    case Operation::OP_SUB:       return "-";
    case Operation::OP_MUL:       return "*";
    case Operation::OP_DIV:       return "/";
    case Operation::OP_LT:        return "<";
    case Operation::OP_LE:        return "<=";
    case Operation::OP_EQ:        return "==";
    case Operation::OP_NE:        return "!=";
    case Operation::OP_AND:       return "&&";
    case Operation::OP_OR:        return "||";
    case Operation::OP_IF_ELSE_ZERO: return "?";
    default:           return ",";
    }
  }

  template<typename T>
  inline std::string casadi_math<T>::post(Operation op) {
    switch (op) {
    case Operation::OP_ASSIGN:       return "";
    case Operation::OP_IF_ELSE_ZERO: return ":0)";
    default:              return ")";
    }
  }

#endif // SWIG

} // namespace casadi

/// \endcond

#endif // CASADI_CALCULUS_HPP
