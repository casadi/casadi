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
#include "exception.hpp"

// Define pi if the compiler fails to do so
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI

/// \cond INTERNAL

namespace casadi {

  /// Enum for quick access to any node
  enum Operation {
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

    // Embedded function call in parallel
    OP_MAP,

    // Matrix multiplication
    OP_MATMUL,

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

    // Nonzero addition
    OP_ADDNONZEROS,

    // Nonzero assignment
    OP_SETNONZEROS,

    // Set sparse
    OP_PROJECT,

    // Assertion
    OP_ASSERTION,

    // Monitor
    OP_MONITOR,

    // Norms
    OP_NORM2, OP_NORM1, OP_NORMINF, OP_NORMF,

    // Horizontal repeat
    OP_HORZREPMAT,

    // Horizontal repeat sum
    OP_HORZREPSUM,

    OP_ERFINV,
    OP_PRINTME,
    OP_LIFT,
    NUM_BUILT_IN_OPS
  };

#ifndef SWIG

  // Get GCC version if GCC is used
#ifdef __GNUC__
#ifdef __GNUC_MINOR__
#ifdef __GNUC_PATCHLEVEL__
#define GCC_VERSION (__GNUC__ * 10000 +__GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#endif // __GNUC_PATCHLEVEL__
#endif // __GNUC_MINOR__
#endif // __GNUC__

  // Disable some Visual studio warnings
#ifdef _MSC_VER

#pragma warning(disable:4996)

  // warning C4018: '<' : signed/unsigned mismatch
#pragma warning(disable:4018)

  // warning C4800: 'int' : forcing value to bool 'true'or 'false'(performance warning)
#pragma warning(disable:4800)
#endif

  ///@{
  /** \brief Enable using elementary numerical operations without std:: prefix */
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
  ///@}

  ///@{
  // Implement "missing" operations
  inline double atanh(double x) throw() {
    if (x==-1) return -std::numeric_limits<double>::infinity();
    if (x==1) return std::numeric_limits<double>::infinity();
    return 0.5*log((1+x)/(1-x));
  }

  inline double asinh(double x) throw() {
    return log(x + sqrt(1+x*x));
  }

  inline double acosh(double x) throw() {
    return log(x + sqrt(1+x)*sqrt(x-1));
  }

  inline int isnan(double x) throw() { return x!=x;}
  inline int isinf(double x) throw() { return isnan(x-x);}

  /// Sign function, note that sign(nan) == nan
  inline double sign(double x) { return x<0 ? -1 : x>0 ? 1 : x;}

  /// fmin, fmax and erf should be available if C99 and/or C++11 required
  inline double fmin(double x, double y) throw() { return std::min(x, y);}
  inline int fmin(int x, int y) throw() { return std::min(x, y);}
  inline double fmax(double x, double y) throw() { return std::max(x, y);}
  inline int fmax(int x, int y) throw() { return std::max(x, y);}

  /// fabs(int) was added in C++11
  inline int fabs(int x) throw() { return std::abs(x);}
  ///@}

#ifdef HAS_ERF
  using ::erf;
#else // HAS ERF
  inline double erf(double x) throw() {
    // Approximation found in Sourceforge and modified,
    // originally from numerical recipes in Fortran
    double sx = x<0 ? -1 : x>0 ? 1 : x;
    double z = sx*x;
    double t = 1.0/(1.0+0.5*z);
    return 1.-sx*(t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
                                       t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
                                                           t*(-0.82215223+t*0.17087277))))))))));
  }
#endif // HAS ERF
  ///@}

  ///@{
  /** \brief  CasADi additions */
  template<class T> T constpow(const T &x, const T &n) { return x.constpow(n);}
  template<class T> T printme(const T &x, const T &y) { return x.printme(y);}
  inline double printme(double x, double y) {
    std::cout << "|> " << y << " : " << x << std::endl;
    return x;
  }

  #ifdef HAS_COPYSIGN
  using std::copysign;
  #else
  /// copysign function
  inline double copysign(double x, double y) { return y>=0 ? fabs(x) : -fabs(x);}
  #endif //HAS_COPYSIGN

  /// Conditional assignment
  inline double if_else_zero(double x, double y) { return x ? y : 0;}

  /// Inverse of the error function
  template<class T> T erfinv(const T &x) {return x.zz_erfinv();}
#ifdef HAS_ERFINV
  using ::erfinv;
#else // HAS ERFINV
  inline double erfinv(double x) throw() {
    // Approximation found in Sourceforge and modified: Not very efficient
    if (x>=1) {
      return x==1 ? std::numeric_limits<double>::infinity() :
          std::numeric_limits<double>::quiet_NaN();
    } else if (x<=-1) {
      return x==-1 ? -std::numeric_limits<double>::infinity() :
          std::numeric_limits<double>::quiet_NaN();
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
      y = y - (erf(y) - x) / (2.0/sqrt(M_PI) * exp(-y*y));
      y = y - (erf(y) - x) / (2.0/sqrt(M_PI) * exp(-y*y));
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

  template<int I>
  struct UnaryOperation {
    /// Function evaluation
    template<typename T> static inline void fcn(const T& x, T& f);

    /// Partial derivatives
    template<typename T> static inline void der(const T& x, const T& f, T* d);
  };

  template<int I>
  struct BinaryOperation {
    /// Function evaluation
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) {
        UnaryOperation<I>::fcn(x, f);}

    /// Partial derivatives - binary function
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        UnaryOperation<I>::der(x, f, d); d[1]=0; }
  };

  template<int I>
  struct BinaryOperationE {
    /// Function evaluation
    template<typename T> static inline T fcn(const T& x, const T& y) {
      T ret;
      BinaryOperation<I>::fcn(x, y, ret);
      return ret;
    }
  };

  /// Calculate function and derivative
  template<int I>
  struct DerBinaryOpertion {
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
  template<int I>
  struct BinaryOperationSS {
    /// Function evaluation
    template<typename T> static inline void fcn(const T& x, const T& y, T& f, int n) {
      BinaryOperation<I>::fcn(x, y, f);
    }

    /// Partial derivatives - binary function
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d, int n) {
      BinaryOperation<I>::der(x, y, f, d);
    }
  };


  /// Perform a binary operation on two vectors
  template<int I>
  struct BinaryOperationVV {
    /// Function evaluation
    template<typename T> static inline void fcn(const T* x, const T* y, T* f, int n) {
      for (int i=0; i<n; ++i) {
        BinaryOperation<I>::fcn(*x++, *y++, *f++);
      }
    }

    /// Partial derivatives - binary function
    template<typename T> static inline void der(const T* x, const T* y, const T* f, T* d, int n) {
      for (int i=0; i<n; ++i, d+=2) {
        BinaryOperation<I>::der(*x++, *y++, *f++, d);
      }
    }
  };

  /// Perform a binary operation on a vector and a scalar
  template<int I>
  struct BinaryOperationVS {
    /// Function evaluation
    template<typename T> static inline void fcn(const T* x, const T& y, T* f, int n) {
      for (int i=0; i<n; ++i) {
        BinaryOperation<I>::fcn(*x++, y, *f++);
      }
    }

    /// Partial derivatives - binary function
    template<typename T> static inline void der(const T* x, const T& y, const T* f, T* d, int n) {
      for (int i=0; i<n; ++i, d+=2) {
        BinaryOperation<I>::der(*x++, y, *f++, d);
      }
    }
  };

  /// Perform a binary operation on a scalar and a vector
  template<int I>
  struct BinaryOperationSV {
    /// Function evaluation
    template<typename T> static inline void fcn(const T& x, const T* y, T* f, int n) {
      for (int i=0; i<n; ++i) {
        BinaryOperation<I>::fcn(x, *y++, *f++);
      }
    }

    /// Partial derivatives - binary function
    template<typename T> static inline void der(const T& x, const T* y, const T* f, T* d, int n) {
      for (int i=0; i<n; ++i, d+=2) {
        BinaryOperation<I>::der(x, *y++, *f++, d);
      }
    }
  };

  ///@{
  /// Smoothness (by default true)
  template<int I> struct SmoothChecker { static const bool check=true;};
  template<>      struct SmoothChecker<OP_LT>{ static const bool check=false;};
  template<>      struct SmoothChecker<OP_LE>{ static const bool check=false;};
  template<>      struct SmoothChecker<OP_FLOOR>{ static const bool check=false;};
  template<>      struct SmoothChecker<OP_CEIL>{ static const bool check=false;};
  template<>      struct SmoothChecker<OP_FMOD>{ static const bool check=false;};
  template<>      struct SmoothChecker<OP_EQ>{ static const bool check=false;};
  template<>      struct SmoothChecker<OP_NE>{ static const bool check=false;};
  template<>      struct SmoothChecker<OP_SIGN>{ static const bool check=false;};
  template<>      struct SmoothChecker<OP_COPYSIGN>{ static const bool check=false;};
  template<>      struct SmoothChecker<OP_NOT>{ static const bool check=false;};
  template<>      struct SmoothChecker<OP_AND>{ static const bool check=false;};
  template<>      struct SmoothChecker<OP_OR>{ static const bool check=false;};
  template<>      struct SmoothChecker<OP_IF_ELSE_ZERO>{ static const bool check=false;};
  ///@}

  ///@{
  /// If evaluated with the first argument zero, is the result zero?
  template<int I> struct F0XChecker { static const bool check=false;};
  template<>      struct F0XChecker<OP_ASSIGN>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_MUL>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_DIV>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_NEG>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_POW>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_CONSTPOW>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_SQRT>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_SQ>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_TWICE>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_SIN>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_TAN>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_ATAN>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_ASIN>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_FLOOR>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_CEIL>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_FMOD>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_FABS>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_SIGN>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_COPYSIGN>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_ERF>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_SINH>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_TANH>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_ASINH>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_ATANH>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_ERFINV>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_AND>{ static const bool check=true;};
  template<>      struct F0XChecker<OP_IF_ELSE_ZERO>{ static const bool check=true;};
  ///@}

  ///@{
  /// If evaluated with the second argument zero, is the result zero?
  template<int I> struct Function0Checker { static const bool check=false;};
  template<>      struct Function0Checker<OP_MUL>{ static const bool check=true;};
  template<>      struct Function0Checker<OP_AND>{ static const bool check=true;};
  template<>      struct Function0Checker<OP_IF_ELSE_ZERO>{ static const bool check=true;};
  ///@}

  ///@{
  /// If evaluated with both arguments zero, is the result zero?
  template<int I> struct F00Checker { static const bool check=F0XChecker<I>::check;};
  template<>      struct F00Checker<OP_ADD>{ static const bool check=true;};
  template<>      struct F00Checker<OP_SUB>{ static const bool check=true;};
  template<>      struct F00Checker<OP_FMIN>{ static const bool check=true;};
  template<>      struct F00Checker<OP_FMAX>{ static const bool check=true;};
  template<>      struct F00Checker<OP_AND>{ static const bool check=true;};
  template<>      struct F00Checker<OP_OR>{ static const bool check=true;};
  template<>      struct F00Checker<OP_COPYSIGN>{ static const bool check=true;};
  ///@}

  ///@{
  /// Is commutative
  template<int I> struct CommChecker { static const bool check=false;};
  template<>      struct CommChecker<OP_ADD>{ static const bool check=true;};
  template<>      struct CommChecker<OP_MUL>{ static const bool check=true;};
  template<>      struct CommChecker<OP_EQ>{ static const bool check=true;};
  template<>      struct CommChecker<OP_NE>{ static const bool check=true;};
  template<>      struct CommChecker<OP_AND>{ static const bool check=true;};
  template<>      struct CommChecker<OP_OR>{ static const bool check=true;};
  ///@}

  ///@{
  /// Always non-negative (false by default)
  template<int I> struct NonnegativeChecker { static const bool check=false;};
  template<>      struct NonnegativeChecker<OP_SQRT>{ static const bool check=true;};
  template<>      struct NonnegativeChecker<OP_SQ>{ static const bool check=true;};
  template<>      struct NonnegativeChecker<OP_EXP>{ static const bool check=true;};
  template<>      struct NonnegativeChecker<OP_LT>{ static const bool check=true;};
  template<>      struct NonnegativeChecker<OP_LE>{ static const bool check=true;};
  template<>      struct NonnegativeChecker<OP_EQ>{ static const bool check=true;};
  template<>      struct NonnegativeChecker<OP_NE>{ static const bool check=true;};
  template<>      struct NonnegativeChecker<OP_NOT>{ static const bool check=true;};
  template<>      struct NonnegativeChecker<OP_AND>{ static const bool check=true;};
  template<>      struct NonnegativeChecker<OP_OR>{ static const bool check=true;};
  ///@}

  ///@{
  /// Is the operation binary as opposed to unary
  template<int I> struct BinaryChecker { static const bool check=false;};
  template<>      struct BinaryChecker<OP_ADD>{ static const bool check=true;};
  template<>      struct BinaryChecker<OP_SUB>{ static const bool check=true;};
  template<>      struct BinaryChecker<OP_MUL>{ static const bool check=true;};
  template<>      struct BinaryChecker<OP_DIV>{ static const bool check=true;};
  template<>      struct BinaryChecker<OP_POW>{ static const bool check=true;};
  template<>      struct BinaryChecker<OP_CONSTPOW>{ static const bool check=true;};
  template<>      struct BinaryChecker<OP_EQ>{ static const bool check=true;};
  template<>      struct BinaryChecker<OP_NE>{ static const bool check=true;};
  template<>      struct BinaryChecker<OP_AND>{ static const bool check=true;};
  template<>      struct BinaryChecker<OP_OR>{ static const bool check=true;};
  template<>      struct BinaryChecker<OP_FMIN>{ static const bool check=true;};
  template<>      struct BinaryChecker<OP_FMAX>{ static const bool check=true;};
  template<>      struct BinaryChecker<OP_PRINTME>{ static const bool check=true;};
  template<>      struct BinaryChecker<OP_ATAN2>{ static const bool check=true;};
  template<>      struct BinaryChecker<OP_IF_ELSE_ZERO>{ static const bool check=true;};
  template<>      struct BinaryChecker<OP_COPYSIGN>{ static const bool check=true;};
  ///@}

  /// Simple assignment
  template<>
  struct UnaryOperation<OP_ASSIGN>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = x;}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0] = 1; }
  };

  /// Addition
  template<>
  struct BinaryOperation<OP_ADD>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = x+y;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=d[1]=1;}
  };

  /// Subtraction
  template<>
  struct BinaryOperation<OP_SUB>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = x-y;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=1; d[1]=-1;}
  };

  /// Multiplication
  template<>
  struct BinaryOperation<OP_MUL>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = x*y;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=y; d[1]=x;}
  };

  /// Division
  template<>
  struct BinaryOperation<OP_DIV>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = x/y;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=1/y; d[1]=-f/y;}
  };

  /// Negation
  template<>
  struct UnaryOperation<OP_NEG>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = -x;}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0]=-1;}
  };

  /// Natural exponent
  template<>
  struct UnaryOperation<OP_EXP>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = exp(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0]=f;}
  };

  /// Natural logarithm
  template<>
  struct UnaryOperation<OP_LOG>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = log(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0]=1/x;}
  };

  /// Power, defined only for x>=0
  template<>
  struct BinaryOperation<OP_POW>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = pow(x, y);}
    // See issue #104 why d[0] is no longer y*f/x
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=y*pow(x, y-1); d[1]=log(x)*f;}
  };

  /// Power, defined only for y constant
  template<>
  struct BinaryOperation<OP_CONSTPOW>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = pow(x, y);}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=y*pow(x, y-1); d[1]=0;}
  };

  /// Square root
  template<>
  struct UnaryOperation<OP_SQRT>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = sqrt(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0]=1/(twice(f));}
  };

  /// Square
  template<>
  struct UnaryOperation<OP_SQ>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = sq(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0]=twice(x);}
  };

  /// Times two
  template<>
  struct UnaryOperation<OP_TWICE>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = 2.*x;}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0] = 2; }
  };

  /// Sine
  template<>
  struct UnaryOperation<OP_SIN>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = sin(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0]=cos(x);}
  };

  /// Cosine
  template<>
  struct UnaryOperation<OP_COS>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = cos(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0]=-sin(x);}
  };

  /// Tangent
  template<>
  struct UnaryOperation<OP_TAN>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = tan(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d)
    { d[0] = 1/sq(cos(x));}
  };

  /// Arcus sine
  template<>
  struct UnaryOperation<OP_ASIN>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = asin(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0]=1/sqrt(1-x*x);}
  };

  /// Arcus cosine
  template<>
  struct UnaryOperation<OP_ACOS>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = acos(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d)
    { d[0]=-1/sqrt(1-x*x);}
  };

  /// Arcus tangent
  template<>
  struct UnaryOperation<OP_ATAN>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = atan(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0] = 1/(1+x*x);}
  };

  /// Less than
  template<>
  struct BinaryOperation<OP_LT>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = x < y;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=d[1]=0;}
  };

  /// Less or equal to
  template<>
  struct BinaryOperation<OP_LE>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = x <= y;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=d[1]=0;}
  };

  /// Floor function
  template<>
  struct UnaryOperation<OP_FLOOR>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = floor(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0] = 0;}
  };

  /// Ceil function
  template<>
  struct UnaryOperation<OP_CEIL>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = ceil(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0] = 0;}
  };

  /// Remainder of division
  template<>
  struct BinaryOperation<OP_FMOD>{
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = fmod(x, y);}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
      d[0]=1; d[1]=(f-x)/y;}
  };

  /// Equal to
  template<>
  struct BinaryOperation<OP_EQ>{
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = x==y;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=d[1]=0;}
  };

  /// Not equal to
  template<>
  struct BinaryOperation<OP_NE>{
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = x!=y;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=d[1]=0;}
  };

  /// Logical not
  template<>
  struct UnaryOperation<OP_NOT>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f) { f = !x;}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0] = 0;}
  };

  /// Logical and
  template<>
  struct BinaryOperation<OP_AND>{
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = x && y;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=d[1]=0;}
  };

  /// Logical or
  template<>
  struct BinaryOperation<OP_OR>{
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = x || y;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=d[1]=0;}
  };

  /// Error function
  template<>
  struct UnaryOperation<OP_ERF>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = erf(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) {
        d[0] = (2/sqrt(M_PI))*exp(-x*x);}
  };

  /// Absolute value
  template<>
  struct UnaryOperation<OP_FABS>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = fabs(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) {
        d[0]=sign(x);}
  };

  /// Sign
  template<>
  struct UnaryOperation<OP_SIGN>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = sign(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0]=0;}
  };

  /// Copysign
  template<>
  struct BinaryOperation<OP_COPYSIGN>{
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = copysign(x, y);}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        T e = 1; d[0]=copysign(e, y); d[1]=0;}
  };

  /// Minimum
  template<>
  struct BinaryOperation<OP_FMIN>{
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = fmin(x, y);}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=x<=y; d[1]=!d[0];}
  };

  /// Maximum
  template<>
  struct BinaryOperation<OP_FMAX>{
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = fmax(x, y);}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=x>=y; d[1]=!d[0];}
  };

  /// Elementwise inverse
  template<>
  struct UnaryOperation<OP_INV>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = 1./x;}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0] = -f*f; }
  };

  /// Hyperbolic sine
  template<>
  struct UnaryOperation<OP_SINH>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = sinh(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0] = cosh(x); }
  };

  /// Hyperbolic cosine
  template<>
  struct UnaryOperation<OP_COSH>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = cosh(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0] = -sinh(x); }
  };

  /// Hyperbolic tangent
  template<>
  struct UnaryOperation<OP_TANH>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = tanh(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0] = 1-f*f; }
  };

  /// Inverse hyperbolic sine
  template<>
  struct UnaryOperation<OP_ASINH>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = asinh(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) {
        d[0] = 1/sqrt(1+x*x); }
  };

  /// Inverse hyperbolic cosine
  template<>
  struct UnaryOperation<OP_ACOSH>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = acosh(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) {
        d[0] = 1/sqrt(x-1)/sqrt(x+1); }
  };

  /// Inverse hyperbolic tangent
  template<>
  struct UnaryOperation<OP_ATANH>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = atanh(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) { d[0] = 1/(1-x*x); }
  };

  /// Inverse of error function
  template<>
  struct UnaryOperation<OP_ERFINV>{
    template<typename T> static inline void fcn(const T& x, T& f) { f = erfinv(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d) {
        d[0] = (sqrt(M_PI)/2)*exp(f*f); }
  };

  /// Identity operator with the side effect of printing
  template<>
  struct BinaryOperation<OP_PRINTME>{
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) {f = printme(x, y); }
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=1; d[1]=0;}
  };

  /// Arctan2
  template<>
  struct BinaryOperation<OP_ATAN2>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = atan2(x, y);}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        T t = x*x+y*y; d[0]=y/t; d[1]=-x/t;}
  };

  /// Conditional assignment
  template<>
  struct BinaryOperation<OP_IF_ELSE_ZERO>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) {
        f = if_else_zero(x, y);}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0]=0; d[1]=x;}
  };

  /// Inverse of error function
  template<>
  struct BinaryOperation<OP_LIFT>{
    template<typename T> static inline void fcn(const T& x, const T& y, T& f) { f = x;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d) {
        d[0] = 1; d[1] = 0; }
  };

  template<template<int> class F>
  bool operation_checker(unsigned int op) {
    switch (op) {
    case OP_ASSIGN:        return F<OP_ASSIGN>::check;
    case OP_ADD:           return F<OP_ADD>::check;
    case OP_SUB:           return F<OP_SUB>::check;
    case OP_MUL:           return F<OP_MUL>::check;
    case OP_DIV:           return F<OP_DIV>::check;
    case OP_NEG:           return F<OP_NEG>::check;
    case OP_EXP:           return F<OP_EXP>::check;
    case OP_LOG:           return F<OP_LOG>::check;
    case OP_POW:           return F<OP_POW>::check;
    case OP_CONSTPOW:      return F<OP_CONSTPOW>::check;
    case OP_SQRT:          return F<OP_SQRT>::check;
    case OP_SQ:            return F<OP_SQ>::check;
    case OP_TWICE:         return F<OP_TWICE>::check;
    case OP_SIN:           return F<OP_SIN>::check;
    case OP_COS:           return F<OP_COS>::check;
    case OP_TAN:           return F<OP_TAN>::check;
    case OP_ASIN:          return F<OP_ASIN>::check;
    case OP_ACOS:          return F<OP_ACOS>::check;
    case OP_ATAN:          return F<OP_ATAN>::check;
    case OP_LT:            return F<OP_LT>::check;
    case OP_LE:            return F<OP_LE>::check;
    case OP_EQ:            return F<OP_EQ>::check;
    case OP_NE:            return F<OP_NE>::check;
    case OP_NOT:           return F<OP_NOT>::check;
    case OP_AND:           return F<OP_AND>::check;
    case OP_OR:            return F<OP_OR>::check;
    case OP_IF_ELSE_ZERO:   return F<OP_IF_ELSE_ZERO>::check;
    case OP_FLOOR:         return F<OP_FLOOR>::check;
    case OP_CEIL:          return F<OP_CEIL>::check;
    case OP_FMOD:          return F<OP_FMOD>::check;
    case OP_FABS:          return F<OP_FABS>::check;
    case OP_SIGN:          return F<OP_SIGN>::check;
    case OP_COPYSIGN:      return F<OP_COPYSIGN>::check;
    case OP_ERF:           return F<OP_ERF>::check;
    case OP_FMIN:          return F<OP_FMIN>::check;
    case OP_FMAX:          return F<OP_FMAX>::check;
    case OP_INV:           return F<OP_INV>::check;
    case OP_SINH:          return F<OP_SINH>::check;
    case OP_COSH:          return F<OP_COSH>::check;
    case OP_TANH:          return F<OP_TANH>::check;
    case OP_ASINH:         return F<OP_ASINH>::check;
    case OP_ACOSH:         return F<OP_ACOSH>::check;
    case OP_ATANH:         return F<OP_ATANH>::check;
    case OP_CONST:         return F<OP_CONST>::check;
    case OP_CALL:          return F<OP_CALL>::check;
    case OP_INPUT:         return F<OP_INPUT>::check;
    case OP_OUTPUT:        return F<OP_OUTPUT>::check;
    case OP_ERFINV:        return F<OP_ERFINV>::check;
    case OP_PRINTME:       return F<OP_PRINTME>::check;
    case OP_ATAN2:         return F<OP_ATAN2>::check;
    case OP_LIFT:          return F<OP_LIFT>::check;
    }

    // False by default
    return false;
  }


  /// Easy access to all the functions for a particular type
  template<typename T>
  struct casadi_math {

    /** \brief Evaluate a built in function (scalar-scalar) */
    static inline void fun(unsigned char op, const T& x, const T& y, T& f);

    /** \brief Evaluate a built in function (vector-vector) */
    static inline void fun(unsigned char op, const T* x, const T* y, T* f, int n);

    /** \brief Evaluate a built in function (vector-scalar) */
    static inline void fun(unsigned char op, const T* x, const T& y, T* f, int n);

    /** \brief Evaluate a built in function (scalar-vector) */
    static inline void fun(unsigned char op, const T& x, const T* y, T* f, int n);

    /** \brief Evaluate a built in derivative function */
    static inline void der(unsigned char op, const T& x, const T& y, const T& f, T* d);

    /** \brief Evaluate the function and the derivative function */
    static inline void derF(unsigned char op, const T& x, const T& y, T& f, T* d);

    /** \brief Number of dependencies */
    static inline int ndeps(unsigned char op);

    /** \brief Print */
    static inline void print(unsigned char op, std::ostream &stream, const std::string& x,
                             const std::string& y);
    static inline void printName(unsigned char op, std::ostream &stream);
    static inline void printPre(unsigned char op, std::ostream &stream);
    static inline void printSep(unsigned char op, std::ostream &stream);
    static inline void printPost(unsigned char op, std::ostream &stream);
  };

  /// Specialize the class so that it can be used with integer type
  template<>
  struct casadi_math<int>{

    /** \brief Evaluate a built in function */
    static inline void fun(unsigned char op, const int& x, const int& y, int& f) {
      double ff(0);
      casadi_math<double>::fun(op, static_cast<double>(x), static_cast<double>(y), ff);
      f = static_cast<int>(ff);
    }

    static inline void fun(unsigned char op, const int* x, const int* y, int* f, int n) {
      for (int i=0; i<n; ++i) {
        double ff(0);
        casadi_math<double>::fun(op, static_cast<double>(*x++), static_cast<double>(*y++), ff);
        *f++ = static_cast<int>(ff);
      }
    }

    static inline void fun(unsigned char op, const int* x, const int& y, int* f, int n) {
      for (int i=0; i<n; ++i) {
        double ff;
        casadi_math<double>::fun(op, static_cast<double>(*x++), static_cast<double>(y), ff);
        *f++ = static_cast<int>(ff);
      }
    }

    static inline void fun(unsigned char op, const int& x, const int* y, int* f, int n) {
      for (int i=0; i<n; ++i) {
        double ff;
        casadi_math<double>::fun(op, static_cast<double>(x), static_cast<double>(*y++), ff);
        *f++ = static_cast<int>(ff);
      }
    }

    /** \brief Evaluate a built in derivative function */
    static inline void der(unsigned char op, const int& x, const int& y, const int& f, int* d) {
      double d_real[2] = {static_cast<double>(d[0]), static_cast<double>(d[1])};
      casadi_math<double>::der(op, static_cast<double>(x), static_cast<double>(y),
                               static_cast<double>(f), d_real);
      d[0] = static_cast<int>(d_real[0]);
      d[1] = static_cast<int>(d_real[1]);
    }

    /** \brief Evaluate the function and the derivative function */
    static inline void derF(unsigned char op, const int& x, const int& y, int& f, int* d) {
      double d_real[2] = {static_cast<double>(d[0]), static_cast<double>(d[1])};
      double f_real(f);
      casadi_math<double>::derF(op, static_cast<double>(x), static_cast<double>(y), f_real, d_real);
      f = static_cast<int>(f_real);
      d[0] = static_cast<int>(d_real[0]);
      d[1] = static_cast<int>(d_real[1]);
    }

    /** \brief Number of dependencies */
    static inline int ndeps(unsigned char op) { return casadi_math<double>::ndeps(op);}

    /** \brief Print */
    static inline void print(unsigned char op, std::ostream &stream, const std::string& x,
                             const std::string& y) { casadi_math<double>::print(op, stream, x, y);}
    static inline void printPre(unsigned char op, std::ostream &stream) {
      casadi_math<double>::printPre(op, stream);}
    static inline void printName(unsigned char op, std::ostream &stream) {
      casadi_math<double>::printName(op, stream);}
    static inline void printSep(unsigned char op, std::ostream &stream) {
      casadi_math<double>::printSep(op, stream);}
    static inline void printPost(unsigned char op, std::ostream &stream) {
      casadi_math<double>::printPost(op, stream);}
  };

  // Template implementations

  template<typename T>
  inline void casadi_math<T>::fun(unsigned char op, const T& x, const T& y, T& f) {
    // NOTE: We define the implementation in a preprocessor macro to be able to force inlining,
    //  and to allow extensions in the VM
#define CASADI_MATH_FUN_BUILTIN_GEN(CNAME, X, Y, F, N)                  \
    case OP_ASSIGN:    CNAME<OP_ASSIGN>::fcn(X, Y, F, N);        break; \
  case OP_ADD:       CNAME<OP_ADD>::fcn(X, Y, F, N);           break;   \
  case OP_SUB:       CNAME<OP_SUB>::fcn(X, Y, F, N);           break;   \
  case OP_MUL:       CNAME<OP_MUL>::fcn(X, Y, F, N);           break;   \
  case OP_DIV:       CNAME<OP_DIV>::fcn(X, Y, F, N);           break;   \
  case OP_NEG:       CNAME<OP_NEG>::fcn(X, Y, F, N);           break;   \
  case OP_EXP:       CNAME<OP_EXP>::fcn(X, Y, F, N);           break;   \
  case OP_LOG:       CNAME<OP_LOG>::fcn(X, Y, F, N);           break;   \
  case OP_POW:       CNAME<OP_POW>::fcn(X, Y, F, N);           break;   \
  case OP_CONSTPOW:  CNAME<OP_CONSTPOW>::fcn(X, Y, F, N);      break;   \
  case OP_SQRT:      CNAME<OP_SQRT>::fcn(X, Y, F, N);          break;   \
  case OP_SQ:        CNAME<OP_SQ>::fcn(X, Y, F, N);            break;   \
  case OP_TWICE:     CNAME<OP_TWICE>::fcn(X, Y, F, N);         break;   \
  case OP_SIN:       CNAME<OP_SIN>::fcn(X, Y, F, N);           break;   \
  case OP_COS:       CNAME<OP_COS>::fcn(X, Y, F, N);           break;   \
  case OP_TAN:       CNAME<OP_TAN>::fcn(X, Y, F, N);           break;   \
  case OP_ASIN:      CNAME<OP_ASIN>::fcn(X, Y, F, N);          break;   \
  case OP_ACOS:      CNAME<OP_ACOS>::fcn(X, Y, F, N);          break;   \
  case OP_ATAN:      CNAME<OP_ATAN>::fcn(X, Y, F, N);          break;   \
  case OP_LT:        CNAME<OP_LT>::fcn(X, Y, F, N);            break;   \
  case OP_LE:        CNAME<OP_LE>::fcn(X, Y, F, N);            break;   \
  case OP_EQ:        CNAME<OP_EQ>::fcn(X, Y, F, N);            break;   \
  case OP_NE:        CNAME<OP_NE>::fcn(X, Y, F, N);            break;   \
  case OP_NOT:       CNAME<OP_NOT>::fcn(X, Y, F, N);           break;   \
  case OP_AND:       CNAME<OP_AND>::fcn(X, Y, F, N);           break;   \
  case OP_OR:        CNAME<OP_OR>::fcn(X, Y, F, N);            break;   \
  case OP_IF_ELSE_ZERO: CNAME<OP_IF_ELSE_ZERO>::fcn(X, Y, F, N); break; \
  case OP_FLOOR:     CNAME<OP_FLOOR>::fcn(X, Y, F, N);         break;   \
  case OP_CEIL:      CNAME<OP_CEIL>::fcn(X, Y, F, N);          break;   \
  case OP_FMOD:      CNAME<OP_FMOD>::fcn(X, Y, F, N);          break;   \
  case OP_FABS:      CNAME<OP_FABS>::fcn(X, Y, F, N);          break;   \
  case OP_SIGN:      CNAME<OP_SIGN>::fcn(X, Y, F, N);          break;   \
  case OP_COPYSIGN:  CNAME<OP_COPYSIGN>::fcn(X, Y, F, N);      break;   \
  case OP_ERF:       CNAME<OP_ERF>::fcn(X, Y, F, N);           break;   \
  case OP_FMIN:      CNAME<OP_FMIN>::fcn(X, Y, F, N);          break;   \
  case OP_FMAX:      CNAME<OP_FMAX>::fcn(X, Y, F, N);          break;   \
  case OP_INV:       CNAME<OP_INV>::fcn(X, Y, F, N);           break;   \
  case OP_SINH:      CNAME<OP_SINH>::fcn(X, Y, F, N);          break;   \
  case OP_COSH:      CNAME<OP_COSH>::fcn(X, Y, F, N);          break;   \
  case OP_TANH:      CNAME<OP_TANH>::fcn(X, Y, F, N);          break;   \
  case OP_ASINH:     CNAME<OP_ASINH>::fcn(X, Y, F, N);         break;   \
  case OP_ACOSH:     CNAME<OP_ACOSH>::fcn(X, Y, F, N);         break;   \
  case OP_ATANH:     CNAME<OP_ATANH>::fcn(X, Y, F, N);         break;   \
  case OP_ATAN2:     CNAME<OP_ATAN2>::fcn(X, Y, F, N);         break;   \
  case OP_ERFINV:    CNAME<OP_ERFINV>::fcn(X, Y, F, N);        break;   \
  case OP_LIFT:      CNAME<OP_LIFT>::fcn(X, Y, F, N);          break;   \
  case OP_PRINTME:   CNAME<OP_PRINTME>::fcn(X, Y, F, N);       break;

#define CASADI_MATH_FUN_BUILTIN(X, Y, F) CASADI_MATH_FUN_BUILTIN_GEN(BinaryOperationSS, X, Y, F, 1)

    switch (op) {
    CASADI_MATH_FUN_BUILTIN(x, y, f)
      }
  }

  template<typename T>
  inline void casadi_math<T>::fun(unsigned char op, const T* x, const T* y, T* f, int n) {
    switch (op) {
      CASADI_MATH_FUN_BUILTIN_GEN(BinaryOperationVV, x, y, f, n)
        }
  }

  template<typename T>
  inline void casadi_math<T>::fun(unsigned char op, const T* x, const T& y, T* f, int n) {
    switch (op) {
      CASADI_MATH_FUN_BUILTIN_GEN(BinaryOperationVS, x, y, f, n)
        }
  }

  template<typename T>
  inline void casadi_math<T>::fun(unsigned char op, const T& x, const T* y, T* f, int n) {
    switch (op) {
      CASADI_MATH_FUN_BUILTIN_GEN(BinaryOperationSV, x, y, f, n)
        }
  }


  template<typename T>
  inline void casadi_math<T>::der(unsigned char op, const T& x, const T& y, const T& f, T* d) {
    // NOTE: We define the implementation in a preprocessor macro to be able to force inlining,
    // and to allow extensions in the VM
#define CASADI_MATH_DER_BUILTIN(X, Y, F, D)                             \
    case OP_ASSIGN:    BinaryOperation<OP_ASSIGN>::der(X, Y, F, D);     break; \
  case OP_ADD:       BinaryOperation<OP_ADD>::der(X, Y, F, D);        break; \
  case OP_SUB:       BinaryOperation<OP_SUB>::der(X, Y, F, D);        break; \
  case OP_MUL:       BinaryOperation<OP_MUL>::der(X, Y, F, D);        break; \
  case OP_DIV:       BinaryOperation<OP_DIV>::der(X, Y, F, D);        break; \
  case OP_NEG:       BinaryOperation<OP_NEG>::der(X, Y, F, D);        break; \
  case OP_EXP:       BinaryOperation<OP_EXP>::der(X, Y, F, D);        break; \
  case OP_LOG:       BinaryOperation<OP_LOG>::der(X, Y, F, D);        break; \
  case OP_POW:       BinaryOperation<OP_POW>::der(X, Y, F, D);        break; \
  case OP_CONSTPOW:  BinaryOperation<OP_CONSTPOW>::der(X, Y, F, D);   break; \
  case OP_SQRT:      BinaryOperation<OP_SQRT>::der(X, Y, F, D);       break; \
  case OP_SQ:        BinaryOperation<OP_SQ>::der(X, Y, F, D);         break; \
  case OP_TWICE:     BinaryOperation<OP_TWICE>::der(X, Y, F, D);      break; \
  case OP_SIN:       BinaryOperation<OP_SIN>::der(X, Y, F, D);        break; \
  case OP_COS:       BinaryOperation<OP_COS>::der(X, Y, F, D);        break; \
  case OP_TAN:       BinaryOperation<OP_TAN>::der(X, Y, F, D);        break; \
  case OP_ASIN:      BinaryOperation<OP_ASIN>::der(X, Y, F, D);       break; \
  case OP_ACOS:      BinaryOperation<OP_ACOS>::der(X, Y, F, D);       break; \
  case OP_ATAN:      BinaryOperation<OP_ATAN>::der(X, Y, F, D);       break; \
  case OP_LT:        BinaryOperation<OP_LT>::der(X, Y, F, D);         break; \
  case OP_LE:        BinaryOperation<OP_LE>::der(X, Y, F, D);         break; \
  case OP_EQ:        BinaryOperation<OP_EQ>::der(X, Y, F, D);         break; \
  case OP_NE:        BinaryOperation<OP_NE>::der(X, Y, F, D);         break; \
  case OP_NOT:       BinaryOperation<OP_NOT>::der(X, Y, F, D);        break; \
  case OP_AND:       BinaryOperation<OP_AND>::der(X, Y, F, D);        break; \
  case OP_OR:        BinaryOperation<OP_OR>::der(X, Y, F, D);         break; \
  case OP_IF_ELSE_ZERO: BinaryOperation<OP_IF_ELSE_ZERO>::der(X, Y, F, D);         break; \
  case OP_FLOOR:     BinaryOperation<OP_FLOOR>::der(X, Y, F, D);      break; \
  case OP_CEIL:      BinaryOperation<OP_CEIL>::der(X, Y, F, D);       break; \
  case OP_FMOD:      BinaryOperation<OP_FMOD>::der(X, Y, F, D);       break; \
  case OP_FABS:      BinaryOperation<OP_FABS>::der(X, Y, F, D);       break; \
  case OP_SIGN:      BinaryOperation<OP_SIGN>::der(X, Y, F, D);       break; \
  case OP_COPYSIGN:  BinaryOperation<OP_COPYSIGN>::der(X, Y, F, D);   break; \
  case OP_ERF:       BinaryOperation<OP_ERF>::der(X, Y, F, D);        break; \
  case OP_FMIN:      BinaryOperation<OP_FMIN>::der(X, Y, F, D);       break; \
  case OP_FMAX:      BinaryOperation<OP_FMAX>::der(X, Y, F, D);       break; \
  case OP_INV:       BinaryOperation<OP_INV>::der(X, Y, F, D);        break; \
  case OP_SINH:      BinaryOperation<OP_SINH>::der(X, Y, F, D);       break; \
  case OP_COSH:      BinaryOperation<OP_COSH>::der(X, Y, F, D);       break; \
  case OP_TANH:      BinaryOperation<OP_TANH>::der(X, Y, F, D);       break; \
  case OP_ASINH:     BinaryOperation<OP_ASINH>::der(X, Y, F, D);      break; \
  case OP_ACOSH:     BinaryOperation<OP_ACOSH>::der(X, Y, F, D);      break; \
  case OP_ATANH:     BinaryOperation<OP_ATANH>::der(X, Y, F, D);      break; \
  case OP_ATAN2:     BinaryOperation<OP_ATAN2>::der(X, Y, F, D);      break; \
  case OP_ERFINV:    BinaryOperation<OP_ERFINV>::der(X, Y, F, D);     break; \
  case OP_LIFT:      BinaryOperation<OP_LIFT>::der(X, Y, F, D);       break; \
  case OP_PRINTME:   BinaryOperation<OP_PRINTME>::der(X, Y, F, D);    break;

    switch (op) {
    CASADI_MATH_DER_BUILTIN(x, y, f, d)
      }
  }


    template<typename T>
      inline void casadi_math<T>::derF(unsigned char op, const T& x, const T& y, T& f, T* d) {
    // NOTE: We define the implementation in a preprocessor macro to be able to force inlining,
    // and to allow extensions in the VM
#define CASADI_MATH_DERF_BUILTIN(X, Y, F, D)                            \
    case OP_ASSIGN:    DerBinaryOpertion<OP_ASSIGN>::derf(X, Y, F, D);        break; \
  case OP_ADD:       DerBinaryOpertion<OP_ADD>::derf(X, Y, F, D);        break; \
  case OP_SUB:       DerBinaryOpertion<OP_SUB>::derf(X, Y, F, D);        break; \
  case OP_MUL:       DerBinaryOpertion<OP_MUL>::derf(X, Y, F, D);        break; \
  case OP_DIV:       DerBinaryOpertion<OP_DIV>::derf(X, Y, F, D);        break; \
  case OP_NEG:       DerBinaryOpertion<OP_NEG>::derf(X, Y, F, D);        break; \
  case OP_EXP:       DerBinaryOpertion<OP_EXP>::derf(X, Y, F, D);        break; \
  case OP_LOG:       DerBinaryOpertion<OP_LOG>::derf(X, Y, F, D);        break; \
  case OP_POW:       DerBinaryOpertion<OP_POW>::derf(X, Y, F, D);        break; \
  case OP_CONSTPOW:  DerBinaryOpertion<OP_CONSTPOW>::derf(X, Y, F, D);   break; \
  case OP_SQRT:      DerBinaryOpertion<OP_SQRT>::derf(X, Y, F, D);       break; \
  case OP_SQ:        DerBinaryOpertion<OP_SQ>::derf(X, Y, F, D);         break; \
  case OP_TWICE:     DerBinaryOpertion<OP_TWICE>::derf(X, Y, F, D);      break; \
  case OP_SIN:       DerBinaryOpertion<OP_SIN>::derf(X, Y, F, D);        break; \
  case OP_COS:       DerBinaryOpertion<OP_COS>::derf(X, Y, F, D);        break; \
  case OP_TAN:       DerBinaryOpertion<OP_TAN>::derf(X, Y, F, D);        break; \
  case OP_ASIN:      DerBinaryOpertion<OP_ASIN>::derf(X, Y, F, D);       break; \
  case OP_ACOS:      DerBinaryOpertion<OP_ACOS>::derf(X, Y, F, D);       break; \
  case OP_ATAN:      DerBinaryOpertion<OP_ATAN>::derf(X, Y, F, D);       break; \
  case OP_LT:        DerBinaryOpertion<OP_LT>::derf(X, Y, F, D);         break; \
  case OP_LE:        DerBinaryOpertion<OP_LE>::derf(X, Y, F, D);         break; \
  case OP_EQ:        DerBinaryOpertion<OP_EQ>::derf(X, Y, F, D);         break; \
  case OP_NE:        DerBinaryOpertion<OP_NE>::derf(X, Y, F, D);         break; \
  case OP_NOT:       DerBinaryOpertion<OP_NOT>::derf(X, Y, F, D);        break; \
  case OP_AND:       DerBinaryOpertion<OP_AND>::derf(X, Y, F, D);        break; \
  case OP_OR:        DerBinaryOpertion<OP_OR>::derf(X, Y, F, D);         break; \
  case OP_IF_ELSE_ZERO: DerBinaryOpertion<OP_IF_ELSE_ZERO>::derf(X, Y, F, D);         break; \
  case OP_FLOOR:     DerBinaryOpertion<OP_FLOOR>::derf(X, Y, F, D);      break; \
  case OP_CEIL:      DerBinaryOpertion<OP_CEIL>::derf(X, Y, F, D);       break; \
  case OP_FMOD:      DerBinaryOpertion<OP_FMOD>::derf(X, Y, F, D);       break; \
  case OP_FABS:      DerBinaryOpertion<OP_FABS>::derf(X, Y, F, D);       break; \
  case OP_SIGN:      DerBinaryOpertion<OP_SIGN>::derf(X, Y, F, D);       break; \
  case OP_COPYSIGN:  DerBinaryOpertion<OP_COPYSIGN>::derf(X, Y, F, D);   break; \
  case OP_ERF:       DerBinaryOpertion<OP_ERF>::derf(X, Y, F, D);        break; \
  case OP_FMIN:      DerBinaryOpertion<OP_FMIN>::derf(X, Y, F, D);       break; \
  case OP_FMAX:      DerBinaryOpertion<OP_FMAX>::derf(X, Y, F, D);       break; \
  case OP_INV:       DerBinaryOpertion<OP_INV>::derf(X, Y, F, D);        break; \
  case OP_SINH:      DerBinaryOpertion<OP_SINH>::derf(X, Y, F, D);       break; \
  case OP_COSH:      DerBinaryOpertion<OP_COSH>::derf(X, Y, F, D);       break; \
  case OP_TANH:      DerBinaryOpertion<OP_TANH>::derf(X, Y, F, D);       break; \
  case OP_ASINH:     DerBinaryOpertion<OP_ASINH>::derf(X, Y, F, D);      break; \
  case OP_ACOSH:     DerBinaryOpertion<OP_ACOSH>::derf(X, Y, F, D);      break; \
  case OP_ATANH:     DerBinaryOpertion<OP_ATANH>::derf(X, Y, F, D);      break; \
  case OP_ATAN2:     DerBinaryOpertion<OP_ATAN2>::derf(X, Y, F, D);      break; \
  case OP_ERFINV:    DerBinaryOpertion<OP_ERFINV>::derf(X, Y, F, D);     break; \
  case OP_LIFT:      DerBinaryOpertion<OP_LIFT>::derf(X, Y, F, D);       break; \
  case OP_PRINTME:   DerBinaryOpertion<OP_PRINTME>::derf(X, Y, F, D);    break;

    switch (op) {
      CASADI_MATH_DERF_BUILTIN(x, y, f, d)
        }
  }

  template<typename T>
  inline int casadi_math<T>::ndeps(unsigned char op) {
#define CASADI_MATH_BINARY_BUILTIN              \
    case OP_ADD:                                \
  case OP_SUB:                                  \
  case OP_MUL:                                  \
  case OP_DIV:                                  \
  case OP_POW:                                  \
  case OP_CONSTPOW:                             \
  case OP_LT:                                   \
  case OP_LE:                                   \
  case OP_EQ:                                   \
  case OP_NE:                                   \
  case OP_AND:                                  \
  case OP_OR:                                   \
  case OP_IF_ELSE_ZERO:                         \
  case OP_COPYSIGN:                             \
  case OP_FMIN:                                 \
  case OP_FMAX:                                 \
  case OP_ATAN2:                                \
  case OP_PRINTME:                              \
  case OP_LIFT:

    switch (op) {
  case OP_CONST:
  case OP_PARAMETER:
  case OP_INPUT:
    return 0;
    CASADI_MATH_BINARY_BUILTIN
      return 2;
  default:
    return 1;
  }
  }

    template<typename T>
      inline void casadi_math<T>::print(unsigned char op, std::ostream &stream, const std::string& x,
      const std::string& y) {
    if (ndeps(op)==2) {
    printPre(op, stream);
    stream << x;
    printSep(op, stream);
    stream << y;
    printPost(op, stream);
  } else {
    printPre(op, stream);
    stream << x;
    printPost(op, stream);
  }
  }

    template<typename T>
      inline void casadi_math<T>::printName(unsigned char op, std::ostream &stream) {
    switch (op) {
  case OP_ASSIGN:              stream << "assign";              break;
  case OP_ADD:                 stream << "add";                 break;
  case OP_SUB:                 stream << "sub";                 break;
  case OP_MUL:                 stream << "mul";                 break;
  case OP_DIV:                 stream << "div";                 break;
  case OP_NEG:                 stream << "neg";                 break;
  case OP_EXP:                 stream << "exp";                 break;
  case OP_LOG:                 stream << "log";                 break;
  case OP_CONSTPOW:
  case OP_POW:                 stream << "pow";                 break;
  case OP_SQRT:                stream << "sqrt";                break;
  case OP_SQ:                  stream << "sq";                  break;
  case OP_TWICE:               stream << "twice";               break;
  case OP_SIN:                 stream << "sin";                 break;
  case OP_COS:                 stream << "cos";                 break;
  case OP_TAN:                 stream << "tan";                 break;
  case OP_ASIN:                stream << "asin";                break;
  case OP_ACOS:                stream << "acos";                break;
  case OP_ATAN:                stream << "atan";                break;
  case OP_LT:                  stream << "lt";                  break;
  case OP_LE:                  stream << "le";                  break;
  case OP_EQ:                  stream << "eq";                  break;
  case OP_NE:                  stream << "ne";                  break;
  case OP_NOT:                 stream << "not";                 break;
  case OP_AND:                 stream << "and";                 break;
  case OP_OR:                  stream << "or";                  break;
  case OP_FLOOR:               stream << "floor";               break;
  case OP_CEIL:                stream << "ceil";                break;
  case OP_FMOD:                stream << "fmod";                break;
  case OP_FABS:                stream << "fabs";                break;
  case OP_SIGN:                stream << "sign";                break;
  case OP_COPYSIGN:            stream << "copysign";            break;
  case OP_IF_ELSE_ZERO:        stream << "if_else_zero";        break;
  case OP_ERF:                 stream << "erf";                 break;
  case OP_FMIN:                stream << "fmin";                break;
  case OP_FMAX:                stream << "fmax";                break;
  case OP_INV:                 stream << "inv";                 break;
  case OP_SINH:                stream << "sinh";                break;
  case OP_COSH:                stream << "cosh";                break;
  case OP_TANH:                stream << "tanh";                break;
  case OP_ASINH:               stream << "asinh";               break;
  case OP_ACOSH:               stream << "acosh";               break;
  case OP_ATANH:               stream << "atanh";               break;
  case OP_ATAN2:               stream << "atan2";               break;
  case OP_CONST:               stream << "const";               break;
  case OP_INPUT:               stream << "input";               break;
  case OP_OUTPUT:              stream << "output";              break;
  case OP_PARAMETER:           stream << "parameter";           break;
  case OP_CALL:                stream << "call";                break;
  case OP_MATMUL:              stream << "matmul";              break;
  case OP_SOLVE:               stream << "solve";               break;
  case OP_TRANSPOSE:           stream << "transpose";           break;
  case OP_DETERMINANT:         stream << "determinant";         break;
  case OP_INVERSE:             stream << "inverse";             break;
  case OP_DOT:          stream << "dot";          break;
  case OP_HORZCAT:             stream << "horzcat";             break;
  case OP_VERTCAT:             stream << "vertcat";             break;
  case OP_DIAGCAT:             stream << "diagcat";             break;
  case OP_HORZSPLIT:           stream << "horzsplit";           break;
  case OP_VERTSPLIT:           stream << "vertsplit";           break;
  case OP_DIAGSPLIT:           stream << "diagsplit";           break;
  case OP_RESHAPE:             stream << "reshape";             break;
  case OP_SUBREF:              stream << "subref";              break;
  case OP_SUBASSIGN:           stream << "subassign";           break;
  case OP_GETNONZEROS:         stream << "getnonzeros";         break;
  case OP_ADDNONZEROS:         stream << "addnonzeros";         break;
  case OP_SETNONZEROS:         stream << "setnonzeros";         break;
  case OP_PROJECT:             stream << "project";             break;
  case OP_ASSERTION:           stream << "assertion";           break;
  case OP_NORM2:               stream << "norm2";               break;
  case OP_NORM1:               stream << "norm1";               break;
  case OP_NORMINF:             stream << "norminf";             break;
  case OP_NORMF:               stream << "normf";               break;
  case OP_ERFINV:              stream << "erfinv";              break;
  case OP_PRINTME:             stream << "printme";             break;
  case OP_LIFT:                stream << "lift";                break;
  }
    }

  template<typename T>
  inline void casadi_math<T>::printPre(unsigned char op, std::ostream &stream) {
    switch (op) {
    case OP_ASSIGN:                          break;
    case OP_ADD:       stream << "(";        break;
    case OP_SUB:       stream << "(";        break;
    case OP_MUL:       stream << "(";        break;
    case OP_DIV:       stream << "(";        break;
    case OP_NEG:       stream << "(-";       break;
    case OP_TWICE:     stream << "(2.*";     break;
    case OP_LT:        stream << "(";        break;
    case OP_LE:        stream << "(";        break;
    case OP_EQ:        stream << "(";        break;
    case OP_NE:        stream << "(";        break;
    case OP_NOT:       stream << "(!";       break;
    case OP_AND:       stream << "(";        break;
    case OP_OR:        stream << "(";        break;
    case OP_IF_ELSE_ZERO: stream << "(";        break;
    case OP_INV:       stream << "(1./";     break;
    default:           printName(op, stream); stream << "("; break;
    }
  }

  template<typename T>
  inline void casadi_math<T>::printSep(unsigned char op, std::ostream &stream) {
    switch (op) {
    case OP_ADD:       stream << "+";        break;
    case OP_SUB:       stream << "-";        break;
    case OP_MUL:       stream << "*";        break;
    case OP_DIV:       stream << "/";        break;
    case OP_LT:        stream << "<";        break;
    case OP_LE:        stream << "<=";       break;
    case OP_EQ:        stream << "==";       break;
    case OP_NE:        stream << "!=";       break;
    case OP_AND:       stream << "&&";       break;
    case OP_OR:        stream << "||";       break;
    case OP_IF_ELSE_ZERO: stream << "?";     break;
    default:           stream << ",";        break;
    }
  }

  template<typename T>
  inline void casadi_math<T>::printPost(unsigned char op, std::ostream &stream) {
    switch (op) {
    case OP_ASSIGN:                        break;
    case OP_IF_ELSE_ZERO: stream << ":0)"; break;
    default:        stream << ")";         break;
    }
  }

#endif // SWIG

} // namespace casadi

/// \endcond

#endif // CASADI_CALCULUS_HPP
