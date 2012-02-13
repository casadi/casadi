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

#ifndef CASADI_CALCULUS_HPP
#define CASADI_CALCULUS_HPP

#include <iostream>
#include <string>
#include <cmath>
#include <limits>
#include <algorithm>
#include "casadi_exception.hpp"

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

#pragma warning (disable:4996)

// warning C4018: '<' : signed/unsigned mismatch
#pragma warning (disable:4018)

// warning C4800: 'int' : forcing value to bool 'true'or 'false'(performance warning)
#pragma warning (disable:4800)
#endif

namespace CasADi{

  //@{
  /** \brief  Pre-C99 elementary functions from the math.h (cmath) header */
  template<class T> T sqrt(const T &x){return x.sqrt();}
  using std::sqrt;

  template<class T> T sin(const T &x){return x.sin();}
  using std::sin;
  
  template<class T> T cos(const T &x){return x.cos();}
  using std::cos;

  template<class T> T tan(const T &x){return x.tan();}
  using std::tan;

  template<class T> T atan(const T &x){return x.arctan();}
  using std::atan;

  template<class T> T asin(const T &x){return x.arcsin();}
  using std::asin;

  template<class T> T acos(const T &x){return x.arccos();}
  using std::acos;

  template<class T> T sinh(const T &x){return x.sinh();}
  using std::sinh;

  template<class T> T cosh(const T &x){return x.cosh();}
  using std::cosh;

  template<class T> T tanh(const T &x){return x.tanh();}
  using std::tanh;

  template<class T> T exp(const T &x){return x.exp();}
  using std::exp;

  template<class T> T log(const T &x){return x.log();}
  using std::log;

  template<class T> T log10(const T &x){return x.log10();}
  using std::log10;

  template<class T> T pow(const T &x, const T &n){ return x.__pow__(n);}
  template<class T> T pow(const T &x,   double n){ return x.__pow__(n);}
  template<class T> T pow(double   x, const T &n){ return T(x).__pow__(n);}
  using std::pow;

  template<class T> T abs(const T &x){return x.fabs();}
  using std::abs;
    
  template<class T> T fabs(const T &x){return x.fabs();}
  using std::fabs;
  
  template<class T> T floor(const T &x){return x.floor();}
  using std::floor;
  
  template<class T> T ceil(const T &x){return x.ceil();}
  using std::ceil;
  //@}

  //@{
  /** \brief  C99 elementary functions from the math.h header */
  template<class T> T erf(const T &x){return x.erf();}
  #ifdef HAS_ERF
  using ::erf;
  #else // HAS ERF
  inline double erf(double x) throw(){
    // Approximation found in sourceforge and modified, originally from numerical recepies in fortran
    double sx = x<0 ? -1 : x>0 ? 1 : x;
    double z = sx*x;
    double t = 1.0/(1.0+0.5*z);
    return 1.-sx*(t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
    t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
    t*(-0.82215223+t*0.17087277))))))))));
  }
  #endif // HAS ERF
  
  template<class T> T fmin(const T &x, const T &n){ return x.fmin(n);}
  template<class T> T fmin(const T &x,   double n){ return x.fmin(n);}
  template<class T> T fmin(double   x, const T &n){ return T(x).fmin(n);}
  inline double fmin(double x, double y) throw(){ return std::min(x,y);}

  template<class T> T fmax(const T &x, const T &n){ return x.fmax(n);}
  template<class T> T fmax(const T &x,   double n){ return x.fmax(n);}
  template<class T> T fmax(double   x, const T &n){ return T(x).fmax(n);}
  inline double fmax(double x, double y) throw(){ return std::max(x,y);}

  inline int isnan(double x) throw(){ return x!=x;}
  inline int isinf(double x) throw(){ return isnan(x-x);}
  //@}

  //@{
  /** \brief  CasADi additions */
  template<class T> T constpow(const T &x, const T &n){ return x.constpow(n);}
  
  template<class T> T printme(const T &x, const T &y){ return x.printme(y);}
  inline double printme(double x, double y){ 
    std::cout << "|> " << y << " : " << x << std::endl;
    return x;
  }
  
  /// Sign function, note that sign(nan) == nan
  template<class T> T sign(const T &x){return x.sign();}

  /// Sign function, note that sign(nan) == nan
  inline double sign(double x){ return x<0 ? -1 : x>0 ? 1 : x;}

  /// Inverse of the error function
  template<class T> T erfinv(const T &x){return x.erfinv();}
  #ifdef HAS_ERFINV
  using ::erfinv;
  #else // HAS ERFINV
  inline double erfinv(double x) throw(){
    // Approximation found in sourceforge and modified: Not very efficent
    if(x>=1){
      return x==1 ? std::numeric_limits<double>::infinity() : std::numeric_limits<double>::quiet_NaN();
    } else if(x<=-1){
      return x==-1 ? -std::numeric_limits<double>::infinity() : std::numeric_limits<double>::quiet_NaN();
    } else if(x<-0.7){
        double z = sqrt(-log((1.0+x)/2.0));
        return -(((1.641345311*z+3.429567803)*z-1.624906493)*z-1.970840454)/((1.637067800*z+3.543889200)*z+1.0);
    } else {
      double y;
      if(x<0.7){
        double z = x*x;
        y = x*(((-0.140543331*z+0.914624893)*z-1.645349621)*z+0.886226899)/((((-0.329097515*z+0.012229801)*z+1.442710462)*z-2.118377725)*z+1.0);
      } else {
        double z = sqrt(-log((1.0-x)/2.0));
        y = (((1.641345311*z+3.429567803)*z-1.624906493)*z-1.970840454)/((1.637067800*z+3.543889200)*z+1.0);
      }
      
      //polish x to full accuracy
      y = y - (erf(y) - x) / (2.0/sqrt(M_PI) * exp(-y*y));
      y = y - (erf(y) - x) / (2.0/sqrt(M_PI) * exp(-y*y));
      return y;
    }
  }
  #endif // HAS_ERFINV
  //@}
}

namespace CasADi{

template<typename T>
T timesTwo(const T& x){
  return x+x;
}
  
template<typename T>
T square(const T& x){
  return x*x;
}
  
template<int I>
class UnaryOperation{
  public:
    /// Function evaluation
    template<typename T> static inline void fcn(const T& x, T& f);
    
    /// Partial derivatives
    template<typename T> static inline void der(const T& x, const T& f, T* d);
};

template<int I>
class BinaryOperation{
  public:
    /// Function evaluation
    template<typename T> static inline void fcn(const T& x, const T& y, T& f){ UnaryOperation<I>::fcn(x,f);}
    
    /// Partial derivatives - binary function
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d){ UnaryOperation<I>::der(x,f,d); d[1]=0; }
};

template<int I>
class BinaryOperationE{
  public:
    /// Function evaluation
    template<typename T> static inline T fcn(const T& x, const T& y){ 
      T ret;
      BinaryOperation<I>::fcn(x,y,ret);
      return ret;
    }
};

template<int I>
class AddBinaryOperation{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f){ f+= BinaryOperationE<I>::fcn(x,y);}
};

template<int I>
class SubBinaryOperation{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f){ f-= BinaryOperationE<I>::fcn(x,y);}
};

template<int I>
class MulBinaryOperation{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f){ f*= BinaryOperationE<I>::fcn(x,y);}
};

template<int I>
class DivBinaryOperation{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f){ f/= BinaryOperationE<I>::fcn(x,y);}
};

/// Enum for quick access to any node
enum Operation{
  ASSIGN,
  ADD,  SUB,  MUL,  DIV,
  NEG,  EXP,  LOG,  POW, CONSTPOW,
  SQRT,  SIN,  COS,  TAN,  
  ASIN,  ACOS,  ATAN,  
  STEP,  
  FLOOR,  CEIL,  EQUALITY, FABS, SIGN, 
  ERF,  FMIN,  FMAX,
  INV,
  SINH,  COSH,  TANH,
  ERFINV,
  PRINTME,
  NUM_BUILT_IN_OPS
};

/// Simple assignment
template<>
class UnaryOperation<ASSIGN>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f){ f = x;}
    template<typename T> static inline void der(const T& x, const T& f, T* d){ d[0] = 1; }
};

/// Addition
template<>
class BinaryOperation<ADD>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f){ f = x+y;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d){ d[0]=d[1]=1;}
};

/// Subtraction
template<>
class BinaryOperation<SUB>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f){ f = x-y;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d){ d[0]=1; d[1]=-1;}
};

/// Multiplication
template<>
class BinaryOperation<MUL>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f){ f = x*y;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d){ d[0]=y; d[1]=x;}
};

/// Division
template<>
class BinaryOperation<DIV>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f){ f = x/y;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d){ d[0]=1/y; d[1]=-f/y;}
};

/// Negation
template<>
class UnaryOperation<NEG>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f){ f = -x;}
    template<typename T> static inline void der(const T& x, const T& f, T* d){ d[0]=-1;}
};

/// Natural exponent
template<>
class UnaryOperation<EXP>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f){ f = exp(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d){ d[0]=f;}
};

/// Natural logarithm
template<>
class UnaryOperation<LOG>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f){ f = log(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d){ d[0]=1/x;}
};

/// Power, defined only for x>=0
template<>
class BinaryOperation<POW>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f){ f = pow(x,y);}
    // See issue #104 why d[0] is no longer y*f/x
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d){ d[0]=y*pow(x,y-1); d[1]=log(x)*f;}
};

/// Power, defined only for y constant
template<>
class BinaryOperation<CONSTPOW>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f){ f = pow(x,y);}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d){ d[0]=y*pow(x,y-1); d[1]=0;}
};

/// Square root
template<>
class UnaryOperation<SQRT>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f){ f = sqrt(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d){ d[0]=1/(timesTwo(f));}
};

/// Sine
template<>
class UnaryOperation<SIN>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f){ f = sin(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d){ d[0]=cos(x);}
};

/// Cosine
template<>
class UnaryOperation<COS>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f){ f = cos(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d){ d[0]=-sin(x);}
};

/// Tangent
template<>
class UnaryOperation<TAN>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f){ f = tan(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d){ d[0] = 1/square(cos(x));}
};

/// Arcus sine
template<>
class UnaryOperation<ASIN>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f){ f = asin(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d){ d[0]=1/sqrt(1-x*x);}
};

/// Arcus cosine
template<>
class UnaryOperation<ACOS>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f){ f = acos(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d){ d[0]=-1/sqrt(1-x*x);}
};

/// Arcus tangent
template<>
class UnaryOperation<ATAN>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f){ f = atan(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d){ d[0] = 1/(1+x*x);}
};

/// Step function
template<>
class UnaryOperation<STEP>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f){ f = x >= T(0.);}
    template<typename T> static inline void der(const T& x, const T& f, T* d){ d[0] = 0;}
};

/// Floor function
template<>
class UnaryOperation<FLOOR>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f){ f = floor(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d){ d[0] = 0;}
};

/// Ceil function
template<>
class UnaryOperation<CEIL>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f){ f = ceil(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d){ d[0] = 0;}
};

/// Equality
template<>
class BinaryOperation<EQUALITY>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f){ f = x==y;}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d){ d[0]=d[1]=0;}
};

/// Error function
template<>
class UnaryOperation<ERF>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f){ f = erf(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d){ d[0] = (2/sqrt(M_PI))*exp(-x*x);}
};

/// Absolute value
template<>
class UnaryOperation<FABS>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f){ f = fabs(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d){ d[0]=sign(x);}
};

/// Sign
template<>
class UnaryOperation<SIGN>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f){ f = sign(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d){ d[0]=0;}
};

/// Minimum
template<>
class BinaryOperation<FMIN>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f){ f = fmin(x,y);}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d){ d[0]=x<=y; d[1]=!d[0];}
};

/// Maximum
template<>
class BinaryOperation<FMAX>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f){ f = fmax(x,y);}
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d){ d[0]=x>=y; d[1]=!d[0];}
};

/// Elementwise inverse
template<>
class UnaryOperation<INV>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f){ f = 1./x;}
    template<typename T> static inline void der(const T& x, const T& f, T* d){ d[0] = -f*f; }
};

/// Hyperbolic sine
template<>
class UnaryOperation<SINH>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f){ f = sinh(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d){ d[0] = cosh(x); }
};

/// Hyperbolic cosine
template<>
class UnaryOperation<COSH>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f){ f = cosh(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d){ d[0] = -sinh(x); }
};

/// Hyperbolic tangent
template<>
class UnaryOperation<TANH>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f){ f = tanh(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d){ d[0] = 1-f*f; }
};

/// Inverse of error function
template<>
class UnaryOperation<ERFINV>{
  public:
    template<typename T> static inline void fcn(const T& x, T& f){ f = erfinv(x);}
    template<typename T> static inline void der(const T& x, const T& f, T* d){ d[0] = (sqrt(M_PI)/2)*exp(f*f); }
};

/// Identity operator with the side effect of printing
template<>
class BinaryOperation<PRINTME>{
  public:
    template<typename T> static inline void fcn(const T& x, const T& y, T& f){f = printme(x,y); }
    template<typename T> static inline void der(const T& x, const T& y, const T& f, T* d){ d[0]=1; d[1]=0;}
};

} // namespace CasADi

#endif //CASADI_CALCULUS_HPP
