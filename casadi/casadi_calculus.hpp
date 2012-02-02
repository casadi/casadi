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
#include "casadi_exception.hpp"
#include "pre_c99_support.hpp"

namespace CasADi{

  //@{
  /** \brief  Pre-C99 elementary functions from the math.h (cmath) header */
  template<class T> T sqrt(const T &x){return x.sqrt();}
  using std::sqrt;
  inline double sqrt(int x){return sqrt(double(x));}

  template<class T> T sin(const T &x){return x.sin();}
  using std::sin;
  inline double sin(int x){return sin(double(x));}
  
  template<class T> T cos(const T &x){return x.cos();}
  using std::cos;
  inline double cos(int x){return cos(double(x));}

  template<class T> T tan(const T &x){return x.tan();}
  using std::tan;
  inline double tan(int x){return tan(double(x));}

  template<class T> T atan(const T &x){return x.arctan();}
  using std::atan;
  inline double atan(int x){return atan(double(x));}

  template<class T> T asin(const T &x){return x.arcsin();}
  using std::asin;
  inline double asin(int x){return asin(double(x));}

  template<class T> T acos(const T &x){return x.arccos();}
  using std::acos;
  inline double acos(int x){return acos(double(x));}

  template<class T> T sinh(const T &x){return x.sinh();}
  using std::sinh;
  inline double sinh(int x){return sinh(double(x));}

  template<class T> T cosh(const T &x){return x.cosh();}
  using std::cosh;
  inline double cosh(int x){return cosh(double(x));}

  template<class T> T tanh(const T &x){return x.tanh();}
  using std::tanh;
  inline double tanh(int x){return tanh(double(x));}

  template<class T> T exp(const T &x){return x.exp();}
  using std::exp;
  inline double exp(int x){return exp(double(x));}

  template<class T> T log(const T &x){return x.log();}
  using std::log;
  inline double log(int x){return log(double(x));}

  template<class T> T log10(const T &x){return x.log10();}
  using std::log10;
  inline double log10(int x){return log10(double(x));}

  template<class T> T pow(const T &x, const T &n){ return x.__pow__(n);}
  template<class T> T pow(const T &x,   double n){ return x.__pow__(n);}
  template<class T> T pow(double   x, const T &n){ return T(x).__pow__(n);}
  using std::pow;
  inline double pow(int x, int y){return pow(double(x),double(y));}

  template<class T> T abs(const T &x){return x.fabs();}
  using std::abs;
    
  template<class T> T fabs(const T &x){return x.fabs();}
  using std::fabs;
  inline int fabs(int x){return x>=0 ? x : -x; }
  
  template<class T> T floor(const T &x){return x.floor();}
  using std::floor;
  inline int floor(int x){return x; }
  
  template<class T> T ceil(const T &x){return x.ceil();}
  using std::ceil;
  inline int ceil(int x){return x; }
  //@}

  //@{
  /** \brief  C99 elementary functions from the math.h header */
  template<class T> T erf(const T &x){return x.erf();}
  using ::erf;
  inline double erf(int x){return erf(double(x));}
  
  template<class T> T fmin(const T &x, const T &n){ return x.fmin(n);}
  template<class T> T fmin(const T &x,   double n){ return x.fmin(n);}
  template<class T> T fmin(double   x, const T &n){ return T(x).fmin(n);}
  using ::fmin;
  inline int fmin(int x, int y){return x<=y ? x : y;}

  template<class T> T fmax(const T &x, const T &n){ return x.fmax(n);}
  template<class T> T fmax(const T &x,   double n){ return x.fmax(n);}
  template<class T> T fmax(double   x, const T &n){ return T(x).fmax(n);}
  using ::fmax;
  inline int fmax(int x, int y){return x>=y ? x : y;}
  //@}

  //@{
  /** \brief  CasADi additions */
  template<class T> T constpow(const T &x, const T &n){ return x.constpow(n);}
  
  template<class T> T printme(const T &x, const T &y){ return x.printme(y);}
  inline int printme(int x, int y){ return x;}

  template<class T> T sign(const T &x){return x.sign();}
  inline double sign(double x){ return x<0 ? -1 : x>0 ? 1 : x;} // NOTE: sign(nan) == nan
  inline int sign(int x){ return x<0 ? -1 : x>0 ? 1 : x; }

  template<class T> T erfinv(const T &x){return x.erfinv();}
  inline double erfinv(int x){ return erfinv(double(x));}
  
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
    template<typename T> inline static void fcn(const T& x, T& f);
    
    /// Partial derivatives
    template<typename T> inline static void der(const T& x, const T& f, T* d);
};

template<int I>
class BinaryOperation{
  public:
    /// Function evaluation
    template<typename T> inline static void fcn(const T& x, const T& y, T& f){ UnaryOperation<I>::fcn(x,f);}
    
    /// Partial derivatives - binary function
    template<typename T> inline static void der(const T& x, const T& y, const T& f, T* d){ UnaryOperation<I>::der(x,f,d); d[1]=0; }
};

template<int I>
class BinaryOperationE{
  public:
    /// Function evaluation
    template<typename T> inline static T fcn(const T& x, const T& y){ 
      T ret;
      BinaryOperation<I>::fcn(x,y,ret);
      return ret;
    }
};

template<int I>
class AddBinaryOperation{
  public:
    template<typename T> inline static void fcn(const T& x, const T& y, T& f){ f+= BinaryOperationE<I>::fcn(x,y);}
};

template<int I>
class SubBinaryOperation{
  public:
    template<typename T> inline static void fcn(const T& x, const T& y, T& f){ f-= BinaryOperationE<I>::fcn(x,y);}
};

template<int I>
class MulBinaryOperation{
  public:
    template<typename T> inline static void fcn(const T& x, const T& y, T& f){ f*= BinaryOperationE<I>::fcn(x,y);}
};

template<int I>
class DivBinaryOperation{
  public:
    template<typename T> inline static void fcn(const T& x, const T& y, T& f){ f/= BinaryOperationE<I>::fcn(x,y);}
};

/// Enum for quick access to any node
enum Operation{
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

/// Addition
template<>
class BinaryOperation<ADD>{
  public:
    template<typename T> inline static void fcn(const T& x, const T& y, T& f){ f = x+y;}
    template<typename T> inline static void der(const T& x, const T& y, const T& f, T* d){ d[0]=d[1]=1;}
};

/// Subtraction
template<>
class BinaryOperation<SUB>{
  public:
    template<typename T> inline static void fcn(const T& x, const T& y, T& f){ f = x-y;}
    template<typename T> inline static void der(const T& x, const T& y, const T& f, T* d){ d[0]=1; d[1]=-1;}
};

/// Multiplication
template<>
class BinaryOperation<MUL>{
  public:
    template<typename T> inline static void fcn(const T& x, const T& y, T& f){ f = x*y;}
    template<typename T> inline static void der(const T& x, const T& y, const T& f, T* d){ d[0]=y; d[1]=x;}
};

/// Division
template<>
class BinaryOperation<DIV>{
  public:
    template<typename T> inline static void fcn(const T& x, const T& y, T& f){ f = x/y;}
    template<typename T> inline static void der(const T& x, const T& y, const T& f, T* d){ d[0]=1/y; d[1]=-f/y;}
};

/// Negation
template<>
class UnaryOperation<NEG>{
  public:
    template<typename T> inline static void fcn(const T& x, T& f){ f = -x;}
    template<typename T> inline static void der(const T& x, const T& f, T* d){ d[0]=-1;}
};

/// Natural exponent
template<>
class UnaryOperation<EXP>{
  public:
    template<typename T> inline static void fcn(const T& x, T& f){ f = exp(x);}
    template<typename T> inline static void der(const T& x, const T& f, T* d){ d[0]=f;}
};

/// Natural logarithm
template<>
class UnaryOperation<LOG>{
  public:
    template<typename T> inline static void fcn(const T& x, T& f){ f = log(x);}
    template<typename T> inline static void der(const T& x, const T& f, T* d){ d[0]=1/x;}
};

/// Power, defined only for x>=0
template<>
class BinaryOperation<POW>{
  public:
    template<typename T> inline static void fcn(const T& x, const T& y, T& f){ f = pow(x,y);}
    // See issue #104 why d[0] is no longer y*f/x
    template<typename T> inline static void der(const T& x, const T& y, const T& f, T* d){ d[0]=y*pow(x,y-1); d[1]=log(x)*f;}
};

/// Power, defined only for y constant
template<>
class BinaryOperation<CONSTPOW>{
  public:
    template<typename T> inline static void fcn(const T& x, const T& y, T& f){ f = pow(x,y);}
    template<typename T> inline static void der(const T& x, const T& y, const T& f, T* d){ d[0]=y*pow(x,y-1); d[1]=0;}
};

/// Square root
template<>
class UnaryOperation<SQRT>{
  public:
    template<typename T> inline static void fcn(const T& x, T& f){ f = sqrt(x);}
    template<typename T> inline static void der(const T& x, const T& f, T* d){ d[0]=1/(timesTwo(f));}
};

/// Sine
template<>
class UnaryOperation<SIN>{
  public:
    template<typename T> inline static void fcn(const T& x, T& f){ f = sin(x);}
    template<typename T> inline static void der(const T& x, const T& f, T* d){ d[0]=cos(x);}
};

/// Cosine
template<>
class UnaryOperation<COS>{
  public:
    template<typename T> inline static void fcn(const T& x, T& f){ f = cos(x);}
    template<typename T> inline static void der(const T& x, const T& f, T* d){ d[0]=-sin(x);}
};

/// Tangent
template<>
class UnaryOperation<TAN>{
  public:
    template<typename T> inline static void fcn(const T& x, T& f){ f = tan(x);}
    template<typename T> inline static void der(const T& x, const T& f, T* d){ d[0] = 1/square(cos(x));}
};

/// Arcus sine
template<>
class UnaryOperation<ASIN>{
  public:
    template<typename T> inline static void fcn(const T& x, T& f){ f = asin(x);}
    template<typename T> inline static void der(const T& x, const T& f, T* d){ d[0]=1/sqrt(1-x*x);}
};

/// Arcus cosine
template<>
class UnaryOperation<ACOS>{
  public:
    template<typename T> inline static void fcn(const T& x, T& f){ f = acos(x);}
    template<typename T> inline static void der(const T& x, const T& f, T* d){ d[0]=-1/sqrt(1-x*x);}
};

/// Arcus tangent
template<>
class UnaryOperation<ATAN>{
  public:
    template<typename T> inline static void fcn(const T& x, T& f){ f = atan(x);}
    template<typename T> inline static void der(const T& x, const T& f, T* d){ d[0] = 1/(1+x*x);}
};

/// Step function
template<>
class UnaryOperation<STEP>{
  public:
    template<typename T> inline static void fcn(const T& x, T& f){ f = x >= T(0.);}
    template<typename T> inline static void der(const T& x, const T& f, T* d){ d[0] = 0;}
};

/// Floor function
template<>
class UnaryOperation<FLOOR>{
  public:
    template<typename T> inline static void fcn(const T& x, T& f){ f = floor(x);}
    template<typename T> inline static void der(const T& x, const T& f, T* d){ d[0] = 0;}
};

/// Ceil function
template<>
class UnaryOperation<CEIL>{
  public:
    template<typename T> inline static void fcn(const T& x, T& f){ f = ceil(x);}
    template<typename T> inline static void der(const T& x, const T& f, T* d){ d[0] = 0;}
};

/// Equality
template<>
class BinaryOperation<EQUALITY>{
  public:
    template<typename T> inline static void fcn(const T& x, const T& y, T& f){ f = x==y;}
    template<typename T> inline static void der(const T& x, const T& y, const T& f, T* d){ d[0]=d[1]=0;}
};

/// Error function
template<>
class UnaryOperation<ERF>{
  public:
    template<typename T> inline static void fcn(const T& x, T& f){ f = erf(x);}
    template<typename T> inline static void der(const T& x, const T& f, T* d){ d[0] = (2/sqrt(M_PI))*exp(-x*x);}
};

/// Absolute value
template<>
class UnaryOperation<FABS>{
  public:
    template<typename T> inline static void fcn(const T& x, T& f){ f = fabs(x);}
    template<typename T> inline static void der(const T& x, const T& f, T* d){ d[0]=sign(x);}
};

/// Sign
template<>
class UnaryOperation<SIGN>{
  public:
    template<typename T> inline static void fcn(const T& x, T& f){ f = sign(x);}
    template<typename T> inline static void der(const T& x, const T& f, T* d){ d[0]=0;}
};

/// Minimum
template<>
class BinaryOperation<FMIN>{
  public:
    template<typename T> inline static void fcn(const T& x, const T& y, T& f){ f = fmin(x,y);}
    template<typename T> inline static void der(const T& x, const T& y, const T& f, T* d){ d[0]=x<=y; d[1]=!d[0];}
};

/// Maximum
template<>
class BinaryOperation<FMAX>{
  public:
    template<typename T> inline static void fcn(const T& x, const T& y, T& f){ f = fmax(x,y);}
    template<typename T> inline static void der(const T& x, const T& y, const T& f, T* d){ d[0]=x>=y; d[1]=!d[0];}
};

/// Elementwise inverse
template<>
class UnaryOperation<INV>{
  public:
    template<typename T> inline static void fcn(const T& x, T& f){ f = 1./x;}
    template<typename T> inline static void der(const T& x, const T& f, T* d){ d[0] = -f*f; }
};

/// Hyperbolic sine
template<>
class UnaryOperation<SINH>{
  public:
    template<typename T> inline static void fcn(const T& x, T& f){ f = sinh(x);}
    template<typename T> inline static void der(const T& x, const T& f, T* d){ d[0] = cosh(x); }
};

/// Hyperbolic cosine
template<>
class UnaryOperation<COSH>{
  public:
    template<typename T> inline static void fcn(const T& x, T& f){ f = cosh(x);}
    template<typename T> inline static void der(const T& x, const T& f, T* d){ d[0] = -sinh(x); }
};

/// Hyperbolic tangent
template<>
class UnaryOperation<TANH>{
  public:
    template<typename T> inline static void fcn(const T& x, T& f){ f = tanh(x);}
    template<typename T> inline static void der(const T& x, const T& f, T* d){ d[0] = 1-f*f; }
};

/// Inverse of error function
template<>
class UnaryOperation<ERFINV>{
  public:
    template<typename T> inline static void fcn(const T& x, T& f){ f = erfinv(x);}
    template<typename T> inline static void der(const T& x, const T& f, T* d){ d[0] = (sqrt(M_PI)/2)*exp(f*f); }
};

/// Identity operator with the side effect of printing
template<>
class BinaryOperation<PRINTME>{
  public:
    template<typename T> inline static void fcn(const T& x, const T& y, T& f){f = printme(x,y); }
    template<typename T> inline static void der(const T& x, const T& y, const T& f, T* d){ d[0]=1; d[1]=0;}
};

template<>
inline void BinaryOperation<PRINTME>::fcn<double>(const double& x, const double& y, double& f) {
  f = x; 
  std::cout << "|> " << y << " : " << x << std::endl;
}

} // namespace CasADi

#endif //CASADI_CALCULUS_HPP
