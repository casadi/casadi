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

#ifndef CASADI_MATH_HPP
#define CASADI_MATH_HPP

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include "casadi_exception.hpp"
#include "pre_c99_support.hpp"

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
    /// Print
    static void print(std::ostream &stream, const std::string& x);

    /// Does evaluating the function with a zero give a zero
    static bool f0_is_zero(){return false;}

    /// Function evaluation
    template<typename T> static void fcn(const T& x, T& f);
    
    /// Partial derivatives
    template<typename T> static void der(const T& x, const T& f, T* d);
};

template<int I>
class BinaryOperation{
  public:
    /// Print
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ UnaryOperation<I>::print(stream,x); }

    /// Does evaluating the function with a zero give a zero
    static bool f00_is_zero(){return UnaryOperation<I>::f0_is_zero();}
    static bool f0x_is_zero(){return UnaryOperation<I>::f0_is_zero();}
    static bool fx0_is_zero(){return false;}

    /// Function evaluation
    template<typename T> static void fcn(const T& x, const T& y, T& f){ UnaryOperation<I>::fcn(x,f);}
    
    /// Partial derivatives - binary function
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ UnaryOperation<I>::der(x,f,d); d[1]=0; }
};

template<int I>
class BinaryOperationE{
  public:
    /// Function evaluation
    template<typename T> static T fcn(const T& x, const T& y){ 
      T ret;
      BinaryOperation<I>::fcn(x,y,ret);
      return ret;
    }
};

/// Enum for quick access to any node
enum Operation{
  ADD,  SUB,  MUL,  DIV,
  NEG,  EXP,  LOG,  POW, CONSTPOW,
  SQRT,  SIN,  COS,  TAN,  
  ASIN,  ACOS,  ATAN,  
  STEP,  
  FLOOR,  CEIL,  EQUALITY,  
  ERF,  FMIN,  FMAX,
  NUM_BUILT_IN_OPS
};

/// Addition
template<>
class BinaryOperation<ADD>{
  public:
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ stream << "(" << x << "+" << y << ")"; }
    static bool f00_is_zero(){return true;}
    static bool f0x_is_zero(){return false;}
    static bool fx0_is_zero(){return false;}
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = x+y;}
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=d[1]=1;}
};

/// Subtraction
template<>
class BinaryOperation<SUB>{
  public:
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ stream << "(" << x << "-" << y << ")"; }
    static bool f00_is_zero(){return true;}
    static bool f0x_is_zero(){return false;}
    static bool fx0_is_zero(){return false;}
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = x-y;}
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=1; d[1]=-1;}
};

/// Multiplication
template<>
class BinaryOperation<MUL>{
  public:
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ stream << "(" << x << "*" << y << ")"; }
    static bool f00_is_zero(){return true;}
    static bool f0x_is_zero(){return true;}
    static bool fx0_is_zero(){return true;}
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = x*y;}
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=y; d[1]=x;}
};

/// Division
template<>
class BinaryOperation<DIV>{
  public:
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ stream << "(" << x << "/" << y << ")"; }
    static bool f00_is_zero(){return false;}
    static bool f0x_is_zero(){return true;}
    static bool fx0_is_zero(){return false;}
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = x/y;}
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=1/y; d[1]=-f/y;}
};

/// Negation
template<>
class UnaryOperation<NEG>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "(-" << x << ")"; }
    static bool f0_is_zero(){return true;}
    template<typename T> static void fcn(const T& x, T& f){ f = -x;}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0]=-1;}
};

/// Natural exponent
template<>
class UnaryOperation<EXP>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "exp(" << x << ")"; }
    static bool f0_is_zero(){return false;}
    template<typename T> static void fcn(const T& x, T& f){ f = std::exp(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0]=f;}
};

/// Natural logarithm
template<>
class UnaryOperation<LOG>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "log(" << x << ")"; }
    static bool f0_is_zero(){return false;}
    template<typename T> static void fcn(const T& x, T& f){ f = std::log(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0]=1/x;}
};

/// Power, defined only for x>=0
template<>
class BinaryOperation<POW>{
  public:
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ stream << "pow(" << x << "," << y << ")"; }
    static bool f00_is_zero(){return false;}
    static bool f0x_is_zero(){return false;}
    static bool fx0_is_zero(){return false;}
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = std::pow(x,y);}
    // See issue #104 why d[0] is no longer y*f/x
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=y*std::pow(x,y-1); d[1]=std::log(x)*f;}
};

/// Power, defined only for y constant
template<>
class BinaryOperation<CONSTPOW>{
  public:
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ stream << "constpow(" << x << "," << y << ")"; }
    static bool f00_is_zero(){return false;}
    static bool f0x_is_zero(){return false;}
    static bool fx0_is_zero(){return false;}
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = std::pow(x,y);}
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=y*std::pow(x,y-1); d[1]=0;}
};

/// Square root
template<>
class UnaryOperation<SQRT>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "sqrt(" << x << ")"; }
    static bool f0_is_zero(){return true;}
    template<typename T> static void fcn(const T& x, T& f){ f = std::sqrt(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0]=1/(timesTwo(f));}
};

/// Sine
template<>
class UnaryOperation<SIN>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "sin(" << x << ")"; }
    static bool f0_is_zero(){return true;}
    template<typename T> static void fcn(const T& x, T& f){ f = std::sin(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0]=std::cos(x);}
};

/// Cosine
template<>
class UnaryOperation<COS>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "cos(" << x << ")"; }
    static bool f0_is_zero(){return false;}
    template<typename T> static void fcn(const T& x, T& f){ f = std::cos(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0]=-std::sin(x);}
};

/// Tangent
template<>
class UnaryOperation<TAN>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "tan(" << x << ")"; }
    static bool f0_is_zero(){return true;}
    template<typename T> static void fcn(const T& x, T& f){ f = std::tan(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0] = 1/square(std::cos(x));}
};

/// Arcus sine
template<>
class UnaryOperation<ASIN>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "asin(" << x << ")"; }
    static bool f0_is_zero(){return true;}
    template<typename T> static void fcn(const T& x, T& f){ f = std::asin(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0]=1/std::sqrt(1-x*x);}
};

/// Arcus cosine
template<>
class UnaryOperation<ACOS>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "acos(" << x << ")"; }
    static bool f0_is_zero(){return false;}
    template<typename T> static void fcn(const T& x, T& f){ f = std::acos(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0]=-1/std::sqrt(1-x*x);}
};

/// Arcus tangent
template<>
class UnaryOperation<ATAN>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "atan(" << x << ")"; }
    static bool f0_is_zero(){return true;}
    template<typename T> static void fcn(const T& x, T& f){ f = std::atan(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0] = 1/(1+x*x);}
};

/// Step function
template<>
class UnaryOperation<STEP>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "(" << x << ">=0)"; }
    static bool f0_is_zero(){return false;}
    template<typename T> static void fcn(const T& x, T& f){ f = x >= T(0.);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0] = 0;}
};

/// Floor function
template<>
class UnaryOperation<FLOOR>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "floor(" << x << ")"; }
    static bool f0_is_zero(){return true;}
    template<typename T> static void fcn(const T& x, T& f){ f = std::floor(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0] = 0;}
};

/// Ceil function
template<>
class UnaryOperation<CEIL>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "ceil(" << x << ")"; }
    static bool f0_is_zero(){return true;}
    template<typename T> static void fcn(const T& x, T& f){ f = std::ceil(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0] = 0;}
};

/// Equality
template<>
class BinaryOperation<EQUALITY>{
  public:
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ stream << "(" << x << "==" << y << ")"; }
    static bool f00_is_zero(){return false;}
    static bool f0x_is_zero(){return false;}
    static bool fx0_is_zero(){return false;}
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = x==y;}
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=d[1]=0;}
};

/// Minimum
template<>
class BinaryOperation<FMIN>{
  public:
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ stream << "fmin(" << x << "," << y << ")"; }
    static bool f00_is_zero(){return true;}
    static bool f0x_is_zero(){return false;}
    static bool fx0_is_zero(){return false;}
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = fmin(x,y);}
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=x<=y; d[1]=!d[0];}
};

/// Maximum
template<>
class BinaryOperation<FMAX>{
  public:
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ stream << "fmax(" << x << "," << y << ")"; }
    static bool f00_is_zero(){return true;}
    static bool f0x_is_zero(){return false;}
    static bool fx0_is_zero(){return false;}
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = fmax(x,y);}
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=x>=y; d[1]=!d[0];}
};

/// Error function
template<>
class UnaryOperation<ERF>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "erf(" << x << ")"; }
    static bool f0_is_zero(){return true;}
    template<typename T> static void fcn(const T& x, T& f){ f = erf(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0] = (2/std::sqrt(M_PI))*std::exp(-x*x);}
};

/// Easy access to all the functions for a particular type
template<typename T>
class casadi_math{
  public:

    /** \brief Printing operation typedef */
    typedef void (*printFunT)(std::ostream &stream, const std::string& x, const std::string& y);

    /** \brief Function typedef */
    typedef void (*funT)(const T&, const T&, T&);
    typedef T (*funTE)(const T&, const T&);
        
    /** \brief Derivative typedef */
    typedef void (*derT)(const T& x, const T& y, const T& f, T* d);

    /** \brief Vector of printing functions */
    static std::vector<printFunT> print;

    /** \brief Vector of function pointers to all the built in functions */
    static std::vector<funT> fun;
    static std::vector<funTE> funE;
    
    /** \brief Vector of function derivative pointers to all the built in functions */
    static std::vector<derT> der;
    
    /** \brief Vector of booleans indicating which functions are zero when evaluating with both arguments zero */
    static std::vector<bool> f00_is_zero;
    
    /** \brief Vector of booleans indicating which functions are zero when evaluating with the first arguments zero */
    static std::vector<bool> f0x_is_zero;
    
    /** \brief Vector of booleans indicating which functions are zero when evaluating with the second arguments zero */
    static std::vector<bool> fx0_is_zero;
    
  protected:
    
    /** \brief Create print */
    static std::vector<printFunT> getPrintFun();

    /** \brief Create fun */
    static std::vector<funT> getFun();
    static std::vector<funTE> getFunE();
  
    /** \brief Create der */
    static std::vector<derT> getDer();
    
    /** \brief Create zero indicator vectors */
    static std::vector<bool> getF00_is_zero();
    static std::vector<bool> getF0x_is_zero();
    static std::vector<bool> getFx0_is_zero();
    
    
};


// Template implementations

template<typename T>
std::vector<typename casadi_math<T>::printFunT> casadi_math<T>::print = casadi_math<T>::getPrintFun();

template<typename T>
std::vector<typename casadi_math<T>::funT> casadi_math<T>::fun = casadi_math<T>::getFun();

template<typename T>
std::vector<typename casadi_math<T>::funTE> casadi_math<T>::funE = casadi_math<T>::getFunE();

template<typename T>
std::vector<typename casadi_math<T>::derT> casadi_math<T>::der = casadi_math<T>::getDer();

template<typename T>
std::vector<bool> casadi_math<T>::f00_is_zero = casadi_math<T>::getF00_is_zero();

template<typename T>
std::vector<bool> casadi_math<T>::f0x_is_zero = casadi_math<T>::getF0x_is_zero();

template<typename T>
std::vector<bool> casadi_math<T>::fx0_is_zero = casadi_math<T>::getFx0_is_zero();

template<typename T>
std::vector<typename casadi_math<T>::printFunT> casadi_math<T>::getPrintFun(){
  // Create return object
  std::vector<typename casadi_math<T>::printFunT> ret(NUM_BUILT_IN_OPS,0);
  
  // Specify operations
  ret[ADD] = BinaryOperation<ADD>::print;
  ret[SUB] = BinaryOperation<SUB>::print;
  ret[MUL] = BinaryOperation<MUL>::print;
  ret[DIV] = BinaryOperation<DIV>::print;
  
  ret[NEG] = BinaryOperation<NEG>::print;
  ret[EXP] = BinaryOperation<EXP>::print;
  ret[LOG] = BinaryOperation<LOG>::print;
  ret[POW] = BinaryOperation<POW>::print;
  ret[CONSTPOW] = BinaryOperation<CONSTPOW>::print;

  ret[SQRT] = BinaryOperation<SQRT>::print;
  ret[SIN] = BinaryOperation<SIN>::print;
  ret[COS] = BinaryOperation<COS>::print;
  ret[TAN] = BinaryOperation<TAN>::print;

  ret[ASIN] = BinaryOperation<ASIN>::print;
  ret[ACOS] = BinaryOperation<ACOS>::print;
  ret[ATAN] = BinaryOperation<ATAN>::print;

  ret[STEP] = BinaryOperation<STEP>::print;
  ret[FLOOR] = BinaryOperation<FLOOR>::print;
  ret[CEIL] = BinaryOperation<CEIL>::print;

  ret[EQUALITY] = BinaryOperation<EQUALITY>::print;

  ret[ERF] = BinaryOperation<ERF>::print;
  ret[FMIN] = BinaryOperation<FMIN>::print;
  ret[FMAX] = BinaryOperation<FMAX>::print;

  // Make sure that all functions were specified
  for(int i=0; i<ret.size(); ++i){
    casadi_assert(ret[i]!=0);
  }
  
  return ret;
  
}

template<typename T>
std::vector<typename casadi_math<T>::funT> casadi_math<T>::getFun(){
  
  // Create return object
  std::vector<typename casadi_math<T>::funT> ret(NUM_BUILT_IN_OPS,0);
  
#ifdef _MSC_VER
  // Specify operations
  ret[ADD] = reinterpret_cast<funT>(BinaryOperation<ADD>::fcn);
  ret[SUB] = reinterpret_cast<funT>(BinaryOperation<SUB>::fcn);
  ret[MUL] = reinterpret_cast<funT>(BinaryOperation<MUL>::fcn);
  ret[DIV] = reinterpret_cast<funT>(BinaryOperation<DIV>::fcn);
    
  ret[NEG] = reinterpret_cast<funT>(BinaryOperation<NEG>::fcn);
  ret[EXP] = reinterpret_cast<funT>(BinaryOperation<EXP>::fcn);
  ret[LOG] = reinterpret_cast<funT>(BinaryOperation<LOG>::fcn);
  ret[POW] = reinterpret_cast<funT>(BinaryOperation<POW>::fcn);
  ret[CONSTPOW] = reinterpret_cast<funT>(BinaryOperation<CONSTPOW>::fcn);

  ret[SQRT] = reinterpret_cast<funT>(BinaryOperation<SQRT>::fcn);
  ret[SIN] = reinterpret_cast<funT>(BinaryOperation<SIN>::fcn);
  ret[COS] = reinterpret_cast<funT>(BinaryOperation<COS>::fcn);
  ret[TAN] = reinterpret_cast<funT>(BinaryOperation<TAN>::fcn);

  ret[ASIN] = reinterpret_cast<funT>(BinaryOperation<ASIN>::fcn);
  ret[ACOS] = reinterpret_cast<funT>(BinaryOperation<ACOS>::fcn);
  ret[ATAN] = reinterpret_cast<funT>(BinaryOperation<ATAN>::fcn);

  ret[STEP] = reinterpret_cast<funT>(BinaryOperation<STEP>::fcn);
  ret[FLOOR] = reinterpret_cast<funT>(BinaryOperation<FLOOR>::fcn);
  ret[CEIL] = reinterpret_cast<funT>(BinaryOperation<CEIL>::fcn);

  ret[EQUALITY] = reinterpret_cast<funT>(BinaryOperation<EQUALITY>::fcn);

  ret[ERF] = reinterpret_cast<funT>(BinaryOperation<ERF>::fcn);
  ret[FMIN] = reinterpret_cast<funT>(BinaryOperation<FMIN>::fcn);
  ret[FMAX] = reinterpret_cast<funT>(BinaryOperation<FMAX>::fcn);
#else // _MSC_VER
  
  // Specify operations
  ret[ADD] = BinaryOperation<ADD>::fcn;
  ret[SUB] = BinaryOperation<SUB>::fcn;
  ret[MUL] = BinaryOperation<MUL>::fcn;
  ret[DIV] = BinaryOperation<DIV>::fcn;
    
  ret[NEG] = BinaryOperation<NEG>::fcn;
  ret[EXP] = BinaryOperation<EXP>::fcn;
  ret[LOG] = BinaryOperation<LOG>::fcn;
  ret[POW] = BinaryOperation<POW>::fcn;
  ret[CONSTPOW] = BinaryOperation<CONSTPOW>::fcn;

  ret[SQRT] = BinaryOperation<SQRT>::fcn;
  ret[SIN] = BinaryOperation<SIN>::fcn;
  ret[COS] = BinaryOperation<COS>::fcn;
  ret[TAN] = BinaryOperation<TAN>::fcn;

  ret[ASIN] = BinaryOperation<ASIN>::fcn;
  ret[ACOS] = BinaryOperation<ACOS>::fcn;
  ret[ATAN] = BinaryOperation<ATAN>::fcn;

  ret[STEP] = BinaryOperation<STEP>::fcn;
  ret[FLOOR] = BinaryOperation<FLOOR>::fcn;
  ret[CEIL] = BinaryOperation<CEIL>::fcn;

  ret[EQUALITY] = BinaryOperation<EQUALITY>::fcn;

  ret[ERF] = BinaryOperation<ERF>::fcn;
  ret[FMIN] = BinaryOperation<FMIN>::fcn;
  ret[FMAX] = BinaryOperation<FMAX>::fcn;
  
#endif // _MSC_VER

  // Make sure that all functions were specified
  for(int i=0; i<ret.size(); ++i){
    casadi_assert(ret[i]!=0);
  }
  
  return ret;
}

template<typename T>
std::vector<typename casadi_math<T>::funTE> casadi_math<T>::getFunE(){
  
  // Create return object
  std::vector<typename casadi_math<T>::funTE> ret(NUM_BUILT_IN_OPS,0);
  
#ifdef _MSC_VER
  // Specify operations
  ret[ADD] = reinterpret_cast<funT>(BinaryOperationE<ADD>::fcn);
  ret[SUB] = reinterpret_cast<funT>(BinaryOperationE<SUB>::fcn);
  ret[MUL] = reinterpret_cast<funT>(BinaryOperationE<MUL>::fcn);
  ret[DIV] = reinterpret_cast<funT>(BinaryOperationE<DIV>::fcn);
    
  ret[NEG] = reinterpret_cast<funT>(BinaryOperationE<NEG>::fcn);
  ret[EXP] = reinterpret_cast<funT>(BinaryOperationE<EXP>::fcn);
  ret[LOG] = reinterpret_cast<funT>(BinaryOperationE<LOG>::fcn);
  ret[POW] = reinterpret_cast<funT>(BinaryOperationE<POW>::fcn);
  ret[CONSTPOW] = reinterpret_cast<funT>(BinaryOperationE<CONSTPOW>::fcn);

  ret[SQRT] = reinterpret_cast<funT>(BinaryOperationE<SQRT>::fcn);
  ret[SIN] = reinterpret_cast<funT>(BinaryOperationE<SIN>::fcn);
  ret[COS] = reinterpret_cast<funT>(BinaryOperationE<COS>::fcn);
  ret[TAN] = reinterpret_cast<funT>(BinaryOperationE<TAN>::fcn);

  ret[ASIN] = reinterpret_cast<funT>(BinaryOperationE<ASIN>::fcn);
  ret[ACOS] = reinterpret_cast<funT>(BinaryOperationE<ACOS>::fcn);
  ret[ATAN] = reinterpret_cast<funT>(BinaryOperationE<ATAN>::fcn);

  ret[STEP] = reinterpret_cast<funT>(BinaryOperationE<STEP>::fcn);
  ret[FLOOR] = reinterpret_cast<funT>(BinaryOperationE<FLOOR>::fcn);
  ret[CEIL] = reinterpret_cast<funT>(BinaryOperationE<CEIL>::fcn);

  ret[EQUALITY] = reinterpret_cast<funT>(BinaryOperationE<EQUALITY>::fcn);

  ret[ERF] = reinterpret_cast<funT>(BinaryOperationE<ERF>::fcn);
  ret[FMIN] = reinterpret_cast<funT>(BinaryOperationE<FMIN>::fcn);
  ret[FMAX] = reinterpret_cast<funT>(BinaryOperationE<FMAX>::fcn);
#else // _MSC_VER
  
  // Specify operations
  ret[ADD] = BinaryOperationE<ADD>::fcn;
  ret[SUB] = BinaryOperationE<SUB>::fcn;
  ret[MUL] = BinaryOperationE<MUL>::fcn;
  ret[DIV] = BinaryOperationE<DIV>::fcn;
    
  ret[NEG] = BinaryOperationE<NEG>::fcn;
  ret[EXP] = BinaryOperationE<EXP>::fcn;
  ret[LOG] = BinaryOperationE<LOG>::fcn;
  ret[POW] = BinaryOperationE<POW>::fcn;
  ret[CONSTPOW] = BinaryOperationE<CONSTPOW>::fcn;

  ret[SQRT] = BinaryOperationE<SQRT>::fcn;
  ret[SIN] = BinaryOperationE<SIN>::fcn;
  ret[COS] = BinaryOperationE<COS>::fcn;
  ret[TAN] = BinaryOperationE<TAN>::fcn;

  ret[ASIN] = BinaryOperationE<ASIN>::fcn;
  ret[ACOS] = BinaryOperationE<ACOS>::fcn;
  ret[ATAN] = BinaryOperationE<ATAN>::fcn;

  ret[STEP] = BinaryOperationE<STEP>::fcn;
  ret[FLOOR] = BinaryOperationE<FLOOR>::fcn;
  ret[CEIL] = BinaryOperationE<CEIL>::fcn;

  ret[EQUALITY] = BinaryOperationE<EQUALITY>::fcn;

  ret[ERF] = BinaryOperationE<ERF>::fcn;
  ret[FMIN] = BinaryOperationE<FMIN>::fcn;
  ret[FMAX] = BinaryOperationE<FMAX>::fcn;
  
#endif // _MSC_VER

  // Make sure that all functions were specified
  for(int i=0; i<ret.size(); ++i){
    casadi_assert(ret[i]!=0);
  }
  
  return ret;
}

template<typename T>
std::vector<typename casadi_math<T>::derT> casadi_math<T>::getDer(){
  // Create return object
  std::vector<typename casadi_math<T>::derT> ret(NUM_BUILT_IN_OPS,0);
  
  // Specify operations
#ifdef _MSC_VER
  ret[ADD] = reinterpret_cast<derT>(BinaryOperation<ADD>::der);
  ret[SUB] = reinterpret_cast<derT>(BinaryOperation<SUB>::der);
  ret[MUL] = reinterpret_cast<derT>(BinaryOperation<MUL>::der);
  ret[DIV] = reinterpret_cast<derT>(BinaryOperation<DIV>::der);
    
  ret[NEG] = reinterpret_cast<derT>(BinaryOperation<NEG>::der);
  ret[EXP] = reinterpret_cast<derT>(BinaryOperation<EXP>::der);
  ret[LOG] = reinterpret_cast<derT>(BinaryOperation<LOG>::der);
  ret[POW] = reinterpret_cast<derT>(BinaryOperation<POW>::der);
  ret[CONSTPOW] = reinterpret_cast<derT>(BinaryOperation<CONSTPOW>::der);

  ret[SQRT] = reinterpret_cast<derT>(BinaryOperation<SQRT>::der);
  ret[SIN] = reinterpret_cast<derT>(BinaryOperation<SIN>::der);
  ret[COS] = reinterpret_cast<derT>(BinaryOperation<COS>::der);
  ret[TAN] = reinterpret_cast<derT>(BinaryOperation<TAN>::der);

  ret[ASIN] = reinterpret_cast<derT>(BinaryOperation<ASIN>::der);
  ret[ACOS] = reinterpret_cast<derT>(BinaryOperation<ACOS>::der);
  ret[ATAN] = reinterpret_cast<derT>(BinaryOperation<ATAN>::der);

  ret[STEP] = reinterpret_cast<derT>(BinaryOperation<STEP>::der);
  ret[FLOOR] = reinterpret_cast<derT>(BinaryOperation<FLOOR>::der);
  ret[CEIL] = reinterpret_cast<derT>(BinaryOperation<CEIL>::der);

  ret[EQUALITY] = reinterpret_cast<derT>(BinaryOperation<EQUALITY>::der);

  ret[ERF] = reinterpret_cast<derT>(BinaryOperation<ERF>::der);
  ret[FMIN] = reinterpret_cast<derT>(BinaryOperation<FMIN>::der);
  ret[FMAX] = reinterpret_cast<derT>(BinaryOperation<FMAX>::der);

#else // _MSC_VER
  ret[ADD] = BinaryOperation<ADD>::der;
  ret[SUB] = BinaryOperation<SUB>::der;
  ret[MUL] = BinaryOperation<MUL>::der;
  ret[DIV] = BinaryOperation<DIV>::der;
    
  ret[NEG] = BinaryOperation<NEG>::der;
  ret[EXP] = BinaryOperation<EXP>::der;
  ret[LOG] = BinaryOperation<LOG>::der;
  ret[POW] = BinaryOperation<POW>::der;
  ret[CONSTPOW] = BinaryOperation<CONSTPOW>::der;

  ret[SQRT] = BinaryOperation<SQRT>::der;
  ret[SIN] = BinaryOperation<SIN>::der;
  ret[COS] = BinaryOperation<COS>::der;
  ret[TAN] = BinaryOperation<TAN>::der;

  ret[ASIN] = BinaryOperation<ASIN>::der;
  ret[ACOS] = BinaryOperation<ACOS>::der;
  ret[ATAN] = BinaryOperation<ATAN>::der;

  ret[STEP] = BinaryOperation<STEP>::der;
  ret[FLOOR] = BinaryOperation<FLOOR>::der;
  ret[CEIL] = BinaryOperation<CEIL>::der;

  ret[EQUALITY] = BinaryOperation<EQUALITY>::der;

  ret[ERF] = BinaryOperation<ERF>::der;
  ret[FMIN] = BinaryOperation<FMIN>::der;
  ret[FMAX] = BinaryOperation<FMAX>::der;
#endif // _MSC_VER

  // Make sure that all functions were specified
  for(int i=0; i<ret.size(); ++i){
    casadi_assert(ret[i]!=0);
  }
  
  return ret;
}

template<typename T>
std::vector<bool> casadi_math<T>::getF00_is_zero(){
   // Create return object
  std::vector<bool> ret(NUM_BUILT_IN_OPS,false);
  
  // Specify operations
  ret[ADD] = BinaryOperation<ADD>::f00_is_zero();
  ret[SUB] = BinaryOperation<SUB>::f00_is_zero();
  ret[MUL] = BinaryOperation<MUL>::f00_is_zero();
  ret[DIV] = BinaryOperation<DIV>::f00_is_zero();
  
  ret[NEG] = BinaryOperation<NEG>::f00_is_zero();
  ret[EXP] = BinaryOperation<EXP>::f00_is_zero();
  ret[LOG] = BinaryOperation<LOG>::f00_is_zero();
  ret[POW] = BinaryOperation<POW>::f00_is_zero();
  ret[CONSTPOW] = BinaryOperation<CONSTPOW>::f00_is_zero();

  ret[SQRT] = BinaryOperation<SQRT>::f00_is_zero();
  ret[SIN] = BinaryOperation<SIN>::f00_is_zero();
  ret[COS] = BinaryOperation<COS>::f00_is_zero();
  ret[TAN] = BinaryOperation<TAN>::f00_is_zero();

  ret[ASIN] = BinaryOperation<ASIN>::f00_is_zero();
  ret[ACOS] = BinaryOperation<ACOS>::f00_is_zero();
  ret[ATAN] = BinaryOperation<ATAN>::f00_is_zero();

  ret[STEP] = BinaryOperation<STEP>::f00_is_zero();
  ret[FLOOR] = BinaryOperation<FLOOR>::f00_is_zero();
  ret[CEIL] = BinaryOperation<CEIL>::f00_is_zero();

  ret[EQUALITY] = BinaryOperation<EQUALITY>::f00_is_zero();

  ret[ERF] = BinaryOperation<ERF>::f00_is_zero();
  ret[FMIN] = BinaryOperation<FMIN>::f00_is_zero();
  ret[FMAX] = BinaryOperation<FMAX>::f00_is_zero();

  return ret;
}

template<typename T>
std::vector<bool> casadi_math<T>::getF0x_is_zero(){
   // Create return object
  std::vector<bool> ret(NUM_BUILT_IN_OPS,false);
  
  // Specify operations
  ret[ADD] = BinaryOperation<ADD>::f0x_is_zero();
  ret[SUB] = BinaryOperation<SUB>::f0x_is_zero();
  ret[MUL] = BinaryOperation<MUL>::f0x_is_zero();
  ret[DIV] = BinaryOperation<DIV>::f0x_is_zero();
  
  ret[NEG] = BinaryOperation<NEG>::f0x_is_zero();
  ret[EXP] = BinaryOperation<EXP>::f0x_is_zero();
  ret[LOG] = BinaryOperation<LOG>::f0x_is_zero();
  ret[POW] = BinaryOperation<POW>::f0x_is_zero();
  ret[CONSTPOW] = BinaryOperation<CONSTPOW>::f0x_is_zero();

  ret[SQRT] = BinaryOperation<SQRT>::f0x_is_zero();
  ret[SIN] = BinaryOperation<SIN>::f0x_is_zero();
  ret[COS] = BinaryOperation<COS>::f0x_is_zero();
  ret[TAN] = BinaryOperation<TAN>::f0x_is_zero();

  ret[ASIN] = BinaryOperation<ASIN>::f0x_is_zero();
  ret[ACOS] = BinaryOperation<ACOS>::f0x_is_zero();
  ret[ATAN] = BinaryOperation<ATAN>::f0x_is_zero();

  ret[STEP] = BinaryOperation<STEP>::f0x_is_zero();
  ret[FLOOR] = BinaryOperation<FLOOR>::f0x_is_zero();
  ret[CEIL] = BinaryOperation<CEIL>::f0x_is_zero();

  ret[EQUALITY] = BinaryOperation<EQUALITY>::f0x_is_zero();

  ret[ERF] = BinaryOperation<ERF>::f0x_is_zero();
  ret[FMIN] = BinaryOperation<FMIN>::f0x_is_zero();
  ret[FMAX] = BinaryOperation<FMAX>::f0x_is_zero();

  return ret;
}

template<typename T>
std::vector<bool> casadi_math<T>::getFx0_is_zero(){
   // Create return object
  std::vector<bool> ret(NUM_BUILT_IN_OPS,false);
  
  // Specify operations
  ret[ADD] = BinaryOperation<ADD>::fx0_is_zero();
  ret[SUB] = BinaryOperation<SUB>::fx0_is_zero();
  ret[MUL] = BinaryOperation<MUL>::fx0_is_zero();
  ret[DIV] = BinaryOperation<DIV>::fx0_is_zero();
  
  ret[NEG] = BinaryOperation<NEG>::fx0_is_zero();
  ret[EXP] = BinaryOperation<EXP>::fx0_is_zero();
  ret[LOG] = BinaryOperation<LOG>::fx0_is_zero();
  ret[POW] = BinaryOperation<POW>::fx0_is_zero();
  ret[CONSTPOW] = BinaryOperation<CONSTPOW>::fx0_is_zero();

  ret[SQRT] = BinaryOperation<SQRT>::fx0_is_zero();
  ret[SIN] = BinaryOperation<SIN>::fx0_is_zero();
  ret[COS] = BinaryOperation<COS>::fx0_is_zero();
  ret[TAN] = BinaryOperation<TAN>::fx0_is_zero();

  ret[ASIN] = BinaryOperation<ASIN>::fx0_is_zero();
  ret[ACOS] = BinaryOperation<ACOS>::fx0_is_zero();
  ret[ATAN] = BinaryOperation<ATAN>::fx0_is_zero();

  ret[STEP] = BinaryOperation<STEP>::fx0_is_zero();
  ret[FLOOR] = BinaryOperation<FLOOR>::fx0_is_zero();
  ret[CEIL] = BinaryOperation<CEIL>::fx0_is_zero();

  ret[EQUALITY] = BinaryOperation<EQUALITY>::fx0_is_zero();

  ret[ERF] = BinaryOperation<ERF>::fx0_is_zero();
  ret[FMIN] = BinaryOperation<FMIN>::fx0_is_zero();
  ret[FMAX] = BinaryOperation<FMAX>::fx0_is_zero();

  return ret;
}


} // namespace CasADi

#endif //CASADI_MATH_HPP
