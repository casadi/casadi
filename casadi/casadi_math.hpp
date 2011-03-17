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

namespace CasADi{
  
template<int I>
class UnaryOperation{
  public:
    /// Print
    static void print(std::ostream &stream, const std::string& x);

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

    /// Function evaluation
    template<typename T> static void fcn(const T& x, const T& y, T& f){ UnaryOperation<I>::fcn(x,f);}
    
    /// Partial derivatives - binary function
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ UnaryOperation<I>::der(x,f,d); d[1]=0; }
};

/// Enum for quick access to any node
enum Operation{
  ADD,  SUB,  MUL,  DIV,
  NEG,  EXP,  LOG,  POW,  
  SQRT,  SIN,  COS,  TAN,  
  ASIN,  ACOS,  ATAN,  
  STEP,  
  FLOOR,  CEIL,  
  EQUALITY,  ERF,  FMIN,  FMAX,
  NUM_BUILT_IN_OPS
};

/// Addition
template<>
class BinaryOperation<ADD>{
  public:
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ stream << "(" << x << "+" << y << ")"; }
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = x+y;}
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=d[1]=1;}
};

/// Subtraction
template<>
class BinaryOperation<SUB>{
  public:
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ stream << "(" << x << "-" << y << ")"; }
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = x-y;}
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=1; d[1]=-1;}
};

/// Multiplication
template<>
class BinaryOperation<MUL>{
  public:
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ stream << "(" << x << "*" << y << ")"; }
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = x*y;}
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=y; d[1]=x;}
};

/// Division
template<>
class BinaryOperation<DIV>{
  public:
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ stream << "(" << x << "/" << y << ")"; }
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = x/y;}
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=1/y; d[1]=-f/y;}
};

/// Negation
template<>
class UnaryOperation<NEG>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "(-" << x << ")"; }
    template<typename T> static void fcn(const T& x, T& f){ f = -x;}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0]=-1;}
};

/// Natural exponent
template<>
class UnaryOperation<EXP>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "exp(" << x << ")"; }
    template<typename T> static void fcn(const T& x, T& f){ f = std::exp(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0]=f;}
};

/// Natural logarithm
template<>
class UnaryOperation<LOG>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "log(" << x << ")"; }
    template<typename T> static void fcn(const T& x, T& f){ f = std::log(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0]=1/x;}
};

/// Power
template<>
class BinaryOperation<POW>{
  public:
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ stream << "pow(" << x << "," << y << ")"; }
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = std::pow(x,y);}
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=y*f/x; d[1]=std::log(x)*f;}
};

/// Square root
template<>
class UnaryOperation<SQRT>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "sqrt(" << x << ")"; }
    template<typename T> static void fcn(const T& x, T& f){ f = std::sqrt(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0]=1/(2*f);}
};

/// Sine
template<>
class UnaryOperation<SIN>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "sin(" << x << ")"; }
    template<typename T> static void fcn(const T& x, T& f){ f = std::sin(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0]=std::cos(x);}
};

/// Cosine
template<>
class UnaryOperation<COS>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "cos(" << x << ")"; }
    template<typename T> static void fcn(const T& x, T& f){ f = std::cos(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0]=-std::sin(x);}
};

/// Tangent
template<>
class UnaryOperation<TAN>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "tan(" << x << ")"; }
    template<typename T> static void fcn(const T& x, T& f){ f = std::tan(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ T cosx = std::cos(x); d[0] = 1/(cosx*cosx);}
};

/// Arcus sine
template<>
class UnaryOperation<ASIN>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "asin(" << x << ")"; }
    template<typename T> static void fcn(const T& x, T& f){ f = std::asin(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0]=1/std::sqrt(1-x*x);}
};

/// Arcus cosine
template<>
class UnaryOperation<ACOS>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "acos(" << x << ")"; }
    template<typename T> static void fcn(const T& x, T& f){ f = std::acos(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0]=-1/std::sqrt(1-x*x);}
};

/// Arcus tangent
template<>
class UnaryOperation<ATAN>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "atan(" << x << ")"; }
    template<typename T> static void fcn(const T& x, T& f){ f = std::atan(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0] = 1/(1+x*x);}
};

/// Step function
template<>
class UnaryOperation<STEP>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "(" << x << ">=0)"; }
    template<typename T> static void fcn(const T& x, T& f){ f = x >= 0;}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0] = 0;}
};

/// Floor function
template<>
class UnaryOperation<FLOOR>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "floor(" << x << ")"; }
    template<typename T> static void fcn(const T& x, T& f){ f = std::floor(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0] = 0;}
};

/// Ceil function
template<>
class UnaryOperation<CEIL>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "ceil(" << x << ")"; }
    template<typename T> static void fcn(const T& x, T& f){ f = std::ceil(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0] = 0;}
};

/// Equality
template<>
class BinaryOperation<EQUALITY>{
  public:
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ stream << "(" << x << "==" << y << ")"; }
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = x==y;}
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=d[1]=0;}
};

/// Minimum
template<>
class BinaryOperation<FMIN>{
  public:
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ stream << "fmin(" << x << "," << y << ")"; }
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = fmin(x,y);}
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=x<=y; d[1]=!d[0];}
};

/// Maximum
template<>
class BinaryOperation<FMAX>{
  public:
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ stream << "fmax(" << x << "," << y << ")"; }
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = fmax(x,y);}
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=x>=y; d[1]=!d[0];}
};

/// Error function
template<>
class UnaryOperation<ERF>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "erf(" << x << ")"; }
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
    
    /** \brief Derivative typedef */
    typedef void (*derT)(const T& x, const T& y, const T& f, T* d);

    /** \brief Vector of printing functions */
    static std::vector<printFunT> print;

    /** \brief Vector of function pointers to all the built in functions */
    static std::vector<funT> fun;
    
    /** \brief Vector of function derivative pointers to all the built in functions */
    static std::vector<derT> der;
    
  protected:
    
    /** \brief Create print */
    static std::vector<printFunT> getPrintFun();

    /** \brief Create fun */
    static std::vector<funT> getFun();
  
    /** \brief Create der */
    static std::vector<derT> getDer();
    
};


// Template implementations

template<typename T>
std::vector<typename casadi_math<T>::printFunT> casadi_math<T>::print = casadi_math<T>::getPrintFun();

template<typename T>
std::vector<typename casadi_math<T>::funT> casadi_math<T>::fun = casadi_math<T>::getFun();

template<typename T>
std::vector<typename casadi_math<T>::derT> casadi_math<T>::der = casadi_math<T>::getDer();

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
  
  // Specify operations
  ret[ADD] = BinaryOperation<ADD>::fcn;
  ret[SUB] = BinaryOperation<SUB>::fcn;
  ret[MUL] = BinaryOperation<MUL>::fcn;
  ret[DIV] = BinaryOperation<DIV>::fcn;
    
  ret[NEG] = BinaryOperation<NEG>::fcn;
  ret[EXP] = BinaryOperation<EXP>::fcn;
  ret[LOG] = BinaryOperation<LOG>::fcn;
  ret[POW] = BinaryOperation<POW>::fcn;

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
  ret[ADD] = BinaryOperation<ADD>::der;
  ret[SUB] = BinaryOperation<SUB>::der;
  ret[MUL] = BinaryOperation<MUL>::der;
  ret[DIV] = BinaryOperation<DIV>::der;
    
  ret[NEG] = BinaryOperation<NEG>::der;
  ret[EXP] = BinaryOperation<EXP>::der;
  ret[LOG] = BinaryOperation<LOG>::der;
  ret[POW] = BinaryOperation<POW>::der;

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
  
  // Make sure that all functions were specified
  for(int i=0; i<ret.size(); ++i){
    casadi_assert(ret[i]!=0);
  }
  
  return ret;
}


} // namespace CasADi

#endif //CASADI_MATH_HPP
