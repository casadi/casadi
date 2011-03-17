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

#ifndef ELEMENTARY_FUNCTIONS_HPP
#define ELEMENTARY_FUNCTIONS_HPP

#include <iostream>
#include <string>
#include <cmath>

namespace CasADi{
  
template<int I>
class UnaryOperation{
  public:
    /// Print
    static void print(std::ostream &stream, const std::string& x);

    /// Function evaluation
    template<typename T> static T fcn(const T& x);
    
    /// Partial derivatives
    template<typename T> static void der(const T& x, T& f, T* d);
};

template<int I>
class BinaryOperation{
  public:
    /// Print
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ UnaryOperation<I>::print(stream,x); }

    /// Function evaluation
    template<typename T> static T fcn(const T& x, const T& y){ return UnaryOperation<I>::fcn(x);}
    
    /// Partial derivatives - binary function
    template<typename T> static void der(const T& x, const T& y, T& f, T* d){ UnaryOperation<I>::der(x,f,d); d[1]=0; }
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
    template<typename T> static T fcn(const T& x, const T& y){ return x+y;}
    template<typename T> static void der(const T& x, const T& y, T& f, T* d){ f = x+y; d[0]=d[1]=1;}
};

/// Subtraction
template<>
class BinaryOperation<SUB>{
  public:
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ stream << "(" << x << "-" << y << ")"; }
    template<typename T> static T fcn(const T& x, const T& y){ return x-y;}
    template<typename T> static void der(const T& x, const T& y, T& f, T* d){ f = x-y; d[0]=1; d[1]=-1;}
};

/// Multiplication
template<>
class BinaryOperation<MUL>{
  public:
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ stream << "(" << x << "*" << y << ")"; }
    template<typename T> static T fcn(const T& x, const T& y){ return x*y;}
    template<typename T> static void der(const T& x, const T& y, T& f, T* d){ f = x*y; d[0]=y; d[1]=x;}
};

/// Division
template<>
class BinaryOperation<DIV>{
  public:
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ stream << "(" << x << "/" << y << ")"; }
    template<typename T> static T fcn(const T& x, const T& y){ return x/y;}
    template<typename T> static void der(const T& x, const T& y, T& f, T* d){ f = x/y; d[0]=1/y; d[1]=-f/y;}
};

/// Negation
template<>
class UnaryOperation<NEG>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "(-" << x << ")"; }
    template<typename T> static T fcn(const T& x){ return -x;}
    template<typename T> static void der(const T& x, T& f, T* d){ f = -x; d[0]=-1;}
};

/// Natural exponent
template<>
class UnaryOperation<EXP>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "exp(" << x << ")"; }
    template<typename T> static T fcn(const T& x){ return std::exp(x);}
    template<typename T> static void der(const T& x, T& f, T* d){ f = std::exp(x); d[0]=f;}
};

/// Natural logarithm
template<>
class UnaryOperation<LOG>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "log(" << x << ")"; }
    template<typename T> static T fcn(const T& x){ return std::log(x);}
    template<typename T> static void der(const T& x, T& f, T* d){ f = std::log(x); d[0]=1/x;}
};

/// Power
template<>
class BinaryOperation<POW>{
  public:
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ stream << "pow(" << x << "," << y << ")"; }
    template<typename T> static T fcn(const T& x, const T& y){ return std::pow(x,y);}
    template<typename T> static void der(const T& x, const T& y, T& f, T* d){ f = std::pow(x,y); d[0]=y*f/x; d[1]=std::log(x)*f;}
};

/// Square root
template<>
class UnaryOperation<SQRT>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "sqrt(" << x << ")"; }
    template<typename T> static T fcn(const T& x){ return std::sqrt(x);}
    template<typename T> static void der(const T& x, T& f, T* d){ f = std::sqrt(x); d[0]=1/(2*f);}
};

/// Sine
template<>
class UnaryOperation<SIN>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "sin(" << x << ")"; }
    template<typename T> static T fcn(const T& x){ return std::sin(x);}
    template<typename T> static void der(const T& x, T& f, T* d){ f = std::sin(x); d[0]=std::cos(x);}
};

/// Cosine
template<>
class UnaryOperation<COS>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "cos(" << x << ")"; }
    template<typename T> static T fcn(const T& x){ return std::cos(x);}
    template<typename T> static void der(const T& x, T& f, T* d){ f = std::cos(x); d[0]=-std::sin(x);}
};

/// Tangent
template<>
class UnaryOperation<TAN>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "tan(" << x << ")"; }
    template<typename T> static T fcn(const T& x){ return std::tan(x);}
    template<typename T> static void der(const T& x, T& f, T* d){ f = std::tan(x); T cosx = std::cos(x); d[0] = 1/(cosx*cosx);}
};

/// Arcus sine
template<>
class UnaryOperation<ASIN>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "asin(" << x << ")"; }
    template<typename T> static T fcn(const T& x){ return std::asin(x);}
    template<typename T> static void der(const T& x, T& f, T* d){ f = std::asin(x); d[0]=1/std::sqrt(1-x*x);}
};

/// Arcus cosine
template<>
class UnaryOperation<ACOS>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "acos(" << x << ")"; }
    template<typename T> static T fcn(const T& x){ return std::acos(x);}
    template<typename T> static void der(const T& x, T& f, T* d){ f = std::acos(x); d[0]=-1/std::sqrt(1-x*x);}
};

/// Arcus tangent
template<>
class UnaryOperation<ATAN>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "atan(" << x << ")"; }
    template<typename T> static T fcn(const T& x){ return std::atan(x);}
    template<typename T> static void der(const T& x, T& f, T* d){ f = std::atan(x); d[0] = 1/(1+x*x);}
};

/// Step function
template<>
class UnaryOperation<STEP>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "(" << x << ">=0)"; }
    template<typename T> static T fcn(const T& x){ return x >= 0;}
    template<typename T> static void der(const T& x, T& f, T* d){ f = x >= 0; d[0] = 0;}
};

/// Floor function
template<>
class UnaryOperation<FLOOR>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "floor(" << x << ")"; }
    template<typename T> static T fcn(const T& x){ return std::floor(x);}
    template<typename T> static void der(const T& x, T& f, T* d){ f = std::floor(x); d[0] = 0;}
};

/// Ceil function
template<>
class UnaryOperation<CEIL>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "ceil(" << x << ")"; }
    template<typename T> static T fcn(const T& x){ return std::ceil(x);}
    template<typename T> static void der(const T& x, T& f, T* d){ f = std::ceil(x); d[0] = 0;}
};

/// Equality
template<>
class BinaryOperation<EQUALITY>{
  public:
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ stream << "(" << x << "==" << y << ")"; }
    template<typename T> static T fcn(const T& x, const T& y){ return x==y;}
    template<typename T> static void der(const T& x, const T& y, T& f, T* d){ f = x==y; d[0]=d[1]=0;}
};

/// Minimum
template<>
class BinaryOperation<FMIN>{
  public:
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ stream << "fmin(" << x << "," << y << ")"; }
    template<typename T> static T fcn(const T& x, const T& y){ return fmin(x,y);}
    template<typename T> static void der(const T& x, const T& y, T& f, T* d){ f = fmin(x,y); d[0]=x<=y; d[1]=!d[0];}
};

/// Maximum
template<>
class BinaryOperation<FMAX>{
  public:
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ stream << "fmax(" << x << "," << y << ")"; }
    template<typename T> static T fcn(const T& x, const T& y){ return fmax(x,y);}
    template<typename T> static void der(const T& x, const T& y, T& f, T* d){ f = fmax(x,y); d[0]=x>=y; d[1]=!d[0];}
};

/// Error function
template<>
class UnaryOperation<ERF>{
  public:
    static void print(std::ostream &stream, const std::string& x){ stream << "erf(" << x << ")"; }
    template<typename T> static T fcn(const T& x){ return erf(x);}
    template<typename T> static void der(const T& x, T& f, T* d){ f = erf(x); d[0] = (2/std::sqrt(M_PI))*std::exp(-x*x);}
};

} // namespace CasADi

#endif //ELEMENTARY_FUNCTIONS_HPP
