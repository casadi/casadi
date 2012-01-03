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

/// Easy access to all the functions for a particular type
template<typename T>
class casadi_math{
  public:

    /** \brief Evaluate a built in function */
    static void fun(unsigned char op, const T& x, const T& y, T& f);
    
    /** \brief Evaluate inplace a built in function (experimental) */
    static void inplacefun(unsigned char op, const T& x, const T& y, T& f);
    
    /** \brief Evaluate a built in derivative function */
    static void der(unsigned char op, const T& x, const T& y, const T& f, T* d);

    /** \brief Evaluate the function and the derivative function */
    static void derF(unsigned char op, const T& x, const T& y, T& f, T* d);
    
    /** \brief Is a function zero when evaluating with both arguments zero */
    static bool f00_is_zero(unsigned char op);
    
    /** \brief Is a function zero when evaluating with the first arguments zero */
    static bool f0x_is_zero(unsigned char op);
    
    /** \brief Is a function zero when evaluating with the second arguments zero */
    static bool fx0_is_zero(unsigned char op);
    
    /** \brief Number of dependencies */
    static int ndeps(unsigned char op);
    
    /** \brief Is a function commutative? */
    static bool isCommutative(unsigned char op);
    
    /** \brief Print */
    static void print(unsigned char op, std::ostream &stream, const std::string& x, const std::string& y);
    static void printPre(unsigned char op, std::ostream &stream);
    static void printSep(unsigned char op, std::ostream &stream);
    static void printPost(unsigned char op, std::ostream &stream);
    
};

// Template implementations

template<typename T>
inline void casadi_math<T>::fun(unsigned char op, const T& x, const T& y, T& f){
  switch(op){
    case ADD:       BinaryOperation<ADD>::fcn(x,y,f);           break;
    case SUB:       BinaryOperation<SUB>::fcn(x,y,f);           break;
    case MUL:       BinaryOperation<MUL>::fcn(x,y,f);           break;
    case DIV:       BinaryOperation<DIV>::fcn(x,y,f);           break;
    case NEG:       BinaryOperation<NEG>::fcn(x,y,f);           break;
    case EXP:       BinaryOperation<EXP>::fcn(x,y,f);           break;
    case LOG:       BinaryOperation<LOG>::fcn(x,y,f);           break;
    case POW:       BinaryOperation<POW>::fcn(x,y,f);           break;
    case CONSTPOW:  BinaryOperation<CONSTPOW>::fcn(x,y,f);      break;
    case SQRT:      BinaryOperation<SQRT>::fcn(x,y,f);          break;
    case SIN:       BinaryOperation<SIN>::fcn(x,y,f);           break;
    case COS:       BinaryOperation<COS>::fcn(x,y,f);           break;
    case TAN:       BinaryOperation<TAN>::fcn(x,y,f);           break;
    case ASIN:      BinaryOperation<ASIN>::fcn(x,y,f);          break;
    case ACOS:      BinaryOperation<ACOS>::fcn(x,y,f);          break;
    case ATAN:      BinaryOperation<ATAN>::fcn(x,y,f);          break;
    case STEP:      BinaryOperation<STEP>::fcn(x,y,f);          break;
    case FLOOR:     BinaryOperation<FLOOR>::fcn(x,y,f);         break;
    case CEIL:      BinaryOperation<CEIL>::fcn(x,y,f);          break;
    case EQUALITY:  BinaryOperation<EQUALITY>::fcn(x,y,f);      break;
    case FABS:      BinaryOperation<FABS>::fcn(x,y,f);          break;
    case SIGN:     BinaryOperation<SIGN>::fcn(x,y,f);         break;
    case ERF:       BinaryOperation<ERF>::fcn(x,y,f);           break;
    case FMIN:      BinaryOperation<FMIN>::fcn(x,y,f);          break;
    case FMAX:      BinaryOperation<FMAX>::fcn(x,y,f);          break;
    case INV:       BinaryOperation<INV>::fcn(x,y,f);           break;
    case SINH:      BinaryOperation<SINH>::fcn(x,y,f);          break;
    case COSH:      BinaryOperation<COSH>::fcn(x,y,f);          break;
    case TANH:      BinaryOperation<TANH>::fcn(x,y,f);          break;
    case PRINTME:   BinaryOperation<PRINTME>::fcn(x,y,f);       break;
  }
}

template<typename T>
inline void casadi_math<T>::inplacefun(unsigned char op, const T& x, const T& y, T& f){
#define FCASES(C,OFF) \
    case ADD+OFF:       C<ADD>::fcn(x,y,f);           break;\
    case SUB+OFF:       C<SUB>::fcn(x,y,f);           break;\
    case MUL+OFF:       C<MUL>::fcn(x,y,f);           break;\
    case DIV+OFF:       C<DIV>::fcn(x,y,f);           break;\
    case NEG+OFF:       C<NEG>::fcn(x,y,f);           break;\
    case EXP+OFF:       C<EXP>::fcn(x,y,f);           break;\
    case LOG+OFF:       C<LOG>::fcn(x,y,f);           break;\
    case POW+OFF:       C<POW>::fcn(x,y,f);           break;\
    case CONSTPOW+OFF:  C<CONSTPOW>::fcn(x,y,f);      break;\
    case SQRT+OFF:      C<SQRT>::fcn(x,y,f);          break;\
    case SIN+OFF:       C<SIN>::fcn(x,y,f);           break;\
    case COS+OFF:       C<COS>::fcn(x,y,f);           break;\
    case TAN+OFF:       C<TAN>::fcn(x,y,f);           break;\
    case ASIN+OFF:      C<ASIN>::fcn(x,y,f);          break;\
    case ACOS+OFF:      C<ACOS>::fcn(x,y,f);          break;\
    case ATAN+OFF:      C<ATAN>::fcn(x,y,f);          break;\
    case STEP+OFF:      C<STEP>::fcn(x,y,f);          break;\
    case FLOOR+OFF:     C<FLOOR>::fcn(x,y,f);         break;\
    case CEIL+OFF:      C<CEIL>::fcn(x,y,f);          break;\
    case EQUALITY+OFF:  C<EQUALITY>::fcn(x,y,f);      break;\
    case FABS+OFF:      C<FABS>::fcn(x,y,f);          break;\
    case SIGN+OFF:     C<SIGN>::fcn(x,y,f);         break;\
    case ERF+OFF:       C<ERF>::fcn(x,y,f);           break;\
    case FMIN+OFF:      C<FMIN>::fcn(x,y,f);          break;\
    case FMAX+OFF:      C<FMAX>::fcn(x,y,f);          break;\
    case INV+OFF:       C<INV>::fcn(x,y,f);           break;\
    case SINH+OFF:      C<SINH>::fcn(x,y,f);          break;\
    case COSH+OFF:      C<COSH>::fcn(x,y,f);          break;\
    case TANH+OFF:      C<TANH>::fcn(x,y,f);          break;\
    case PRINTME+OFF:   C<PRINTME>::fcn(x,y,f);       break;
    
  switch(op){
    FCASES(BinaryOperation,0)
    FCASES(AddBinaryOperation,NUM_BUILT_IN_OPS)
    FCASES(SubBinaryOperation,2*NUM_BUILT_IN_OPS)
    FCASES(MulBinaryOperation,3*NUM_BUILT_IN_OPS)
    FCASES(DivBinaryOperation,4*NUM_BUILT_IN_OPS)
  }
#undef FCASES
}

template<typename T>
inline void casadi_math<T>::der(unsigned char op, const T& x, const T& y, const T& f, T* d){
  switch(op){
    case ADD:       BinaryOperation<ADD>::der(x,y,f,d);        break;
    case SUB:       BinaryOperation<SUB>::der(x,y,f,d);        break;
    case MUL:       BinaryOperation<MUL>::der(x,y,f,d);        break;
    case DIV:       BinaryOperation<DIV>::der(x,y,f,d);        break;
    case NEG:       BinaryOperation<NEG>::der(x,y,f,d);        break;
    case EXP:       BinaryOperation<EXP>::der(x,y,f,d);        break;
    case LOG:       BinaryOperation<LOG>::der(x,y,f,d);        break;
    case POW:       BinaryOperation<POW>::der(x,y,f,d);        break;
    case CONSTPOW:  BinaryOperation<CONSTPOW>::der(x,y,f,d);   break;
    case SQRT:      BinaryOperation<SQRT>::der(x,y,f,d);       break;
    case SIN:       BinaryOperation<SIN>::der(x,y,f,d);        break;
    case COS:       BinaryOperation<COS>::der(x,y,f,d);        break;
    case TAN:       BinaryOperation<TAN>::der(x,y,f,d);        break;
    case ASIN:      BinaryOperation<ASIN>::der(x,y,f,d);       break;
    case ACOS:      BinaryOperation<ACOS>::der(x,y,f,d);       break;
    case ATAN:      BinaryOperation<ATAN>::der(x,y,f,d);       break;
    case STEP:      BinaryOperation<STEP>::der(x,y,f,d);       break;
    case FLOOR:     BinaryOperation<FLOOR>::der(x,y,f,d);      break;
    case CEIL:      BinaryOperation<CEIL>::der(x,y,f,d);       break;
    case EQUALITY:  BinaryOperation<EQUALITY>::der(x,y,f,d);   break;
    case FABS:      BinaryOperation<FABS>::der(x,y,f,d);       break;
    case SIGN:     BinaryOperation<SIGN>::der(x,y,f,d);      break;
    case ERF:       BinaryOperation<ERF>::der(x,y,f,d);        break;
    case FMIN:      BinaryOperation<FMIN>::der(x,y,f,d);       break;
    case FMAX:      BinaryOperation<FMAX>::der(x,y,f,d);       break;
    case INV:       BinaryOperation<INV>::der(x,y,f,d);        break;
    case SINH:      BinaryOperation<SINH>::der(x,y,f,d);       break;
    case COSH:      BinaryOperation<COSH>::der(x,y,f,d);       break;
    case TANH:      BinaryOperation<TANH>::der(x,y,f,d);       break;
    case PRINTME:   BinaryOperation<PRINTME>::der(x,y,f,d);    break;
  }
}


template<typename T>
inline void casadi_math<T>::derF(unsigned char op, const T& x, const T& y, T& f, T* d){
  // Copy result to temp since it might get overwritten if y or x have the same address as f
  T ff;
  
  // A lookup table
  switch(op){
    case ADD:       BinaryOperation<ADD>::fcn(x,y,ff);   BinaryOperation<ADD>::der(x,y,ff,d);        break;
    case SUB:       BinaryOperation<SUB>::fcn(x,y,ff);   BinaryOperation<SUB>::der(x,y,ff,d);        break;
    case MUL:       BinaryOperation<MUL>::fcn(x,y,ff);   BinaryOperation<MUL>::der(x,y,ff,d);        break;
    case DIV:       BinaryOperation<DIV>::fcn(x,y,ff);   BinaryOperation<DIV>::der(x,y,ff,d);        break;
    case NEG:       BinaryOperation<NEG>::fcn(x,y,ff);   BinaryOperation<NEG>::der(x,y,ff,d);        break;
    case EXP:       BinaryOperation<EXP>::fcn(x,y,ff);   BinaryOperation<EXP>::der(x,y,ff,d);        break;
    case LOG:       BinaryOperation<LOG>::fcn(x,y,ff);   BinaryOperation<LOG>::der(x,y,ff,d);        break;
    case POW:       BinaryOperation<POW>::fcn(x,y,ff);   BinaryOperation<POW>::der(x,y,ff,d);        break;
    case CONSTPOW:  BinaryOperation<CONSTPOW>::fcn(x,y,ff);   BinaryOperation<CONSTPOW>::der(x,y,ff,d);   break;
    case SQRT:      BinaryOperation<SQRT>::fcn(x,y,ff);   BinaryOperation<SQRT>::der(x,y,ff,d);       break;
    case SIN:       BinaryOperation<SIN>::fcn(x,y,ff);   BinaryOperation<SIN>::der(x,y,ff,d);        break;
    case COS:       BinaryOperation<COS>::fcn(x,y,ff);   BinaryOperation<COS>::der(x,y,ff,d);        break;
    case TAN:       BinaryOperation<TAN>::fcn(x,y,ff);   BinaryOperation<TAN>::der(x,y,ff,d);        break;
    case ASIN:      BinaryOperation<ASIN>::fcn(x,y,ff);   BinaryOperation<ASIN>::der(x,y,ff,d);       break;
    case ACOS:      BinaryOperation<ACOS>::fcn(x,y,ff);   BinaryOperation<ACOS>::der(x,y,ff,d);       break;
    case ATAN:      BinaryOperation<ATAN>::fcn(x,y,ff);   BinaryOperation<ATAN>::der(x,y,ff,d);       break;
    case STEP:      BinaryOperation<STEP>::fcn(x,y,ff);   BinaryOperation<STEP>::der(x,y,ff,d);       break;
    case FLOOR:     BinaryOperation<FLOOR>::fcn(x,y,ff);   BinaryOperation<FLOOR>::der(x,y,ff,d);      break;
    case CEIL:      BinaryOperation<CEIL>::fcn(x,y,ff);   BinaryOperation<CEIL>::der(x,y,ff,d);       break;
    case EQUALITY:  BinaryOperation<EQUALITY>::fcn(x,y,ff);   BinaryOperation<EQUALITY>::der(x,y,ff,d);   break;
    case FABS:      BinaryOperation<FABS>::fcn(x,y,ff);   BinaryOperation<FABS>::der(x,y,ff,d);        break;
    case SIGN:     BinaryOperation<FABS>::fcn(x,y,ff);   BinaryOperation<SIGN>::der(x,y,ff,d);        break;
    case ERF:       BinaryOperation<ERF>::fcn(x,y,ff);   BinaryOperation<ERF>::der(x,y,ff,d);        break;
    case FMIN:      BinaryOperation<FMIN>::fcn(x,y,ff);   BinaryOperation<FMIN>::der(x,y,ff,d);       break;
    case FMAX:      BinaryOperation<FMAX>::fcn(x,y,ff);   BinaryOperation<FMAX>::der(x,y,ff,d);       break;
    case INV:       BinaryOperation<INV>::fcn(x,y,ff);   BinaryOperation<INV>::der(x,y,ff,d);         break;
    case SINH:      BinaryOperation<SINH>::fcn(x,y,ff);   BinaryOperation<SINH>::der(x,y,ff,d);        break;
    case COSH:      BinaryOperation<COSH>::fcn(x,y,ff);   BinaryOperation<COSH>::der(x,y,ff,d);        break;
    case TANH:      BinaryOperation<TANH>::fcn(x,y,ff);   BinaryOperation<TANH>::der(x,y,ff,d);        break;
    case PRINTME:   BinaryOperation<PRINTME>::fcn(x,y,ff);   BinaryOperation<PRINTME>::der(x,y,ff,d);     break;
  }
  f = ff;
}

template<typename T>
inline bool casadi_math<T>::f00_is_zero(unsigned char op){
  switch(op){
    case ADD:       return true;
    case SUB:       return true;
    case MUL:       return true;
    case NEG:       return true;
    case SQRT:      return true;
    case SIN:       return true;
    case TAN:       return true;
    case ASIN:      return true;
    case ATAN:      return true;
    case FLOOR:     return true;
    case CEIL:      return true;
    case FMIN:      return true;
    case FMAX:      return true;
    case FABS:      return true;
    case SIGN:     return true;
    case ERF:       return true;
    case SINH:      return true;
    case TANH:      return true;
    default:        return false;
  }
}

template<typename T>
inline bool casadi_math<T>::f0x_is_zero(unsigned char op){
  switch(op){
    case MUL:       return true;
    case DIV:       return true;
    case NEG:       return true;
    case SQRT:      return true;
    case SIN:       return true;
    case TAN:       return true;
    case ASIN:      return true;
    case ATAN:      return true;
    case FLOOR:     return true;
    case CEIL:      return true;
    case FABS:      return true;
    case SIGN:     return true;
    case ERF:       return true;
    case SINH:      return true;
    case TANH:      return true;
    default:        return false;
  }
}
    
template<typename T>
inline bool casadi_math<T>::fx0_is_zero(unsigned char op){
  switch(op){
    case MUL:       return true;
    default:        return false;
  }
}
    
template<typename T>
inline bool casadi_math<T>::isCommutative(unsigned char op){
  switch(op){
    case SUB:       return false;
    case DIV:       return false;
    case POW:       return false;
    case CONSTPOW:  return false;
    case EQUALITY:  return false;
    case PRINTME:   return false;
    default:        return true;
  }
}

template<typename T>
inline int casadi_math<T>::ndeps(unsigned char op){
  switch(op){
    case ADD:       return 2;
    case SUB:       return 2;
    case MUL:       return 2;
    case DIV:       return 2;
    case POW:       return 2;
    case CONSTPOW:  return 2;
    case EQUALITY:  return 2;
    case FMIN:      return 2;
    case FMAX:      return 2;
    case PRINTME:   return 2;
    default:        return 1;
  }
}

template<typename T>
inline void casadi_math<T>::print(unsigned char op, std::ostream &stream, const std::string& x, const std::string& y){
  if(ndeps(op)==2){
    printPre(op,stream);
    stream << x;
    printSep(op,stream);
    stream << y;
    printPost(op,stream);
  } else {
    printPre(op,stream);
    stream << x;
    printPost(op,stream);
  }
}

template<typename T>
inline void casadi_math<T>::printPre(unsigned char op, std::ostream &stream){
  switch(op){
    case ADD:       stream << "(";        break;
    case SUB:       stream << "(";        break;
    case MUL:       stream << "(";        break;
    case DIV:       stream << "(";        break;
    case NEG:       stream << "(-";       break;
    case EXP:       stream << "exp(";     break;
    case LOG:       stream << "log(";     break;
    case POW:       stream << "pow(";     break;
    case CONSTPOW:  stream << "pow(";     break;
    case SQRT:      stream << "sqrt(";    break;
    case SIN:       stream << "sin(";     break;
    case COS:       stream << "cos(";     break;
    case TAN:       stream << "tan(";     break;
    case ASIN:      stream << "asin(";    break;
    case ACOS:      stream << "acos(";    break;
    case ATAN:      stream << "atan(";    break;
    case STEP:      stream << "(";        break;
    case FLOOR:     stream << "floor(";   break;
    case CEIL:      stream << "ceil(";    break;
    case EQUALITY:  stream << "(";        break;
    case FABS:      stream << "fabs(";    break;
    case SIGN:     stream << "sign(";   break;
    case ERF:       stream << "erf(";     break;
    case FMIN:      stream << "fmin(";    break;
    case FMAX:      stream << "fmax(";    break;
    case INV:       stream << "(1./";     break;
    case SINH:      stream << "sinh(";    break;
    case COSH:      stream << "cosh(";    break;
    case TANH:      stream << "tanh(";    break;
    case PRINTME:   stream << "printme("; break;
  }
}

template<typename T>
inline void casadi_math<T>::printSep(unsigned char op, std::ostream &stream){
  switch(op){
    case ADD:       stream << "+";        break;
    case SUB:       stream << "-";        break;
    case MUL:       stream << "*";        break;
    case DIV:       stream << "/";        break;
    case EQUALITY:  stream << "==";       break;
    default:        stream << ",";        break;
  }
}

template<typename T>
inline void casadi_math<T>::printPost(unsigned char op, std::ostream &stream){
  switch(op){
    case STEP:      stream << ">=0)";     break;
    default:        stream << ")";        break;
  }
}

} // namespace CasADi

#endif //CASADI_MATH_HPP
