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

    /** \brief Print part
    *
    * print(stream, x) is equivalent to \n
    *
    * printPre(stream) \n
    * stream << x;
    * printPost(stream) \n
    */
    static void printPre(std::ostream &stream);
    static void printPost(std::ostream &stream);
    
    static int ndeps() { return 1;}
    
    /// Does evaluating the function with a zero give a zero
    static bool f0_is_zero(){return false;}

    /// Function evaluation
    template<typename T> static void fcn(const T& x, T& f);
    
    /// Partial derivatives
    template<typename T> static void der(const T& x, const T& f, T* d);
};

#define printRoutinesUnary(pre,post) \
  static void print(std::ostream &stream, const std::string& x){ stream << pre << x << post;  } \
  static void printPre (std::ostream &stream) { stream << pre; } \
  static void printPost(std::ostream &stream) { stream << post; } 
  
#define printRoutines(pre,sep,post) \
  static void print(std::ostream &stream, const std::string& x, const std::string& y){ stream << pre << x << sep << y << post;  } \
  static void printPre (std::ostream &stream) { stream << pre; } \
  static void printSep (std::ostream &stream) { stream << sep; } \
  static void printPost(std::ostream &stream) { stream << post; } 

template<int I>
class BinaryOperation{
  public:
    /// Print
    static void print(std::ostream &stream, const std::string& x, const std::string& y){ UnaryOperation<I>::printPre(stream); stream << x; BinaryOperation<I>::printPost(stream); }
    /** \brief Print part
    *
    * print(stream, x,y) is equivalent to \n
    *
    * printPre(stream) \n
    * stream << x;
    * printSep(stream) \n
    * stream << y;
    * printPost(stream) \n
    */
    static void printPre (std::ostream &stream) { UnaryOperation<I>::printPre(stream); }
    static void printSep (std::ostream &stream) { return; }
    static void printPost(std::ostream &stream) { UnaryOperation<I>::printPost(stream); }
    
    /// Does evaluating the function with a zero give a zero
    static bool f00_is_zero(){return UnaryOperation<I>::f0_is_zero();}
    static bool f0x_is_zero(){return UnaryOperation<I>::f0_is_zero();}
    static bool fx0_is_zero(){return false;}

    static int ndeps() { return UnaryOperation<I>::ndeps();}
    
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
  INV,
  SINH,  COSH,  TANH,
  PRINTME,
  NUM_BUILT_IN_OPS
};

/// Addition
template<>
class BinaryOperation<ADD>{
  public:
    printRoutines("(","+",")")
    static bool f00_is_zero(){return true;}
    static bool f0x_is_zero(){return false;}
    static bool fx0_is_zero(){return false;}
    static int ndeps() { return 2;}
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = x+y;}
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=d[1]=1;}
};

/// Subtraction
template<>
class BinaryOperation<SUB>{
  public:
    printRoutines("(","-",")")
    static bool f00_is_zero(){return true;}
    static bool f0x_is_zero(){return false;}
    static bool fx0_is_zero(){return false;}
    static int ndeps() { return 2;}
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = x-y;}
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=1; d[1]=-1;}
};

/// Multiplication
template<>
class BinaryOperation<MUL>{
  public:
    printRoutines("(","*",")")
    static bool f00_is_zero(){return true;}
    static bool f0x_is_zero(){return true;}
    static bool fx0_is_zero(){return true;}
    static int ndeps() { return 2;}
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = x*y;}
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=y; d[1]=x;}
};

/// Division
template<>
class BinaryOperation<DIV>{
  public:
    printRoutines("(","/",")")
    static bool f00_is_zero(){return false;}
    static bool f0x_is_zero(){return true;}
    static bool fx0_is_zero(){return false;}
    static int ndeps() { return 2;}
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = x/y;}
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=1/y; d[1]=-f/y;}
};

/// Negation
template<>
class UnaryOperation<NEG>{
  public:
    printRoutinesUnary("(-",")")
    static bool f0_is_zero(){return true;}
    static int ndeps() { return 1;}
    template<typename T> static void fcn(const T& x, T& f){ f = -x;}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0]=-1;}
};

/// Natural exponent
template<>
class UnaryOperation<EXP>{
  public:
    printRoutinesUnary("exp(",")")
    static bool f0_is_zero(){return false;}
    static int ndeps() { return 1;}
    template<typename T> static void fcn(const T& x, T& f){ f = std::exp(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0]=f;}
};

/// Natural logarithm
template<>
class UnaryOperation<LOG>{
  public:
    printRoutinesUnary("log(",")")
    static bool f0_is_zero(){return false;}
    static int ndeps() { return 1;}
    template<typename T> static void fcn(const T& x, T& f){ f = std::log(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0]=1/x;}
};

/// Power, defined only for x>=0
template<>
class BinaryOperation<POW>{
  public:
    printRoutines("pow(",",",")")
    static bool f00_is_zero(){return false;}
    static bool f0x_is_zero(){return false;}
    static bool fx0_is_zero(){return false;}
    static int ndeps() { return 2;}
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = std::pow(x,y);}
    // See issue #104 why d[0] is no longer y*f/x
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=y*std::pow(x,y-1); d[1]=std::log(x)*f;}
};

/// Power, defined only for y constant
template<>
class BinaryOperation<CONSTPOW>{
  public:
    printRoutines("pow(",",",")")
    static bool f00_is_zero(){return false;}
    static bool f0x_is_zero(){return false;}
    static bool fx0_is_zero(){return false;}
    static int ndeps() { return 2;}
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = std::pow(x,y);}
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=y*std::pow(x,y-1); d[1]=0;}
};

/// Square root
template<>
class UnaryOperation<SQRT>{
  public:
    printRoutinesUnary("sqrt(",")")
    static bool f0_is_zero(){return true;}
    static int ndeps() { return 1;}
    template<typename T> static void fcn(const T& x, T& f){ f = std::sqrt(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0]=1/(timesTwo(f));}
};

/// Sine
template<>
class UnaryOperation<SIN>{
  public:
    printRoutinesUnary("sin(",")")
    static bool f0_is_zero(){return true;}
    static int ndeps() { return 1;}
    template<typename T> static void fcn(const T& x, T& f){ f = std::sin(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0]=std::cos(x);}
};

/// Cosine
template<>
class UnaryOperation<COS>{
  public:
    printRoutinesUnary("cos(",")")
    static bool f0_is_zero(){return false;}
    static int ndeps() { return 1;}
    template<typename T> static void fcn(const T& x, T& f){ f = std::cos(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0]=-std::sin(x);}
};

/// Tangent
template<>
class UnaryOperation<TAN>{
  public:
    printRoutinesUnary("tan(",")")
    static bool f0_is_zero(){return true;}
    static int ndeps() { return 1;}
    template<typename T> static void fcn(const T& x, T& f){ f = std::tan(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0] = 1/square(std::cos(x));}
};

/// Arcus sine
template<>
class UnaryOperation<ASIN>{
  public:
    printRoutinesUnary("asin(",")")
    static bool f0_is_zero(){return true;}
    static int ndeps() { return 1;}
    template<typename T> static void fcn(const T& x, T& f){ f = std::asin(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0]=1/std::sqrt(1-x*x);}
};

/// Arcus cosine
template<>
class UnaryOperation<ACOS>{
  public:
    printRoutinesUnary("acos(",")")
    static bool f0_is_zero(){return false;}
    static int ndeps() { return 1;}
    template<typename T> static void fcn(const T& x, T& f){ f = std::acos(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0]=-1/std::sqrt(1-x*x);}
};

/// Arcus tangent
template<>
class UnaryOperation<ATAN>{
  public:
    printRoutinesUnary("atan(",")")
    static bool f0_is_zero(){return true;}
    static int ndeps() { return 1;}
    template<typename T> static void fcn(const T& x, T& f){ f = std::atan(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0] = 1/(1+x*x);}
};

/// Step function
template<>
class UnaryOperation<STEP>{
  public:
    printRoutinesUnary("(",">=0)")
    static bool f0_is_zero(){return false;}
    static int ndeps() { return 1;}
    template<typename T> static void fcn(const T& x, T& f){ f = x >= T(0.);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0] = 0;}
};

/// Floor function
template<>
class UnaryOperation<FLOOR>{
  public:
    printRoutinesUnary("floor(",")")
    static bool f0_is_zero(){return true;}
    static int ndeps() { return 1;}
    template<typename T> static void fcn(const T& x, T& f){ f = std::floor(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0] = 0;}
};

/// Ceil function
template<>
class UnaryOperation<CEIL>{
  public:
    printRoutinesUnary("ceil(",")")
    static bool f0_is_zero(){return true;}
    static int ndeps() { return 1;}
    template<typename T> static void fcn(const T& x, T& f){ f = std::ceil(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0] = 0;}
};

/// Equality
template<>
class BinaryOperation<EQUALITY>{
  public:
    printRoutines("(","==",")")
    static bool f00_is_zero(){return false;}
    static bool f0x_is_zero(){return false;}
    static bool fx0_is_zero(){return false;}
    static int ndeps() { return 2;}
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = x==y;}
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=d[1]=0;}
};

/// Minimum
template<>
class BinaryOperation<FMIN>{
  public:
    printRoutines("fmin(",",",")")
    static bool f00_is_zero(){return true;}
    static bool f0x_is_zero(){return false;}
    static bool fx0_is_zero(){return false;}
    static int ndeps() { return 2;}
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = fmin(x,y);}
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=x<=y; d[1]=!d[0];}
};

/// Maximum
template<>
class BinaryOperation<FMAX>{
  public:
    printRoutines("fmax(",",",")")
    static bool f00_is_zero(){return true;}
    static bool f0x_is_zero(){return false;}
    static bool fx0_is_zero(){return false;}
    static int ndeps() { return 2;}
    template<typename T> static void fcn(const T& x, const T& y, T& f){ f = fmax(x,y);}
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=x>=y; d[1]=!d[0];}
};

/// Error function
template<>
class UnaryOperation<ERF>{
  public:
    printRoutinesUnary("erf(",")")
    static bool f0_is_zero(){return true;}
    static int ndeps() { return 1;}
    template<typename T> static void fcn(const T& x, T& f){ f = erf(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0] = (2/std::sqrt(M_PI))*std::exp(-x*x);}
};

/// Elementwise inverse
template<>
class UnaryOperation<INV>{
  public:
    printRoutinesUnary("(1/",")")
    static bool f0_is_zero(){return false;}
    static int ndeps() { return 1;}
    template<typename T> static void fcn(const T& x, T& f){ f = 1./x;}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0] = -f*f; }
};

/// Hyperbolic sine
template<>
class UnaryOperation<SINH>{
  public:
    printRoutinesUnary("sinh(",")")
    static bool f0_is_zero(){return true;}
    static int ndeps() { return 1;}
    template<typename T> static void fcn(const T& x, T& f){ f = std::sinh(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0] = std::cosh(x); }
};

/// Hyperbolic cosine
template<>
class UnaryOperation<COSH>{
  public:
    printRoutinesUnary("cosh(",")")
    static bool f0_is_zero(){return false;}
    static int ndeps() { return 1;}
    template<typename T> static void fcn(const T& x, T& f){ f = std::cosh(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0] = -std::sinh(x); }
};

/// Hyperbolic tangent
template<>
class UnaryOperation<TANH>{
  public:
    printRoutinesUnary("tanh(",")")
    static bool f0_is_zero(){return true;}
    static int ndeps() { return 1;}
    template<typename T> static void fcn(const T& x, T& f){ f = std::tanh(x);}
    template<typename T> static void der(const T& x, const T& f, T* d){ d[0] = 1-f*f; }
};

/// Identity operator with the side effect of printing
template<>
class BinaryOperation<PRINTME>{
  public:
    printRoutines("printme(",",",")")
    static bool f00_is_zero(){return false;}
    static bool f0x_is_zero(){return false;}
    static bool fx0_is_zero(){return false;}
    static int ndeps() { return 2;}
    template<typename T> static void fcn(const T& x, const T& y, T& f){f = x; }
    
    template<typename T> static void der(const T& x, const T& y, const T& f, T* d){ d[0]=1; d[1]=0;}
};


#ifdef WITH_PRINTME 
template<> inline
void BinaryOperation<PRINTME>::fcn<double>(const double& x, const double& y, double& f) {
       f = x; 
       std::cout << "|> " << y << " : " << x << std::endl;
}
#endif //WITH_PRINTME 

/// Easy access to all the functions for a particular type
template<typename T>
class casadi_math{
  public:

    /** \brief Printing operation typedef */
    typedef void (*printFunT)(std::ostream &stream, const std::string& x, const std::string& y);

    /** \brief Printing operation typedef */
    typedef void (*printCompFunT)(std::ostream &stream);
    
    /** \brief Function typedef */
    typedef void (*funT)(const T&, const T&, T&);
    typedef T (*funTE)(const T&, const T&);
        
    /** \brief Derivative typedef */
    typedef void (*derT)(const T& x, const T& y, const T& f, T* d);
    
    /** \brief Vector of printing functions */
    static std::vector<printFunT> print;
    
    static std::vector<printCompFunT> printPre;
    static std::vector<printCompFunT> printSep;
    static std::vector<printCompFunT> printPost;

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
    
    /** \brief Vector of ints indicating the number of dependencies */
    static std::vector<int> ndeps;
    
  protected:
    
    /** \brief Create print */
    static std::vector<printFunT> getPrintFun();
    static std::vector<printCompFunT> getPrintPreFun();
    static std::vector<printCompFunT> getPrintSepFun();
    static std::vector<printCompFunT> getPrintPostFun();
    
    /** \brief Create fun */
    static std::vector<funT> getFun();
    static std::vector<funTE> getFunE();
  
    /** \brief Create der */
    static std::vector<derT> getDer();
    
    /** \brief Create zero indicator vectors */
    static std::vector<bool> getF00_is_zero();
    static std::vector<bool> getF0x_is_zero();
    static std::vector<bool> getFx0_is_zero();
    
    static std::vector<int> getNdeps();
    
};


// Template implementations

template<typename T>
std::vector<typename casadi_math<T>::printFunT> casadi_math<T>::print = casadi_math<T>::getPrintFun();

template<typename T>
std::vector<typename casadi_math<T>::printCompFunT> casadi_math<T>::printPre = casadi_math<T>::getPrintPreFun();

template<typename T>
std::vector<typename casadi_math<T>::printCompFunT> casadi_math<T>::printSep = casadi_math<T>::getPrintSepFun();

template<typename T>
std::vector<typename casadi_math<T>::printCompFunT> casadi_math<T>::printPost = casadi_math<T>::getPrintPostFun();

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
std::vector<int> casadi_math<T>::ndeps = casadi_math<T>::getNdeps();

#define populatePrintFun(getPrintFunName,printFunName,Type) \
template<typename T> \
std::vector<Type> casadi_math<T>::getPrintFunName(){ \
 \
  std::vector<Type> ret(NUM_BUILT_IN_OPS,0); \
   \
  ret[ADD] = BinaryOperation<ADD>::printFunName; \
  ret[SUB] = BinaryOperation<SUB>::printFunName; \
  ret[MUL] = BinaryOperation<MUL>::printFunName; \
  ret[DIV] = BinaryOperation<DIV>::printFunName; \
   \
  ret[NEG] = BinaryOperation<NEG>::printFunName; \
  ret[EXP] = BinaryOperation<EXP>::printFunName; \
  ret[LOG] = BinaryOperation<LOG>::printFunName; \
  ret[POW] = BinaryOperation<POW>::printFunName; \
  ret[CONSTPOW] = BinaryOperation<CONSTPOW>::printFunName; \
 \
  ret[SQRT] = BinaryOperation<SQRT>::printFunName; \
  ret[SIN] = BinaryOperation<SIN>::printFunName; \
  ret[COS] = BinaryOperation<COS>::printFunName; \
  ret[TAN] = BinaryOperation<TAN>::printFunName; \
 \
  ret[ASIN] = BinaryOperation<ASIN>::printFunName; \
  ret[ACOS] = BinaryOperation<ACOS>::printFunName; \
  ret[ATAN] = BinaryOperation<ATAN>::printFunName; \
 \
  ret[STEP] = BinaryOperation<STEP>::printFunName; \
  ret[FLOOR] = BinaryOperation<FLOOR>::printFunName; \
  ret[CEIL] = BinaryOperation<CEIL>::printFunName; \
 \
  ret[EQUALITY] = BinaryOperation<EQUALITY>::printFunName; \
 \
  ret[ERF] = BinaryOperation<ERF>::printFunName; \
  ret[FMIN] = BinaryOperation<FMIN>::printFunName; \
  ret[FMAX] = BinaryOperation<FMAX>::printFunName; \
  ret[INV] = BinaryOperation<INV>::printFunName; \
 \
  ret[SINH] = BinaryOperation<SINH>::printFunName; \
  ret[COSH] = BinaryOperation<COSH>::printFunName; \
  ret[TANH] = BinaryOperation<TANH>::printFunName; \
 \
   ret[PRINTME] = BinaryOperation<PRINTME>::printFunName; \
 \
  for(int i=0; i<ret.size(); ++i){ \
    casadi_assert(ret[i]!=0); \
  } \
   \
  return ret; \
   \
}

populatePrintFun(getPrintFun,print,typename casadi_math<T>::printFunT)
populatePrintFun(getPrintPreFun,printPre,typename casadi_math<T>::printCompFunT)
populatePrintFun(getPrintSepFun,printSep,typename casadi_math<T>::printCompFunT)
populatePrintFun(getPrintPostFun,printPost,typename casadi_math<T>::printCompFunT)

populatePrintFun(getNdeps,ndeps(),int)

template<typename T>
std::vector<typename casadi_math<T>::funT> casadi_math<T>::getFun(){
  
  // Create return object
  std::vector<typename casadi_math<T>::funT> ret(NUM_BUILT_IN_OPS,0);
  
  // Specify operations
  ret[ADD] = BinaryOperation<ADD>::fcn<T>;
  ret[SUB] = BinaryOperation<SUB>::fcn<T>;
  ret[MUL] = BinaryOperation<MUL>::fcn<T>;
  ret[DIV] = BinaryOperation<DIV>::fcn<T>;
    
  ret[NEG] = BinaryOperation<NEG>::fcn<T>;
  ret[EXP] = BinaryOperation<EXP>::fcn<T>;
  ret[LOG] = BinaryOperation<LOG>::fcn<T>;
  ret[POW] = BinaryOperation<POW>::fcn<T>;
  ret[CONSTPOW] = BinaryOperation<CONSTPOW>::fcn<T>;

  ret[SQRT] = BinaryOperation<SQRT>::fcn<T>;
  ret[SIN] = BinaryOperation<SIN>::fcn<T>;
  ret[COS] = BinaryOperation<COS>::fcn<T>;
  ret[TAN] = BinaryOperation<TAN>::fcn<T>;

  ret[ASIN] = BinaryOperation<ASIN>::fcn<T>;
  ret[ACOS] = BinaryOperation<ACOS>::fcn<T>;
  ret[ATAN] = BinaryOperation<ATAN>::fcn<T>;

  ret[STEP] = BinaryOperation<STEP>::fcn<T>;
  ret[FLOOR] = BinaryOperation<FLOOR>::fcn<T>;
  ret[CEIL] = BinaryOperation<CEIL>::fcn<T>;

  ret[EQUALITY] = BinaryOperation<EQUALITY>::fcn<T>;

  ret[ERF] = BinaryOperation<ERF>::fcn<T>;
  ret[FMIN] = BinaryOperation<FMIN>::fcn<T>;
  ret[FMAX] = BinaryOperation<FMAX>::fcn<T>;

  ret[INV] = BinaryOperation<INV>::fcn<T>;

  ret[SINH] = BinaryOperation<SINH>::fcn<T>;
  ret[COSH] = BinaryOperation<COSH>::fcn<T>;
  ret[TANH] = BinaryOperation<TANH>::fcn<T>;

  ret[PRINTME] = BinaryOperation<PRINTME>::fcn<T>;
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
  
  // Specify operations
  ret[ADD] = BinaryOperationE<ADD>::fcn<T>;
  ret[SUB] = BinaryOperationE<SUB>::fcn<T>;
  ret[MUL] = BinaryOperationE<MUL>::fcn<T>;
  ret[DIV] = BinaryOperationE<DIV>::fcn<T>;
    
  ret[NEG] = BinaryOperationE<NEG>::fcn<T>;
  ret[EXP] = BinaryOperationE<EXP>::fcn<T>;
  ret[LOG] = BinaryOperationE<LOG>::fcn<T>;
  ret[POW] = BinaryOperationE<POW>::fcn<T>;
  ret[CONSTPOW] = BinaryOperationE<CONSTPOW>::fcn<T>;

  ret[SQRT] = BinaryOperationE<SQRT>::fcn<T>;
  ret[SIN] = BinaryOperationE<SIN>::fcn<T>;
  ret[COS] = BinaryOperationE<COS>::fcn<T>;
  ret[TAN] = BinaryOperationE<TAN>::fcn<T>;

  ret[ASIN] = BinaryOperationE<ASIN>::fcn<T>;
  ret[ACOS] = BinaryOperationE<ACOS>::fcn<T>;
  ret[ATAN] = BinaryOperationE<ATAN>::fcn<T>;

  ret[STEP] = BinaryOperationE<STEP>::fcn<T>;
  ret[FLOOR] = BinaryOperationE<FLOOR>::fcn<T>;
  ret[CEIL] = BinaryOperationE<CEIL>::fcn<T>;

  ret[EQUALITY] = BinaryOperationE<EQUALITY>::fcn<T>;

  ret[ERF] = BinaryOperationE<ERF>::fcn<T>;
  ret[FMIN] = BinaryOperationE<FMIN>::fcn<T>;
  ret[FMAX] = BinaryOperationE<FMAX>::fcn<T>;

  ret[INV] = BinaryOperationE<INV>::fcn<T>;

  ret[SINH] = BinaryOperationE<SINH>::fcn<T>;
  ret[COSH] = BinaryOperationE<COSH>::fcn<T>;
  ret[TANH] = BinaryOperationE<TANH>::fcn<T>;

  ret[PRINTME] = BinaryOperationE<PRINTME>::fcn<T>;
  
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
  ret[ADD] = BinaryOperation<ADD>::der<T>;
  ret[SUB] = BinaryOperation<SUB>::der<T>;
  ret[MUL] = BinaryOperation<MUL>::der<T>;
  ret[DIV] = BinaryOperation<DIV>::der<T>;
    
  ret[NEG] = BinaryOperation<NEG>::der<T>;
  ret[EXP] = BinaryOperation<EXP>::der<T>;
  ret[LOG] = BinaryOperation<LOG>::der<T>;
  ret[POW] = BinaryOperation<POW>::der<T>;
  ret[CONSTPOW] = BinaryOperation<CONSTPOW>::der<T>;

  ret[SQRT] = BinaryOperation<SQRT>::der<T>;
  ret[SIN] = BinaryOperation<SIN>::der<T>;
  ret[COS] = BinaryOperation<COS>::der<T>;
  ret[TAN] = BinaryOperation<TAN>::der<T>;

  ret[ASIN] = BinaryOperation<ASIN>::der<T>;
  ret[ACOS] = BinaryOperation<ACOS>::der<T>;
  ret[ATAN] = BinaryOperation<ATAN>::der<T>;

  ret[STEP] = BinaryOperation<STEP>::der<T>;
  ret[FLOOR] = BinaryOperation<FLOOR>::der<T>;
  ret[CEIL] = BinaryOperation<CEIL>::der<T>;

  ret[EQUALITY] = BinaryOperation<EQUALITY>::der<T>;

  ret[ERF] = BinaryOperation<ERF>::der<T>;
  ret[FMIN] = BinaryOperation<FMIN>::der<T>;
  ret[FMAX] = BinaryOperation<FMAX>::der<T>;

  ret[INV] = BinaryOperation<INV>::der<T>;

  ret[SINH] = BinaryOperation<SINH>::der<T>;
  ret[COSH] = BinaryOperation<COSH>::der<T>;
  ret[TANH] = BinaryOperation<TANH>::der<T>;

  ret[PRINTME] = BinaryOperation<PRINTME>::der<T>;



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

  ret[INV] = BinaryOperation<INV>::f00_is_zero();

  ret[SINH] = BinaryOperation<SINH>::f00_is_zero();
  ret[COSH] = BinaryOperation<COSH>::f00_is_zero();
  ret[TANH] = BinaryOperation<TANH>::f00_is_zero();

  ret[PRINTME] = BinaryOperation<PRINTME>::f00_is_zero();

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

  ret[INV] = BinaryOperation<INV>::f0x_is_zero();

  ret[SINH] = BinaryOperation<SINH>::f0x_is_zero();
  ret[COSH] = BinaryOperation<COSH>::f0x_is_zero();
  ret[TANH] = BinaryOperation<TANH>::f0x_is_zero();

  ret[PRINTME] = BinaryOperation<PRINTME>::f0x_is_zero();
  
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

  ret[INV] = BinaryOperation<INV>::fx0_is_zero();

  ret[SINH] = BinaryOperation<SINH>::fx0_is_zero();
  ret[COSH] = BinaryOperation<COSH>::fx0_is_zero();
  ret[TANH] = BinaryOperation<TANH>::fx0_is_zero();

  ret[PRINTME] = BinaryOperation<PRINTME>::fx0_is_zero();

  return ret;
}


} // namespace CasADi

#endif //CASADI_MATH_HPP
