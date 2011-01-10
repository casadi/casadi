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

/** \brief  C/C++ */
#include <iostream>
//#include <fstream>
#include <string>
#include <sstream>
#include <limits>
#include <cmath>
#include <vector>

namespace CasADi{

  /** \brief  forward declaration of Node and Matrix */
class SXNode; // include will follow in the end
class SXMatrix; 

/** \brief The basic scalar symbolic class of CasADi
  \author Joel Andersson 
  \date 2010
*/ 

#ifdef SWIG
#ifdef WITH_IMPLICITCONV
%implicitconv SX;
#endif // WITH_IMPLICITCONV
#endif // SWIG


class SX{
  friend class SXNode;

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
 	\param Name of the symbol

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
    explicit SX(const char name[]);  // variable

    explicit SX(SXNode* node); // (must be explicit, otherwise 0/NULL would be ambigous)
    /** \brief Copy constructor */
    SX(const SX& scalar); // copy constructor

    // Destructor
    ~SX();

  // Assignment
  SX& operator=(const SX& scalar);
  SX& operator=(double scalar); // needed since otherwise both a = SX(double) and a = Matrix(double) would be ok
//   SX& operator=(const SXMatrix& scalar);
#endif // SWIG

#ifndef SWIG
  //@{
  /** \brief  Operators that change the object */
  friend SX& operator+=(SX &ex, const SX &scalar);
  friend SX& operator-=(SX &ex, const SX &scalar);
  friend SX& operator*=(SX &ex, const SX &scalar);
  friend SX& operator/=(SX &ex, const SX &scalar);
  //@}
  
  /** \brief  Negation */
  friend SX operator-(const SX &ex);
  
  //@{
  /** \brief  Operators that create new objects (SX on the left hand side) */
  friend SX operator+(const SX &x, const SX &y);
  friend SX operator-(const SX &x, const SX &y);
  friend SX operator*(const SX &x, const SX &y);
  friend SX operator/(const SX &x, const SX &y);
  //@}
  
  //@ {
  /** \brief  Conditional operators */
  friend SX operator<=(const SX &a, const SX &b);
  friend SX operator>=(const SX &a, const SX &b);
  friend SX operator<(const SX &a, const SX &b);
  friend SX operator>(const SX &a, const SX &b);
  friend SX operator&&(const SX &a, const SX &b);
  friend SX operator||(const SX &a, const SX &b);
  friend SX operator==(const SX &a, const SX &b);
  friend SX operator!=(const SX &a, const SX &b);
  friend SX operator!(const SX &a);
  //@}
#endif // SWIG


#ifndef SWIG
  /** \brief  print to stream */
  friend std::ostream& operator<<(std::ostream &stream, const SX &scalar);

  /** \brief  string representation (SWIG workaround) */
  std::string toString() const;
  
  /** \brief  Get a pointer to the node */
  SXNode* const get() const; // note: constant pointer, not pointer to constant object! (to allow access to the counter)

  /** \brief  Access functions of the node */
  const SXNode* operator->() const;
  SXNode* operator->();
#endif // SWIG
  
  /** \brief  Perform operations by ID */
  static SX binary(int op, const SX& x, const SX& y);
  static SX unary(int op, const SX& x);

  bool isConstant() const;
  bool isInteger() const;
  bool isSymbolic() const;
  bool isBinary() const;
  bool isZero() const;
  bool isOne() const;
  bool isMinusOne() const;
  bool isNan() const;
  bool isInf() const;
  bool isMinusInf() const;
  const std::string& getName() const;
  int getOp() const;
  bool isEqual(const SX& scalar) const;
  double getValue() const;
  int getIntValue() const;

  protected:
#ifndef SWIG
SXNode* node;

/** \brief  Function that do not have corresponding c-functions and are therefore not available publically */
friend SX sign(const SX &x);
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

//  all binary operations with a particular right argument
#define binops(T,t) \
T __add__(t b){  return *$self + b;} \
T __radd__(t b){ return b + *$self;} \
T __sub__(t b){  return *$self - b;} \
T __rsub__(t b){ return b - *$self;} \
T __mul__(t b){  return *$self * b;} \
T __rmul__(t b){ return b * *$self;} \
T __div__(t b){  return *$self / b;} \
T __rdiv__(t b){ return b / *$self;} \
T __pow__(t b){  return std::pow(*$self,b);} \
T __rpow__(t b){ return std::pow(b,*$self);} \
T fmin(t b){     return std::fmin(*$self,b);} \
T fmax(t b){     return std::fmax(*$self,b);}

// Binary operations with all right hand sides
binops(SX, const SX&)
binops(SX, double)

// all unary operations
#define unops(T) \
T __neg__(){ return - *$self;}\
T exp(){ return std::exp(*$self);}\
T log(){ return std::log(*$self);}\
T sqrt(){ return std::sqrt(*$self);}\
T sin(){ return std::sin(*$self);}\
T cos(){ return std::cos(*$self);}\
T tan(){ return std::tan(*$self);}\
T arcsin(){ return std::asin(*$self);}\
T arccos(){ return std::acos(*$self);}\
T arctan(){ return std::atan(*$self);}\
T floor(){ return std::floor(*$self);}\
T ceil(){ return std::ceil(*$self);}\
T erf(){ return std::erf(*$self);}

unops(SX)

}
#endif // SWIG

#ifndef SWIG
// Template specialization
template<>
class casadi_limits<SX>{
  public:
    static bool isZero(const SX& val);
    static bool isOne(const SX& val);
    static bool isConstant(const SX& val);
    static bool isInteger(const SX& val);

    static const SX zero;
    static const SX one;
    static const SX two;
    static const SX minus_one;
    static const SX nan;
    static const SX inf; 
    static const SX minus_inf;
};

#endif // SWIG


} // namespace CasADi

#ifndef SWIG

// Template specialization
namespace std{
template<>
class numeric_limits<CasADi::SX>{
  public:
    static CasADi::SX infinity() throw();
    static CasADi::SX quiet_NaN() throw();
    // More to come
};


/** \brief  Global functions with c equivalents: The implementation and syntax mirrors the standard c functions in math.h */
#define SX CasADi::SX
SX sqrt(const SX &x);
SX sin(const SX &x);
SX cos(const SX &x);
SX tan(const SX &x);
SX atan(const SX &x);
SX asin(const SX &x);
SX acos(const SX &x);
SX exp(const SX &x);
SX log(const SX &x);
SX pow(const SX &x, const SX &n);
SX abs(const SX &x);
SX fabs(const SX &x); // same as abs
SX floor(const SX &x);
SX ceil(const SX &x);
SX erf(const SX &x);
SX fmin(const SX &a, const SX &b);
SX fmax(const SX &a, const SX &b);
#undef SX
} // namespace std

/** \brief  The following functions needs the class so they cannot be included in the beginning of the header */
#include "sx_node.hpp"

#endif // SWIG


#endif // SX_HPP
