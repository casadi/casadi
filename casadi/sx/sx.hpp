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
  
  //@{
  /** \brief  constant nodes (to make sure that there are no dublications of them) */
  static const SX zero; // constant 0
  static const SX one;  // constant 1
  static const SX two;  // constant 2
  static const SX mone;  // constant -1
  static const SX nan; // not a number
  static const SX inf; // infinity
  static const SX minf; // minus infinity
  //@}
  
  /** \brief  print to stream */
  friend std::ostream& operator<<(std::ostream &stream, const SX &scalar);

  /** \brief  string representation (SWIG workaround) */
  std::string toString() const;
  
  /** \brief  Get a pointer to the node */
  SXNode* const get() const; // note: constant pointer, not pointer to constant object! (to allow access to the counter)

  /** \brief  Access functions of the node */
  const SXNode* operator->() const;
  SXNode* operator->();

  /** \brief  Perform operations by ID */
  static SX binary(int op, const SX& x, const SX& y);
  static SX unary(int op, const SX& x);

  protected:
SXNode* node;

/** \brief  Function that do not have corresponding c-functions and are therefore not available publically */
friend SX sign(const SX &x);
/** \brief inline if-test */
friend SX if_else(const SX& cond, const SX& if_true, const SX& if_false); // replaces the ternary conditional operator "?:", which cannot be overloaded

};

/** \brief Make a vector/matrix of symbolic variables - dimension 0 */
void make_symbolic(SX& v, const std::string& name);

/** \brief Make a vector/matrix of symbolic variables - higher dimension recursively */
template<typename A>
void make_symbolic(std::vector< A >& v, const std::string& name){
  for(int i=0; i<v.size(); ++i){
    std::stringstream ss;
    ss << name << "_" << i;
    make_symbolic(v[i],ss.str());
  }
}

/** \brief Create a one-dimensional stl vector of length n with symbolic variables */
std::vector<SX> create_symbolic(const std::string& name, int n);

/** \brief Create a two-dimensional stl vector of length n-by-m with symbolic variables */
std::vector< std::vector<SX> > create_symbolic(const std::string& name, int n, int m);

/** \brief Create a three-dimensional stl vector of length n-by-m-by-p with symbolic variables */
std::vector< std::vector< std::vector< SX> > > create_symbolic(const std::string& name, int n, int m, int p);

} // namespace CasADi

/** \brief  Global functions with c equivalents: The implementation and syntax mirrors the standard c functions in math.h */
namespace std{
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


#endif // SX_HPP
