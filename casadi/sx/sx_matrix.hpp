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

#ifndef SX_MATRIX_HPP
#define SX_MATRIX_HPP

#include "sx.hpp"
#include <vector>
#include <string>
#include <algorithm>
#include "../matrix/matrix.hpp"

namespace CasADi{
  
/** \brief General sparse symbolic matrix class
  NOTE NOTE NOTE
  MARKED FOR REMOVAL

  General sparse symbolic matrix class that is designed with the idea that "everything is a matrix", that is, also scalars and vectors.
  This philosophy makes it easy to use and to interface in particularily with Matlab and Python.
  It is likely not as efficient as expression template based matrix classes, since this is not the objective,
  but as symbolic evaluations only need to take place in general (not repeatedly as for numerical evaluations),
  it is very reasonable to trade speed for ease-to-use.
  
  The syntax tries to stay as close as possible to both the SX syntax (which is essentially equivalent to a double), as well ublas syntax 
  when it comes to vector/matrix operations. Vectors are considered to be column vectors, and index starts with 0.
  The syntax is also identical to that of the matrix expression class, MX, a more general class where syntax trees of scalar and matrix
  operations (including branchings) is built.

  It is possible to make the class more efficeint, e.g. by copy on assignment, though if it is motivated is unclear.
  A more reasonable approach in these situations is to use the MX class to build up a syntax tree, and then evaluate symbolically.
  
  The storage format is a (modified) compressed row storage (CRS) format. This way, a vector element can always be accessed in constant time.
  
  \author Joel Andersson 
  \date 2010

  joel.andersson@esat.kuleuven.be
*/

class SXMatrix : public Matrix<SX>{
public:

/** \brief  constructors */
/// empty 0-by-0 matrix constructor
SXMatrix();

/// empty n-by-m matrix constructor
SXMatrix(int n, int m);

/// dense n-by-m matrix filled with val constructor
SXMatrix(int n, int m, const SX& val);    

//@{
/** \brief  These constructors enable implicit type conversion */
SXMatrix(const SX &scalar);      // create a variable from a scalar
SXMatrix(double val);            // constant
//@}

//@{
/** \brief  Create an expression from an stl vector (can be implemented as templates) */
SXMatrix(const std::vector<SX>& x);
SXMatrix(const std::vector<double>& x);
//@}

//@{
/** \brief  Create a non-vector expression from an stl vector */
SXMatrix(const std::vector<SX>& x,  int n, int m);
SXMatrix(const std::vector<double>& x,  int n, int m);
//@}

/** \brief  Create a matrix of symbolic variables  */
explicit SXMatrix(const std::string& name, int n=1, int m=1);

//@{
/** \brief  Operators that changes the object */
friend SXMatrix& operator+=(SXMatrix &ex, const SXMatrix &expr);
friend SXMatrix& operator-=(SXMatrix &ex, const SXMatrix &expr);
friend SXMatrix& operator*=(SXMatrix &ex, const SXMatrix &expr);
friend SXMatrix& operator/=(SXMatrix &ex, const SXMatrix &expr);
//@}

/** \brief  Negation */
friend SXMatrix operator-(SXMatrix &ex);

/** \brief  Operators that create new objects */
friend SXMatrix operator+(const SXMatrix &x, const SXMatrix &y);
friend SXMatrix operator-(const SXMatrix &x, const SXMatrix &y);

/** \brief Element-wise multiplication */
friend SXMatrix operator*(const SXMatrix &x, const SXMatrix &y);
friend SXMatrix operator/(const SXMatrix &x, const SXMatrix &y);

friend std::ostream& operator<<(std::ostream &stream, const SXMatrix &mat);

};

} // namespace CasADi

#endif // SX_MATRIX_HPP


