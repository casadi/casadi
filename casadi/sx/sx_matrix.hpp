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
SXMatrix();                               // 
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
explicit SXMatrix(const std::string& name, int n=1, int m=1);   // variable

/** \brief  destructor */
~SXMatrix();

/** Element class for SXmatrix 
    Element is the return type for operator[] and operator() of the SXMatrix class.
    From several alternative solutions, this was the only one which didn't cause ambigousity.
    Suggestions for improvement are welcome.
  \author Joel Andersson 
  \date 2010

*/
class Element{
  friend class SXMatrix;
  public:
    //@{
    /// Methods that modify a part of the parent obejct (A[i] = ?, A[i] += ?, etc.)
    SX& operator=(const Element &y);
    SX& operator=(const SXMatrix &y);
    void operator-=(const SXMatrix &y){*this = *this - y;}
    void operator+=(const SXMatrix &y){*this = *this + y;}
    void operator*=(const SXMatrix &y){*this = *this * y;}
    void operator/=(const SXMatrix &y){*this = *this / y;}
    //@}
    
//    operator SXMatrix(){ return mat.getElement(i,j);}
    
    /// Get a pointer to the node
    SXNode* const get() const;

    /// Access functions of the node
    const SXNode* operator->() const;
    SXNode* operator->();

    private:
    Element(SXMatrix& mat, int i, int j);
    SXMatrix& mat;
    int i, j;
};

  /// Implicit type conversion from Element to SXMatrix
  SXMatrix(const Element &el);

  /// Access an element 
  SXMatrix::Element operator()(int i, int j=0);             // matrix style access 

  /// Const access an element 
  const SX operator()(int i, int j=0) const; // matrix style access

  /// Const access a non-zero element
  const SX& operator[](int k) const;
  
  /// Access a non-zero element
  SX& operator[](int k);

  /// Get several non-zero elements
  const std::vector<SX> operator[](const std::vector<int>& ind) const;

  /// A submatrix class enabling the syntax "A({1,3},{2,4}) = B;"
  class Sub{
    friend class SXMatrix;
    public:
      void operator=(const Element &y);
      void operator=(const SXMatrix &y);

      operator const SXMatrix() const;
      operator SXMatrix& ();
    
      private:
      Sub(SXMatrix& mat, std::vector<int> i, std::vector<int> j);
      SXMatrix& mat;
      std::vector<int> i, j;
    };

  /// Get a const submatrix (for B = A(I,J))
  const SXMatrix operator()(const std::vector<int>& i, const std::vector<int>& j=std::vector<int>(1,0)) const;

  /// Get a reference to a submatrix (for A(I,J) = B and B = A(I,J))
  Sub operator()(const std::vector<int>& i, const std::vector<int>& j=std::vector<int>(1,0));

void clear();
/// Resize/reshape an SXMatrix in-place
void resize(int n, int m);
void reserve(int nnz);


//@{
/** \brief  Operators that changes the object */
// They do not, in fact change the object do they?
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

//@{
  /** \brief  Perform operations by ID */
friend SXMatrix binary(int op, const SXMatrix &x, const SXMatrix &y);
friend SXMatrix unary(int op, const SXMatrix &x);
friend SXMatrix scalar_matrix(int op, const SX &x, const SXMatrix &y);
friend SXMatrix matrix_scalar(int op, const SXMatrix &x, const SX &y);
friend SXMatrix matrix_matrix(int op, const SXMatrix &x, const SXMatrix &y);
//@}

friend std::ostream& operator<<(std::ostream &stream, const SXMatrix &mat);

/** \brief  Fill the matrix with the value val, make empty sparse if zero */
void fill(const SX& val);

/** \brief  Not implemented */
void push_back(const SX& el); // Add element at the end (public member function)
void pop_back();              // Delete last element (public member function)

};

} // namespace CasADi


/** \brief  Global functions with c equivalents: The implementation and syntax mirrors the standard c functions in math.h */
namespace std{
#define SXMatrix CasADi::SXMatrix
SXMatrix sin(const SXMatrix &x);
SXMatrix cos(const SXMatrix &x);
SXMatrix tan(const SXMatrix &x);
SXMatrix atan(const SXMatrix &x);
SXMatrix asin(const SXMatrix &x);
SXMatrix acos(const SXMatrix &x);
SXMatrix exp(const SXMatrix &x); // natural expontial
SXMatrix log(const SXMatrix &x); // natuaral logarithm
SXMatrix pow(const SXMatrix &x, const SXMatrix &n); // power function
SXMatrix sqrt(const SXMatrix &x); // square root
SXMatrix fmin(const SXMatrix &x, const SXMatrix &y);
SXMatrix fmax(const SXMatrix &x, const SXMatrix &y);
SXMatrix floor(const SXMatrix &x);
SXMatrix ceil(const SXMatrix &x); 
SXMatrix erf(const SXMatrix &x); 
#undef SXMatrix
} // namespace std


#endif // SX_MATRIX_HPP


