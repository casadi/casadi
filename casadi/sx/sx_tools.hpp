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

#ifndef SX_TOOLS_HPP
#define SX_TOOLS_HPP

#include "sx.hpp"
#include "../matrix/matrix_tools.hpp"

#ifdef WITH_UBLAS
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
namespace ublas = boost::numeric::ublas;
#endif

/** Functions using SX */

namespace CasADi{

#ifndef SWIG
/** \brief Make a vector/matrix of symbolic variables - dimension 0 */
void make_symbolic(SX& v, const std::string& name);

/** \brief Make a vector/matrix of symbolic variables - higher dimension recursively */
template<typename A>
void make_symbolic(std::vector< A >& v, const std::string& name){
  for(unsigned int i=0; i<v.size(); ++i){
    std::stringstream ss;
    ss << name << "_" << i;
    make_symbolic(v[i],ss.str());
  }
}
#endif

/** \brief Create an n-by-m matrix with symbolic variables */
Matrix<SX> symbolic(const std::string& name, int n=1, int m=1);

/** \brief Create a vector of length p with n-by-m matrices with symbolic variables */
std::vector<Matrix<SX> > symbolic(const std::string& name, int n, int m, int p);

/** \brief Create an matrix with symbolic variables, given a sparsity pattern */
Matrix<SX> symbolic(const std::string& name, const CRSSparsity& sp);

/** \brief Create a one-dimensional stl vector of length n with symbolic variables */
std::vector<SX> create_symbolic(const std::string& name, int n);

/** \brief Create a two-dimensional stl vector of length n-by-m with symbolic variables */
std::vector< std::vector<SX> > create_symbolic(const std::string& name, int n, int m);

/** \brief Create a three-dimensional stl vector of length n-by-m-by-p with symbolic variables */
std::vector< std::vector< std::vector< SX> > > create_symbolic(const std::string& name, int n, int m, int p);
  
/** \brief  Expand the expression as a weighted sum (with constant weights)  */
void expand(const Matrix<SX>& ex, Matrix<SX> &weights, Matrix<SX>& terms);

/** \brief  Simplify the expression: formulates the expression as and eliminates terms */
void simplify(SX& ex);

#ifndef SWIG
/** \brief Make a vector/matrix of symbolic variables */
template <typename iter_type, typename func_type>
void make_symbolic(iter_type first, func_type last, const std::string& name){
  int count = 0;
  for(iter_type cur=first; cur!=last; ++cur){
    std::stringstream ss;
    ss << name << "_" << count++;
    *cur = SX(ss.str());
  }
}
#endif

/** \brief  Concatenate */
Matrix<SX>& operator<<(Matrix<SX>& expr, const Matrix<SX>& add); // remove when C++0X becomes available

/** \brief  create an n-by-n identity matrix */
Matrix<SX> eyeSX(int n); 

/** \brief  create a matrix with all infinities */
Matrix<SX> infSX(int n=1,int m=1);

/** \brief  create a matrix with all ones */
Matrix<SX> onesSX(int n, int m=1);

/** \brief  create a matrix with all zeros */
Matrix<SX> zerosSX(int n, int m=1);

/** \brief Create a piecewise constant function 
  Create a piecewise constant function with n=val.size() intervals

  Inputs:
  \param t a scalar variable (e.g. time)
  \param tval vector with the discrete values of t at the interval transitions (length n-1)
  \param val vector with the value of the function for each interval (length n)
*/
Matrix<SX> pw_const(const Matrix<SX> &t, const Matrix<SX> &tval, const Matrix<SX> &val);

/** Create a piecewise linear function 
  Create a piecewise linear function:

  Inputs:
  \brief t a scalar variable (e.g. time)
  \brief tval vector with the the discrete values of t (monotonically increasing)
  \brief val vector with the corresponding function values (same length as tval)
*/
Matrix<SX> pw_lin(const SX &t, const Matrix<SX> &tval, const Matrix<SX> &val);

Matrix<SX> if_else(const Matrix<SX> &cond, const Matrix<SX> &if_true, const Matrix<SX> &if_false);
/// Heaviside step function
Matrix<SX> heaviside(const Matrix<SX> &x); // heaviside step function

/// sign function
Matrix<SX> sign(const Matrix<SX> &x);     // sign function

/** \brief  Integrate f from a to b using Gaussian quadrature with n points */
Matrix<SX> gauss_quadrature(Matrix<SX> f, const Matrix<SX> &x, const Matrix<SX> &a, const Matrix<SX> &b, int order=5, const Matrix<SX>& w=Matrix<SX>());

#ifndef SWIG
/**
Returns true if at least one element in list contains the scalar e.
*/
bool contains(const Matrix<SX> &list, const Matrix<SX> &e);
#endif

/** \brief  Simplify an expression */
void simplify(Matrix<SX> &ex);
/// remove identical calculations
void compress(Matrix<SX> &ex, int level=5); 
/// substitute variable var with expression expr
Matrix<SX> substitute(const Matrix<SX> &ex, const Matrix<SX> &var, const Matrix<SX> &expr); 

// /** \brief  Make the expression smooth by replacing non-smooth nodes with binary variables */
//void makeSmooth(Matrix<SX> &ex, Matrix<SX> &bvar, Matrix<SX> &bexpr);

/** \brief  Substitute derivatives with variables */
/** \brief void replaceDerivatives(Matrix<SX> &ex, const Matrix<SX> &var, const Matrix<SX> &dvar); */

#ifndef SWIG
// "operator?:" can not be overloaded
template<typename T>
T if_else(const SX& cond, const T& if_true, const T &if_false){
  return if_false + (if_true-if_false)*cond;
}
#endif

/** \brief  Get the sparsity pattern of a matrix */
Matrix<SX> spy(const Matrix<SX>& A);

#ifndef SWIG
/** \brief  Create a block matrix */
template<int n, int m>
Matrix<SX> blockmatrix(Matrix<SX> array[n][m]){
/** \brief  Return matrix */
  Matrix<SX> ret;

/** \brief  loop over rows */
  for(int i=0; i<n; ++i){
/** \brief  Create a row */
    Matrix<SX> row;
    
/** \brief  append components to the row */
    for(int j=0; j<m; ++j){
      row << array[i][j];
    }
    
/** \brief  append row to matrix */
    ret << trans(row);
  }

  return ret;
}

/** \brief  Create a block matrix (vector) */
template<int n>
Matrix<SX> blockmatrix(Matrix<SX> array[n]){
/** \brief  Return matrix */
  Matrix<SX> ret;

/** \brief  loop over rows */
  for(int i=0; i<n; ++i){
/** \brief  append components */
    ret << array[i];
  }

  return ret;
}
#endif

/** \brief  Get the number of nodes of the tree of a Matrix<SX> */
int numNodes(const Matrix<SX>& A);

/// Check dependency: very inefficient algorithm
bool dependsOn(const Matrix<SX>& f, const Matrix<SX> &arg);

/** \brief  check if smooth */
bool isSmooth(const Matrix<SX>& ex);

/** \brief  check if symbolic */
bool isSymbolic(const Matrix<SX>& ex);

//@{
/** \brief Calculate jacobian via source code transformation

Uses CasADi::SXFunction::jac
 */
Matrix<SX> jacobian(const Matrix<SX> &ex, const Matrix<SX> &arg);
Matrix<SX> gradient(const Matrix<SX> &ex, const Matrix<SX> &arg);
Matrix<SX> hessian(const Matrix<SX> &ex, const Matrix<SX> &arg);
void hessian(const Matrix<SX> &ex, const Matrix<SX> &arg, Matrix<SX> &H, Matrix<SX> &g); // hessian and gradient
//@}

/** \brief  Obtain the values of a constant expression */
double getValue(const Matrix<SX> &ex, int i=0, int j=0);          // for constant expressions only
int getIntValue(const Matrix<SX> &ex, int i=0, int j=0);          // integer version
void getValue(const Matrix<SX> &ex, double *res); // for all constant expressions
void getIntValue(const Matrix<SX> &ex, int *res); // integer version
const std::string& getName(const Matrix<SX> &ex); // get the name (only for scalar variables)

/** \brief  To and from string */
std::istream& operator>>(std::istream &stream, Matrix<SX> &expr);

Matrix<SX> operator<=(const Matrix<SX> &a, const Matrix<SX> &b);
Matrix<SX> operator>=(const Matrix<SX> &a, const Matrix<SX> &b);
Matrix<SX> operator<(const Matrix<SX> &a, const Matrix<SX> &b);
Matrix<SX> operator>(const Matrix<SX> &a, const Matrix<SX> &b);
#ifndef SWIG
Matrix<SX> operator&&(const Matrix<SX> &a, const Matrix<SX> &b);
Matrix<SX> operator||(const Matrix<SX> &a, const Matrix<SX> &b);
Matrix<SX> operator!(const Matrix<SX> &a);
#endif
Matrix<SX> operator==(const Matrix<SX> &a, const Matrix<SX> &b);
Matrix<SX> operator!=(const Matrix<SX> &a, const Matrix<SX> &b);

/** \brief  Fill the matrix with the value val, make empty sparse if zero */
void fill(Matrix<SX>& mat, const SX& val);

} // namespace CasADi

#endif // SX_TOOLS_HPP
