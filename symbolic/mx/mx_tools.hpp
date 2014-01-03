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

#ifndef MX_TOOLS_HPP
#define MX_TOOLS_HPP

#include "mx.hpp"

#include "../matrix/generic_matrix_tools.hpp"
#include "../matrix/generic_expression_tools.hpp"
#include "../fx/linear_solver.hpp"

namespace CasADi{

  /** \brief  concatenate vertically
  *
  *  vertcat(vertsplit(x,...)) = x
  */
  MX vertcat(const std::vector<MX>& x);

  /** \brief  split vertically, retaining groups of rows
  * \param output_offset List of all start rows for each group
  *      the last row group will run to the end.
  * 
  *   vertcat(vertsplit(x,...)) = x
  */
  std::vector<MX> vertsplit(const MX& x, const std::vector<int>& output_offset);

  /** \brief  split vertically, retaining fixed-sized groups of rows
  * \param incr Size of each group of rows
  *
  *  vertcat(vertsplit(x,...)) = x
  */
  std::vector<MX> vertsplit(const MX& x, int incr=1);
  
  /** \brief  concatenate horizontally 
  *
  *   horzcat(horzsplit(x,...)) = x
  */
  MX horzcat(const std::vector<MX>& comp);
  
  /** \brief  split horizontally, retaining groups of columns
  * \param output_offset List of all start columns for each group
  *      the last column group will run to the end.
  *
  *   horzcat(horzsplit(x,...)) = x
  */
  std::vector<MX> horzsplit(const MX& x, const std::vector<int>& output_offset);

  /** \brief  split horizontally, retaining fixed-sized groups of columns
  * \param incr Size of each group of columns
  *
  *   horzcat(horzsplit(x,...)) = x
  */
  std::vector<MX> horzsplit(const MX& x, int incr=1);
  
  /** \brief Construct a matrix from a list of list of blocks.
  *
  *   blockcat(blocksplit(x,...,...)) = x
  */
  MX blockcat(const std::vector< std::vector<MX > > &v);
  
  /** \brief  chop up into blocks
  * \brief vert_offset Defines the boundaries of the block rows
  * \brief horz_offset Defines the boundaries of the block columns
  *
  *   blockcat(blocksplit(x,...,...)) = x
  */
  std::vector< std::vector<MX > > blocksplit(const MX& x, const std::vector<int>& vert_offset, const std::vector<int>& horz_offset);

  /** \brief  chop up into blocks
  * \brief vert_incr Defines the increment for block boundaries in row dimension
  * \brief horz_incr Defines the increment for block boundaries in column dimension
  *
  *   blockcat(blocksplit(x,...,...)) = x
  */
  std::vector< std::vector<MX > > blocksplit(const MX& x, int vert_incr = 1, int horz_incr = 1);

#ifndef SWIG
  /** \brief Construct a matrix from a list of list of blocks.*/
  MX blockcat(const MX &A,const MX &B,const MX &C,const MX &D);
#endif // SWIG

  /** \brief  concatenate vertically while vectorizing all arguments with vec */
  MX veccat(const std::vector<MX>& comp);

  /** \brief  concatenate vertically while vectorizing all arguments with vecNZ */
  MX vecNZcat(const std::vector<MX>& comp);

#ifndef SWIG
  /** \brief  concatenate vertically, two matrices */
  MX vertcat(const MX& a, const MX& b);

  /** \brief  concatenate horizontally, two matrices */
  MX horzcat(const MX& a, const MX& b);
#endif // SWIG

  /** \brief  Frobenius norm  */
  MX norm_F(const MX &x);

  /** \brief  2-norm  */
  MX norm_2(const MX &x);

  /** \brief 1-norm  */
  MX norm_1(const MX &x);

  /** \brief Infinity-norm */
  MX norm_inf(const MX &x);

  /** \brief Transpose an expression */
  MX trans(const MX &x);

  /** \brief  Take the matrix product of 2 MX objects
  *
  * With optional sp_z you can specify the sparsity of the result
  * A typical use case might be where the product is only constructed to 
  * inspect the trace of it. sp_z diagonal will be more efficient then.
  *
  */
  MX mul(const MX &x, const MX &y, const CRSSparsity& sp_z=CRSSparsity());

  /** \brief  Take the matrix product of n MX objects */
  MX mul(const std::vector< MX > &x);

  /** \brief  Take the inner product of two vectors 
      Equals
      \code
      trans(x)*y
      \endcode
      with x and y vectors
  */
  MX inner_prod(const MX &x, const MX &y);

  /** \brief  Take the outer product of two vectors 
      Equals
      \code
      x*trans(y)
      \endcode
      with x and y vectors
  */
  MX outer_prod(const MX &x, const MX &y);

  /** \brief Branching on MX nodes
      Ternary operator, "cond ? if_true : if_false"
  */
  MX if_else(const MX &cond, const MX &if_true, const MX &if_false); 

#ifndef SWIG
  //! \brief Returns a reshaped version of the MX
  MX reshape(const MX &x, int n, int m);
#endif // SWIG

  //! \brief Returns a reshaped version of the MX, dimensions as a vector
  MX reshape(const MX &x, const std::vector<int> sz);

  //! \brief Reshape the MX
  MX reshape(const MX &x, const CRSSparsity& sp);

  /** \brief Returns a vectorized version of the MX
      Vectorizing is an expensive operation, unlike flatten
      Same as reshape(trans(x), x.numel(),1)
    
      a b
      c d 
    
      turns into
    
      a
      c
      b
      d
    
  */
  MX vec(const MX &x);

  /** \brief Returns a flattened version of the MX
      Flattening is a cheap (non-copying) operation
      Same as reshape(x, x.numel(),1)
    
      a b
      c d 
    
      turns into
    
      a
      b
      c
      d
    
  */
  MX flatten(const MX &x);

  /** \brief Returns a flattened version of the MX, preserving only nonzeros
   */
  MX vecNZ(const MX &x);

  /** \brief  Unite two matrices no overlapping sparsity */
  MX unite(const MX& A, const MX& B);

  /** \brief  check if symbolic */
  bool isSymbolic(const MX& ex);

  /** \brief  check if all nonzeros are symbolic (this function is currently identical to isSymbolic) */
  bool isSymbolicSparse(const MX& ex);

  /** \brief  check if identity */
  bool isIdentity(const MX& ex);

  /** \brief  check if zero (note that false negative answers are possible) */
  bool isZero(const MX& ex);

  /** \brief  check if zero (note that false negative answers are possible) */
  bool isOne(const MX& ex);

  /** \brief  Simplify an expression */
  void simplify(MX& ex);

  /** \brief  Is the expression a transpose? */
  bool isTranspose(const MX& ex);

  /** \brief  check if zero (note that false negative answers are possible) */
  bool isMinusOne(const MX& ex);

  /** \brief  check if vector */
  bool isVector(const MX& ex);

  /** \brief  check if vector */
  bool isDense(const MX& ex);
  
  /// Checks if expression does not contain NaN or Inf
  bool isRegular(const MX& ex);

  MX trace(const MX& A);

  /** \brief Repeat matrix A n times vertically and m times horizontally */
  MX repmat(const MX &A, int n, int m); 

  /** \brief create a clipped view into a matrix
      Create a sparse matrix from a dense matrix A, with sparsity pattern sp
  **/
  //MX clip(const MX& A, const CRSSparsity& sp);

  /** \brief  Make the matrix dense if not already*/
  void makeDense(MX& x);

  /** \brief  Make the matrix dense if not already */
  MX densify(const MX& x);

  /** \brief  Create a parent MX on which all given MX's will depend.

      In some sense, this function is the inverse of 
 
      \param deps  Must all be symbolic matrices.
  */
  MX createParent(std::vector<MX> &deps);

  /** \brief  Create a parent MX on which a bunch of MX's (sizes given as argument) will depend
   */
  std::pair<MX, std::vector<MX> > createParent(const std::vector<MX> &deps);


  /** \brief  Create a parent MX on which a bunch of MX's (sizes given as argument) will depend
   */
  std::pair<MX, std::vector<MX> > createParent(const std::vector<CRSSparsity> &deps);

  /** Count number of nodes */
  int countNodes(const MX& A);

  /** \brief  Get the diagonal of a matrix or construct a diagonal

      When the input is square, the diagonal elements are returned.
      If the input is vector-like, a diagonal matrix is constructed with it.
  */
  MX diag(const MX& x);

  /** \brief   Construct a matrix with given blocks on the diagonal */
  MX blkdiag(const std::vector<MX> &A);

#ifndef SWIG
  /** \brief   Construct a matrix with given blocks on the diagonal */
  MX blkdiag(const MX &A, const MX& B);
#endif // SWIG

  /** \brief Return a row-wise summation of elements */
  MX sumRows(const MX &x);

  /** \brief Return a column-wise summation of elements */
  MX sumCols(const MX &x);

  /// Return summation of all elements
  MX sumAll(const MX &x); 

  /** \brief  Evaluate a polynomial with coefficeints p in x */
  MX polyval(const MX& p, const MX& x);

  /** \brief  Construct symbolic arrays and variables using CasADi's MX expression graph representation
      The "msym" function is intended to work in a similar way as "sym" used in the Symbolic Toolbox for Matlab but instead creating an MX object.
      The MX expression graph is more general but also have considerably more overhead than the alternative SX expression graph.
  */
  //@{
  /** \brief Create a matrix symbolic variable of given sparsity */
  MX msym(const std::string& name, int n=1, int m=1);

  /** \brief Create a matrix symbolic variable of given sparsity */
  MX msym(const std::string& name, const std::pair<int,int> & nm);

  /** \brief Create a matrix variable from a constant matrix */
  MX msym(const Matrix<double>& x);

  /** \brief Create a matrix symbolic variable of given sparsity */
  MX msym(const std::string& name, const CRSSparsity& sp);

  /** \brief Create a vector of length p with with matrix symbolic variables of given sparsity */
  std::vector<MX> msym(const std::string& name, const CRSSparsity& sp, int p);

  /** \brief Create a vector of length p with n-by-m matrix symbolic variables */
  std::vector<MX> msym(const std::string& name, int n, int m, int p);

  /** \brief Create a vector of length r of vectors of length p with matrix symbolic variables with given sparsity*/
  std::vector<std::vector<MX> > msym(const std::string& name, const CRSSparsity& sp, int p, int r);

  /** \brief Create a vector of length r of vectors of length p with n-by-m matrices with symbolic variables */
  std::vector<std::vector<MX> > msym(const std::string& name, int n, int m, int p, int r);

  //@}

  /** \brief  Check if two expressions are equal
  *
  *  Might very well give false negatives
  *
  *   Note: does not work when CasadiOptions.setSimplificationOnTheFly(False) was called
   */
  bool isEqual(const MX& ex1,const MX &ex2);

  /** \brief Get a string representation for a binary MX, using custom arguments */
  std::string getOperatorRepresentation(const MX& x, const std::vector<std::string>& args);

  /** \brief  Substitute variable v with expression vdef in an expression ex */
  MX substitute(const MX &ex, const MX& v, const MX& vdef);

  /** \brief  Substitute variable var with expression expr in multiple expressions */
  std::vector<MX> substitute(const std::vector<MX> &ex, const std::vector<MX> &v, const std::vector<MX> &vdef);
  
  /** \brief  Substitute variable v with expression vdef in an expression ex, preserving nodes */
  MX graph_substitute(const MX &ex, const std::vector<MX> &v, const std::vector<MX> &vdef);

  /** \brief  Substitute variable var with expression expr in multiple expressions, preserving nodes  */
  std::vector<MX> graph_substitute(const std::vector<MX> &ex, const std::vector<MX> &v, const std::vector<MX> &vdef);

  /** \brief Inplace substitution
   * Substitute variables v out of the expressions vdef sequentially 
   */
#ifndef SWIG
  void substituteInPlace(const std::vector<MX>& v, std::vector<MX>& vdef, bool reverse=false);
#else // SWIG
  void substituteInPlace(const std::vector<MX>& v, std::vector<MX>& INOUT, bool reverse=false);
#endif // SWIG

  /** \brief Inplace substitution with piggyback expressions
   * Substitute variables v out of the expressions vdef sequentially, as well as out of a number of other expressions piggyback */
#ifndef SWIG
  void substituteInPlace(const std::vector<MX>& v, std::vector<MX>& vdef, std::vector<MX>& ex, bool reverse=false);
#else // SWIG
  void substituteInPlace(const std::vector<MX>& v, std::vector<MX>& INOUT, std::vector<MX>& INOUT, bool reverse=false);
#endif // SWIG

  template<> inline
  MX GenericMatrix<MX>::sym(const std::string& name, const CRSSparsity& sp){ return msym(name,sp);}

  /** \brief Extract shared subexpressions from an set of expressions */
  void extractShared(std::vector<MX>& ex, 
                     std::vector<MX>& v, std::vector<MX>& vdef, 
                     const std::string& v_prefix="v_", const std::string& v_suffix="");

  /** \brief Print compact, introducing new variables for shared subexpressions */
  void printCompact(const MX& ex, std::ostream &stream=std::cout);

  //@{
  /** \brief Calculate jacobian via source code transformation

      Uses CasADi::MXFunction::jac
  */
  MX jacobian(const MX &ex, const MX &arg);
  MX gradient(const MX &ex, const MX &arg);
  MX tangent(const MX &ex, const MX &arg);
  //@}

  /** \brief Matrix determinant (experimental) */
  MX det(const MX& A);

  /** \brief Matrix inverse (experimental) */
  MX inv(const MX& A);

  /** \brief Get all symbols contained in the supplied expression
  * Get all symbols on which the supplied expression depends
  * \see MXFunction::getFree()
  */
  std::vector<MX> getSymbols(const MX& e);
  
  /** \brief Get all symbols contained in the supplied expression
  * Get all symbols on which the supplied expression depends
  * \see MXFunction::getFree()
  */
  std::vector<MX> getSymbols(const std::vector<MX>& e);
  
  /** \brief Expand MX graph to SXFunction call
  *
  *  Expand the given expression e, optionally
  *  supplying expressions contained in it at which expansion should stop.
  *
  */
  MX matrix_expand(const MX& e, const std::vector<MX> &boundary = std::vector<MX>());

  /** \brief Expand MX graph to SXFunction call
  *
  *  Expand the given expression e, optionally
  *  supplying expressions contained in it at which expansion should stop.
  *
  */
  std::vector<MX> matrix_expand(const std::vector<MX>& e, const std::vector<MX> &boundary = std::vector<MX>());

  /** \brief Kronecker tensor product
  *
  * Creates a block matrix in which each element (i,j) is a_ij*b 
  */
  MX kron(const MX& a, const MX& b);

  /** \brief Solve a system of equations: A*x = b 
  */
  MX solve(const MX& A, const MX& b, linearSolverCreator lsolver, const Dictionary& dict = Dictionary());


} // namespace CasADi

#ifdef SWIG
// Template instantiations
%template(Pair_MX_MXVector) std::pair<CasADi::MX, std::vector<CasADi::MX> >;
#endif // SWIG





#endif // MX_TOOLS_HPP
