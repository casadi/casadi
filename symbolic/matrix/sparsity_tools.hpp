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

#ifndef SPARSITY_TOOLS_HPP
#define SPARSITY_TOOLS_HPP

#include "crs_sparsity.hpp"

namespace CasADi{

  /** \brief Hash value of an integer */
  template<typename T>
  inline size_t hash_value(T v){ return size_t(v);}

  /** \brief Generate a hash value incrementally (function taken from boost) */
  template<typename T>
  inline void hash_combine(std::size_t& seed, T v){
    seed ^= hash_value(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }

  /** \brief Generate a hash value incrementally (function taken from boost) */
  inline void hash_combine(std::size_t& seed, const std::vector<int>& v){
    for(std::vector<int>::const_iterator i=v.begin(); i!=v.end(); ++i) hash_combine(seed,*i);
  }

  /** \brief Hash a sparsity pattern */
  std::size_t hash_sparsity(int nrow, int ncol, const std::vector<int>& col, const std::vector<int>& rowind);

  /**
     \brief Create a dense rectangular sparsity pattern
  **/
  CRSSparsity sp_dense(int n, int m=1);
  

  /**
     \brief Create a dense rectangular sparsity pattern
  **/
  CRSSparsity sp_dense(const std::pair<int,int> &nm );
  

  /**
     \brief Create a sparse rectangular sparsity pattern
  **/
  CRSSparsity sp_sparse(int n, int m=1);
  
  /**
     \brief Create a dense rectangular sparsity pattern
  **/
  CRSSparsity sp_sparse(const std::pair<int,int> &nm);
  
  /**
     \brief Create the sparsity pattern for a unit vector of length n and a nonzero on position el
  **/
  CRSSparsity sp_unit(int n, int el);

  /**
     \brief Create a lower triangular square sparsity pattern
     
     \see lowerSparsity
  **/
  CRSSparsity sp_tril(int n);

  /**
     \brief Create diagonal square sparsity pattern
  **/
  CRSSparsity sp_diag(int n);
  
  /**
     \brief Create a single band in a square sparsity pattern
     *
     * sp_band(n,0) is equivalent to sp_diag(n) \n
     * sp_band(n,-1) has a band below the diagonal \n
     * \param p indicate
     **/
  CRSSparsity sp_band(int n, int p);
  
  /**
     \brief Create banded square sparsity pattern
     *
     * sp_band(n,0) is equivalent to sp_diag(n) \n
     * sp_band(n,1) is tri-diagonal matrix \n
     **/
  CRSSparsity sp_banded(int n, int p);
  
  /** \brief Construct a block sparsity pattern from (row,col) vectors
      
   */
  CRSSparsity sp_rowcol(const std::vector<int>& row, const std::vector<int>& col,int nrow, int ncol);
  
  
  /** \brief Get the indices of all non-zero elements as they would appear in a Dense matrix  
      A : DenseMatrix  4 x 3
      B : SparseMatrix 4 x 3 , 5 structural non-zeros
      
      k = getNZDense(A)
      A[k] will contain the elements of A that are non-zero in B         
  */
  std::vector<int> getNZDense(const CRSSparsity& sp);
  
  
  CRSSparsity reshape(const CRSSparsity& a, int n, int m);
  
  CRSSparsity vec(const CRSSparsity& a);
  
  /** \ brief Return the transpose of the sparsity pattern
   */
  CRSSparsity trans(const CRSSparsity& a);
  
  /**
   * \brief Return the lower part of the sparsity pattern
   * 
   * \param includeDiagonal specify wether the diagonal must be part of the result
   *
   * \see sp_tril
   */
  CRSSparsity lowerSparsity(const CRSSparsity& a, bool includeDiagonal = true);
  
  /**
   * \brief Return the non-zero entries that make up the lower part of the provided matrix
   */
  std::vector<int> lowerNZ(const CRSSparsity& a);
  
  /**
     \brief Create a sparsity pattern given the nonzeros in sparse triplet form
  **/
  CRSSparsity sp_triplet(int n, int m, const std::vector<int>& row, const std::vector<int>& col, std::vector<int>& mapping, bool invert_mapping=false);
  
  /**
     \brief Create a sparsity pattern given the nonzeros in sparse triplet form (no nonzero mapping)
     columns_are_sorted==true means that the column entries already in increasing order for each row and without any duplicates
  **/
  CRSSparsity sp_triplet(int n, int m, const std::vector<int>& row, const std::vector<int>& col);
  
  
  /** \brief Get the sparsity resulting from a matrix multiplication
   */
  CRSSparsity mul(const  CRSSparsity& a, const  CRSSparsity &b);
  
  /** \brief Concatenate a list of sparsities vertically
  * Alternative terminology: vertical stack, vstack, vertical append, [a;b]
  */
  CRSSparsity vertcat(const std::vector<CRSSparsity > &v);

  /** \brief Concatenate a list of sparsities horizontally
  * Alternative terminology: horizontal stack, hstack, horizontal append, [a b]
  */
  CRSSparsity horzcat(const std::vector<CRSSparsity > &v);

  /** \brief   Construct a Sparsity with given blocks on the diagonal */
  CRSSparsity blkdiag(const std::vector< CRSSparsity > &v);

  #ifndef SWIG
  CRSSparsity vertcat(const CRSSparsity &x, const CRSSparsity &y);

  CRSSparsity horzcat(const CRSSparsity &x, const CRSSparsity &y);
  
  CRSSparsity blkdiag(const CRSSparsity &x, const CRSSparsity &y);
  #endif // SWIG
  
  /** \brief Represent a sparsity pattern as an array of integers, the most compact way of representing a sparsity pattern
      The format:
      * The first two entries are the number of rows (nrow) and columns (ncol)
      * The next nrow+1 entries are the row offsets (rowind). Note that the last element rowind[nrow] gives the number of nonzeros
      * The last rowind[nrow] entries are the column indices
      **/
  /// @{
  
  /// Compress a sparsity pattern
  std::vector<int> sp_compress(const CRSSparsity& a);
  
  /// Decompress a sparsity pattern
  CRSSparsity sp_compress(const std::vector<int>& v);
  
#ifndef SWIG
  /// Decompress a sparsity pattern (array version)
  CRSSparsity sp_compress(const int* v);
#endif // SWIG  

  /// Obtain the structural rank of a sparsity-pattern
  int rank(const CRSSparsity& a);
  
  /// Check whether the sparsity-pattern inidcates structural singularity
  bool isSingular(const CRSSparsity& a);

  /// @}


}

#endif // SPARSITY_TOOLS_HPP
