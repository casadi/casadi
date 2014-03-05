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

#include "sparsity.hpp"

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
  std::size_t hash_sparsity(int nrow, int ncol, const std::vector<int>& colind, const std::vector<int>& row);

  /**
     \brief Create a dense rectangular sparsity pattern
  **/
  Sparsity sp_dense(int nrow, int ncol=1);
  
  /**
     \brief Create a dense sparsity pattern
  **/
  Sparsity sp_dense(const std::pair<int,int> &nm );

  /**
     \brief Create a sparse sparsity pattern
  **/
  Sparsity sp_sparse(int nrow, int ncol=1);
  
  /**
     \brief Create a dense sparsity pattern
  **/
  Sparsity sp_sparse(const std::pair<int,int> &nm);
  
  /**
     \brief Create the sparsity pattern for a unit vector of length n and a nonzero on position el
  **/
  Sparsity sp_unit(int n, int el);

  /**
     \brief Create a upper triangular square sparsity pattern
     
     \see upperSparsity
  **/
  Sparsity sp_triu(int n);

  /**
     \brief Create a lower triangular square sparsity pattern
     
     \see lowerSparsity
  **/
  Sparsity sp_tril(int n);

  /**
     \brief Create diagonal square sparsity pattern
  **/
  Sparsity sp_diag(int n);
  
  /**
     \brief Create a single band in a square sparsity pattern
     *
     * sp_band(n,0) is equivalent to sp_diag(n) \n
     * sp_band(n,-1) has a band below the diagonal \n
     * \param p indicate
     **/
  Sparsity sp_band(int n, int p);
  
  /**
     \brief Create banded square sparsity pattern
     *
     * sp_band(n,0) is equivalent to sp_diag(n) \n
     * sp_band(n,1) is tri-diagonal matrix \n
     **/
  Sparsity sp_banded(int n, int p);
  
  /** \brief Construct a block sparsity pattern from (row,col) vectors
      
   */
  Sparsity sp_rowcol(const std::vector<int>& row, const std::vector<int>& col, int nrow, int ncol);
  
  
  /** \brief Get the indices of all non-zero elements as they would appear in a Dense matrix  
      A : DenseMatrix  4 x 3
      B : SparseMatrix 4 x 3 , 5 structural non-zeros
      
      k = getNZDense(A)
      A[k] will contain the elements of A that are non-zero in B         
  */
  std::vector<int> getNZDense(const Sparsity& sp);
  
  /** \ brief Reshape the sparsity pattern keeping the relative location of the nonzeros
   */
  Sparsity reshape(const Sparsity& a, int nrow, int ncol);

  /** \ brief Vectorize the pattern */
  Sparsity vec(const Sparsity& a);
  
  Sparsity flatten(const Sparsity& a);
  
  /** \ brief Return the transpose of the sparsity pattern
   */
  Sparsity trans(const Sparsity& a);
  
  /**
   * \brief Return the upper part of the sparsity pattern
   * \param includeDiagonal specify wether the diagonal must be part of the result
   * \see sp_triu
   */
  Sparsity upperSparsity(const Sparsity& a, bool includeDiagonal = true);

  /**
   * \brief Return the lower part of the sparsity pattern
   * \param includeDiagonal specify wether the diagonal must be part of the result
   * \see sp_tril
   */
  Sparsity lowerSparsity(const Sparsity& a, bool includeDiagonal = true);
  
  /** \brief Return the non-zero entries that make up the upper part of the provided matrix */
  std::vector<int> upperNZ(const Sparsity& a);

  /** \brief Return the non-zero entries that make up the lower part of the provided matrix */
  std::vector<int> lowerNZ(const Sparsity& a);
  
  /**
     \brief Create a sparsity pattern given the nonzeros in sparse triplet form
  **/
  Sparsity sp_triplet(int nrow, int ncol, const std::vector<int>& row, const std::vector<int>& col, std::vector<int>& mapping, bool invert_mapping=false);
  
  /**
     \brief Create a sparsity pattern given the nonzeros in sparse triplet form (no nonzero mapping)
     rows_are_sorted==true means that the row entries already in increasing order for each col and without any duplicates
  **/
  Sparsity sp_triplet(int nrow, int ncol, const std::vector<int>& row, const std::vector<int>& col);
  
  
  /** \brief Get the sparsity resulting from a matrix multiplication
   */
  Sparsity mul(const  Sparsity& a, const  Sparsity &b);
  
  /** \brief Concatenate a list of sparsities vertically
  * Alternative terminology: vertical stack, vstack, vertical append, [a;b]
  */
  Sparsity vertcat(const std::vector<Sparsity > &v);

  /** \brief Concatenate a list of sparsities horizontally
  * Alternative terminology: horizontal stack, hstack, horizontal append, [a b]
  */
  Sparsity horzcat(const std::vector<Sparsity > &v);

  /** \brief   Construct a Sparsity with given blocks on the diagonal */
  Sparsity blkdiag(const std::vector< Sparsity > &v);

  #ifndef SWIG
  Sparsity horzcat(const Sparsity &x, const Sparsity &y);

  Sparsity vertcat(const Sparsity &x, const Sparsity &y);
  
  Sparsity blkdiag(const Sparsity &x, const Sparsity &y);
  #endif // SWIG
  
  /** \brief Split up a sparsity pattern horizontally */
  std::vector<Sparsity> horzsplit(const Sparsity& sp, const std::vector<int>& output_offset);

  /** \brief Split up a sparsity pattern vertically */
  std::vector<Sparsity> vertsplit(const Sparsity& sp, const std::vector<int>& output_offset);

  /** \brief Represent a sparsity pattern as an array of integers, the most compact way of representing a sparsity pattern
      The format:
      * The first two entries are the number of rows (nrow) and columns (ncol)
      * The next ncol+1 entries are the col offsets (colind). Note that the last element, colind[ncol], gives the number of nonzeros
      * The last colind[ncol] entries are the row indices
      **/
  /// @{
  
  /// Compress a sparsity pattern
  std::vector<int> sp_compress(const Sparsity& a);
  
  /// Decompress a sparsity pattern
  Sparsity sp_compress(const std::vector<int>& v);
  
#ifndef SWIG
  /// Decompress a sparsity pattern (array version)
  Sparsity sp_compress(const int* v);
#endif // SWIG  

  /// Obtain the structural rank of a sparsity-pattern
  int rank(const Sparsity& a);
  
  /// Check whether the sparsity-pattern inidcates structural singularity
  bool isSingular(const Sparsity& a);

  /// @}


}

#endif // SPARSITY_TOOLS_HPP
