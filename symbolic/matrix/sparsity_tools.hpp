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


#ifndef WITHOUT_PRE_1_9_X
    /** \brief [DEPRECATED] Replaced with static methods in the Sparsity class
     */
    //@{
  inline Sparsity sp_dense(int nrow, int ncol=1){ return Sparsity::dense(nrow,ncol);}
  inline Sparsity sp_dense(const std::pair<int,int> &rc ){ return Sparsity::dense(rc);}
  inline Sparsity sp_sparse(int nrow, int ncol=1){ return Sparsity::sparse(nrow,ncol);}
  inline Sparsity sp_sparse(const std::pair<int,int> &rc){ return Sparsity::sparse(rc);}
  inline Sparsity sp_unit(int n, int el){ return Sparsity::unit(n,el);}
  inline Sparsity sp_triu(int n){ return Sparsity::triu(n);}
  inline Sparsity sp_tril(int n){ return Sparsity::tril(n);}
  inline Sparsity sp_diag(int n){ return Sparsity::diagonal(n);}
  inline Sparsity sp_band(int n, int p){ return Sparsity::band(n,p);}
  inline Sparsity sp_banded(int n, int p){ return Sparsity::banded(n,p);}
  inline Sparsity sp_rowcol(const std::vector<int>& row, const std::vector<int>& col, int nrow, int ncol){ return Sparsity::rowcol(row,col,nrow,ncol);}
  inline Sparsity sp_triplet(int nrow, int ncol, const std::vector<int>& row, const std::vector<int>& col, std::vector<int>& mapping, bool invert_mapping=false){ return Sparsity::triplet(nrow,ncol,row,col,mapping,invert_mapping);}
  inline Sparsity sp_triplet(int nrow, int ncol, const std::vector<int>& row, const std::vector<int>& col){ return Sparsity::triplet(nrow,ncol,row,col);}
  inline Sparsity trans(const Sparsity& a){ return a.transpose();}
  inline Sparsity upperSparsity(const Sparsity& a, bool includeDiagonal = true){ return a.upper(includeDiagonal); }
  inline Sparsity lowerSparsity(const Sparsity& a, bool includeDiagonal = true){ return a.lower(includeDiagonal); }
  inline std::vector<int> upperNZ(const Sparsity& a) { return a.upperNZ(); }
  inline std::vector<int> lowerNZ(const Sparsity& a) { return a.lowerNZ(); }
  //@}
#endif  
  
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
