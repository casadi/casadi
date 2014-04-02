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

/**
\ingroup expression_tools
@{ 
*/

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
  inline Sparsity sp_diag(int n){ return Sparsity::diag(n);}
  inline Sparsity sp_band(int n, int p){ return Sparsity::band(n,p);}
  inline Sparsity sp_banded(int n, int p){ return Sparsity::banded(n,p);}
  inline Sparsity sp_rowcol(const std::vector<int>& row, const std::vector<int>& col, int nrow, int ncol){ return Sparsity::rowcol(row,col,nrow,ncol);}
  inline Sparsity sp_triplet(int nrow, int ncol, const std::vector<int>& row, const std::vector<int>& col, std::vector<int>& mapping, bool invert_mapping=false){ return Sparsity::triplet(nrow,ncol,row,col,mapping,invert_mapping);}
  inline Sparsity sp_triplet(int nrow, int ncol, const std::vector<int>& row, const std::vector<int>& col){ return Sparsity::triplet(nrow,ncol,row,col);}
  inline Sparsity trans(const Sparsity& a){ return a.transpose();}
  inline Sparsity upperSparsity(const Sparsity& a, bool includeDiagonal = true){ return a.getTriu(includeDiagonal); }
  inline Sparsity lowerSparsity(const Sparsity& a, bool includeDiagonal = true){ return a.getTril(includeDiagonal); }
  inline std::vector<int> upperNZ(const Sparsity& a) { return a.getUpperNZ(); }
  inline std::vector<int> lowerNZ(const Sparsity& a) { return a.getLowerNZ(); }
  inline bool isSingular(const Sparsity& a){ return a.isSingular();}
  inline std::vector<int> sp_compress(const Sparsity& a){ return a.compress();}
  inline Sparsity sp_compress(const std::vector<int>& v){ return Sparsity::compressed(v);}
#ifndef SWIG
  inline Sparsity sp_compress(const int* v){ return Sparsity::compressed(v);}
#endif // SWIG
  inline std::vector<int> getNZDense(const Sparsity& sp){ return sp.getElements();}
  //@}
#endif  
    
  /** \brief Reshape the sparsity pattern keeping the relative location of the nonzeros
   */
  Sparsity reshape(const Sparsity& a, int nrow, int ncol);

  /** \brief Transpose the pattern */
  inline Sparsity transpose(const Sparsity& a){ return a.transpose();}

  /** \brief Vectorize the pattern */
  Sparsity vec(const Sparsity& a);
  
  /** \brief Get the sparsity resulting from a matrix multiplication
   */
  Sparsity mul(const Sparsity& a, const Sparsity &b);
  
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

  /// Obtain the structural rank of a sparsity-pattern
  int rank(const Sparsity& a);
  
  /// Get upper triangular part
  inline Sparsity triu(const Sparsity& sp, bool includeDiagonal=true){ return sp.getTriu(includeDiagonal);}

  /// Get lower triangular part
  inline Sparsity tril(const Sparsity& sp, bool includeDiagonal=true){ return sp.getTril(includeDiagonal);}
  
  /*
  @}
  */
}

#endif // SPARSITY_TOOLS_HPP
