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

/**
 Create a dense rectangular sparsity pattern
**/
CRSSparsity sp_dense(int n, int m=1);

/**
 Create a sparse rectangular sparsity pattern
**/
CRSSparsity sp_sparse(int n, int m=1);

/**
 Create an upper triangular square sparsity pattern
**/
CRSSparsity sp_tril(int n);

/**
 Create diagonal square sparsity pattern
**/
CRSSparsity sp_diag(int n);

/** \brief Construct a sparsity pattern from (row,col) vectors

  row and column must be of same length.
  
  If you can guarantee that row is montone, pass the extra argument as true.
*/
CRSSparsity sp_NZ(std::vector<int> row, std::vector<int> col,int nrow, int ncol, bool monotone=false);

/** \brief Construct a block sparsity pattern from (row,col) vectors

*/
CRSSparsity sp_rowcol(std::vector<int> row, std::vector<int> col,int nrow, int ncol);


/** \brief Get the indices of all non-zero elements as they would appear in a Dense matrix  
     A : DenseMatrix  4 x 3
     B : SparseMatrix 4 x 3 , 5 structural non-zeros
     
     k = getNZDense(A)
     A[k] will contain the elements of A that are non-zero in B         
*/
std::vector<int> getNZDense(const CRSSparsity& sp);
    
    
CRSSparsity reshape(const CRSSparsity& a, int n, int m);

CRSSparsity vec(const CRSSparsity& a);

}

#endif // SPARSITY_TOOLS_HPP
