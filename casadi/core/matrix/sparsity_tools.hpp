/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


#ifndef CASADI_SPARSITY_TOOLS_HPP
#define CASADI_SPARSITY_TOOLS_HPP

#include "sparsity.hpp"

namespace casadi {

/**
\ingroup expression_tools
@{
*/

  /** \brief Reshape the sparsity pattern keeping the relative location of the nonzeros
   */
  CASADI_CORE_EXPORT Sparsity reshape(const Sparsity& a, int nrow, int ncol);

  /** \brief Transpose the pattern */
  inline Sparsity transpose(const Sparsity& a) { return a.transpose();}

  /** \brief Vectorize the pattern */
  CASADI_CORE_EXPORT Sparsity vec(const Sparsity& a);

  /** \brief Get the sparsity resulting from a matrix multiplication
   */
  CASADI_CORE_EXPORT Sparsity mul(const Sparsity& a, const Sparsity &b);

  /** \brief Get the sparsity resulting from a series of matrix multiplication
   */
  CASADI_CORE_EXPORT Sparsity mul(const std::vector<Sparsity>& s);

  /** \brief Concatenate a list of sparsities vertically
  * Alternative terminology: vertical stack, vstack, vertical append, [a;b]
  */
  CASADI_CORE_EXPORT Sparsity vertcat(const std::vector<Sparsity > &v);

  /** \brief Concatenate a list of sparsities horizontally
  * Alternative terminology: horizontal stack, hstack, horizontal append, [a b]
  */
  CASADI_CORE_EXPORT Sparsity horzcat(const std::vector<Sparsity > &v);

  /** \brief Construct a sparsity from a list of list of sparsities.
   */
  CASADI_CORE_EXPORT Sparsity blockcat(const std::vector< std::vector< Sparsity > > &v);

  /** \brief   Construct a Sparsity with given blocks on the diagonal */
  CASADI_CORE_EXPORT Sparsity blkdiag(const std::vector< Sparsity > &v);

  #ifndef SWIG
  CASADI_CORE_EXPORT Sparsity horzcat(const Sparsity &x, const Sparsity &y);

  CASADI_CORE_EXPORT Sparsity vertcat(const Sparsity &x, const Sparsity &y);

  CASADI_CORE_EXPORT Sparsity blkdiag(const Sparsity &x, const Sparsity &y);
  #endif // SWIG

  /** \brief Split up a sparsity pattern horizontally */
  CASADI_CORE_EXPORT
    std::vector<Sparsity> horzsplit(const Sparsity& sp, const std::vector<int>& output_offset);

  /** \brief Split up a sparsity pattern vertically */
  CASADI_CORE_EXPORT
    std::vector<Sparsity> vertsplit(const Sparsity& sp, const std::vector<int>& output_offset);

  /** \brief Split up a sparsity pattern diagonally */
  CASADI_CORE_EXPORT
    std::vector<Sparsity> diagsplit(const Sparsity& sp,
      const std::vector<int>& output_offset1,
      const std::vector<int>& output_offset2);

  /** \brief  split diagonally, retaining groups of square matrices
  * \param incr Size of each matrix
  *
  *  diagsplit(diagsplit(x, ...)) = x
  */
  CASADI_CORE_EXPORT std::vector<Sparsity> diagsplit(const Sparsity& x, int incr=1);

  /** \brief  split diagonally, retaining fixed-sized matrices
  * \param incr1 Row dimension of each matrix
  * \param incr2 Column dimension of each matrix
  *
  *  diagsplit(diagsplit(x, ...)) = x
  */
  CASADI_CORE_EXPORT std::vector<Sparsity> diagsplit(const Sparsity& x, int incr1, int incr2);

  /** \brief  split diagonally, retaining square matrices
  * \param output_offset List of all start locations for each group
  *      the last matrix will run to the end.
  *
  *   diagcat(diagsplit(x, ...)) = x
  */
  CASADI_CORE_EXPORT std::vector<Sparsity> diagsplit(const Sparsity& x,
                                                   const std::vector<int>& output_offset);

  /// Obtain the structural rank of a sparsity-pattern
  CASADI_CORE_EXPORT int rank(const Sparsity& a);

  /// Get upper triangular part
  inline Sparsity triu(const Sparsity& sp, bool includeDiagonal=true)
  { return sp.getTriu(includeDiagonal);}

  /// Get lower triangular part
  inline Sparsity tril(const Sparsity& sp, bool includeDiagonal=true)
  { return sp.getTril(includeDiagonal);}

  /*
  @}
  */
} // namespace casadi

#endif // CASADI_SPARSITY_TOOLS_HPP
