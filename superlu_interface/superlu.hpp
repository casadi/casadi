/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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

#ifndef SUPERLU_HPP
#define SUPERLU_HPP

#include "casadi/fx/fx.hpp"

namespace CasADi{
  
/** \brief  Forward declaration of internal class */
class SuperLUInternal;

/// Type of solve call
enum Factorization{DOFACT, SAMEPATTERN, SAMEPATTERN_SAMEROWPERM, FACTORED};

/** \brief  Public class */
class SuperLU : public FX{
public:

  /// Default (empty) constructor
  SuperLU();
  
  /// Create a linear solver given a sparsity pattern
  SuperLU(int nrow, int ncol, const std::vector<int>& rowind, const std::vector<int>& col, int nrhs=1);
  
  /// Solve
  void solve(Factorization fact);
  
  /** \brief  Access functions of the node */
  SuperLUInternal* operator->();
  const SuperLUInternal* operator->() const;
};

} // namespace CasADi

#endif //SUPERLU_HPP

