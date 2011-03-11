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

#ifndef CSPARSE_HPP
#define CSPARSE_HPP

#include "casadi/fx/linear_solver.hpp"

namespace CasADi{
  namespace Interfaces{

  
/** \brief  Forward declaration of internal class */
class CSparseInternal;

/** \brief  Public class */
class CSparse : public LinearSolver{
public:

  /// Default (empty) constructor
  CSparse();
  
  /// Create a linear solver given a sparsity pattern
  CSparse(const CRSSparsity& sp);
  
  /** \brief  Access internal functions and data members */
  CSparseInternal* operator->();
  
  /** \brief  Access internal functions and data members */
  const CSparseInternal* operator->() const;
  
  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
};

  } // namespace Interface
} // namespace CasADi

#endif //CSPARSE_HPP

