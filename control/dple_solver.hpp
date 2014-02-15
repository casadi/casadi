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

#ifndef DPLE_SOLVER_HPP
#define DPLE_SOLVER_HPP

#include "../symbolic/fx/fx.hpp"

/** \defgroup DPLE_doc Discrete periodic Lyapunov Equation solver


    Given matrices A_k and symmetric V_k,  k = 0..K-1
    
    A_k in R^(n x n)
    V_k in R^n
    
    provides all of P_k that satisfy:

    P_0 = A_(K-1)*P_(K-1)*A_(K-1)' + V_k
    P_k+1 = A_k*P_k*A_k' + V_k  for k = 1..K-1
    
    

*/
namespace CasADi{

  /// Input arguments of a dple solver [dpleIn]
  enum DPLEInput{
    /// A matrices (vertcat when const_dim, blkdiag otherwise) [a]
    DPLE_A,
    /// V matrices (vertcat when const_dim, blkdiag otherwise) [v]
    DPLE_V,
    DPLE_NUM_IN
  };

  /// Output arguments of a dple solver [dpleOut]
  enum DPLEOutput{
    /// Lyapunov matrix (vertcat when const_dim, blkdiag otherwise) (cholesky of P if pos_def) [p]
    DPLE_P,
    /// Number of arguments.
    DPLE_NUM_OUT
  };

  /// Forward declaration of internal class
  class DpleInternal;

  /**  @copydoc DPLE_doc
       \author Joris gillis
      \date 2014
      
  */
  class DpleSolver : public FX {
  public:
    /// Default constructor
    DpleSolver();
    
    /// Clone
    DpleSolver clone() const;
  
    /// Print solver statistics
    void printStats(std::ostream &stream=std::cout) const;
  
    /// Access functions of the node
    DpleInternal* operator->();

    /// Access functions of the node
    const DpleInternal* operator->() const;
 
    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;
  
 
  };

} // namespace CasADi

#endif // DPLE_SOLVER_HPP
