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

#ifndef SDP_SDQP_INTERNAL_HPP
#define SDP_SDQP_INTERNAL_HPP

#include "symbolic/function/sdqp_solver_internal.hpp"
#include "symbolic/function/sdp_solver.hpp"
#include "interfaces/csparse/csparse_cholesky.hpp"

/// \cond INTERNAL
namespace CasADi{

  /** \brief Internal class for SDPSDQPInternal
   * 
   @copydoc SDQPSolver_doc
   * */
  class SDPSDQPInternal : public SDQPSolverInternal {
    friend class SDPSDQPSolver;
  public:

    /** \brief Constructor */
    explicit SDPSDQPInternal(const std::vector<Sparsity> &st);

    /** \brief Clone */
    virtual SDPSDQPInternal* clone() const{ return new SDPSDQPInternal(*this);}
  
    /// Deep copy data members
    virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);

    /** \brief Destructor */
    virtual ~SDPSDQPInternal();

    /** \brief Initialize */
    virtual void init();
  
    /** \brief Solve the SDQP */
    virtual void evaluate();
  
    /// Underlying SDP solver
    SDPSolver sdpsolver_;
    
    /// Cholseky Decomposition
    CSparseCholesky cholesky_;
    
    /// Mapping
    Function mapping_;
  };

} // namespace CasADi
/// \endcond
#endif //SDP_SDQP_INTERNAL_HPP

