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

#ifndef SDP_SOCP_INTERNAL_HPP
#define SDP_SOCP_INTERNAL_HPP

#include "symbolic/fx/socp_solver_internal.hpp"
#include "symbolic/fx/sdp_solver.hpp"

namespace CasADi{

  /** \brief Internal class for SDPSOCPInternal
   * 
      @copydoc LPSolver_doc
   * */
class SDPSOCPInternal : public SOCPSolverInternal {
  friend class SDPSOCPSolver;
public:

  /** \brief  Clone */
  virtual SDPSOCPInternal* clone() const;
  
  /** \brief  Create a new Solver */
  explicit SDPSOCPInternal(const std::vector<CRSSparsity> &st);

  /** \brief  Destructor */
  virtual ~SDPSOCPInternal();

  /** \brief  Initialize */
  virtual void init();
  
  virtual void evaluate();
  
  protected:
    SDPSolver sdpsolver_;
    
    /// Mapping from G, H, E, F to P, G
    FX mapping_;

};

} // namespace CasADi

#endif //SDP_SOCP_INTERNAL_HPP

