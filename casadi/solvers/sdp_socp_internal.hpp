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

#ifndef CASADI_SDP_SOCP_INTERNAL_HPP
#define CASADI_SDP_SOCP_INTERNAL_HPP

#include "casadi/core/function/socp_solver_internal.hpp"
#include "casadi/core/function/sdp_solver.hpp"

#include <casadi/solvers/casadi_socpsolver_sdp_export.h>

/** \defgroup plugin_SocpSolver_sdp
   Solve SOCPs using an SdpSolver
*/


/** \pluginsection{SocpSolver,sdp} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{SocpSolver,sdp}

      @copydoc SocpSolver_doc
      @copydoc plugin_SocpSolver_sdp

      \author Joris Gillis
      \date 2013
  */
  class CASADI_SOCPSOLVER_SDP_EXPORT SDPSOCPInternal : public SocpSolverInternal {
    friend class SDPSocpSolver;
  public:

    /** \brief  Create a new Solver */
    explicit SDPSOCPInternal(const std::vector<Sparsity> &st);

    /** \brief  Destructor */
    virtual ~SDPSOCPInternal();

    /** \brief  Clone */
    virtual SDPSOCPInternal* clone() const;
  
    /** \brief  Create a new SOCP Solver */
    static SocpSolverInternal* creator(const SOCPStructure& st)
    { return new SDPSOCPInternal(st);}

    /** \brief  Initialize */
    virtual void init();

    virtual void evaluate();

    /// A documentation string
    static const std::string meta_doc;

  protected:
    SdpSolver sdpsolver_;

    /// Mapping from G, H, E, F to P, G
    Function mapping_;

  };

} // namespace casadi

/// \endcond
#endif // CASADI_SDP_SOCP_INTERNAL_HPP
