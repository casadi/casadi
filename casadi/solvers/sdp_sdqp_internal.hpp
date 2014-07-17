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

#ifndef CASADI_SDP_SDQP_INTERNAL_HPP
#define CASADI_SDP_SDQP_INTERNAL_HPP

#include "casadi/core/function/sdqp_solver_internal.hpp"
#include "casadi/core/function/sdp_solver.hpp"
#include "casadi/core/function/linear_solver.hpp"

#include <casadi/solvers/casadi_sdqpsolver_sdp_export.h>

/// \cond INTERNAL
namespace casadi {

  /** \brief SDP SDQP Solver for quadratic programming
   *
   *  Note: this implementation relies on Cholesky decomposition:
   *        <tt>Chol(H) = L ->  H = LL'</tt> with L lower triangular
   *   This requires Pi, H to be positive definite. Positive semi-definite is not sufficient.
   *    Notably, H==0  will not work.
   *
   *  A better implementation would rely on matrix square root,
   *  but we need singular value decomposition to implement that.
   *
   *
   @copydoc SdqpSolver_doc

   \author Joris Gillis
   \date 2013
  */
  class CASADI_SDQPSOLVER_SDP_EXPORT SDPSDQPInternal : public SdqpSolverInternal {
    friend class SDPSdqpSolver;
  public:

    /** \brief Constructor */
    explicit SDPSDQPInternal(const std::vector<Sparsity> &st);

    /** \brief Clone */
    virtual SDPSDQPInternal* clone() const { return new SDPSDQPInternal(*this);}

    /** \brief  Create a new SDQP Solver */
    static SdqpSolverInternal* creator(const SDQPStructure& st)
    { return new SDPSDQPInternal(st);}

    /// Deep copy data members
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief Destructor */
    virtual ~SDPSDQPInternal();

    /** \brief Initialize */
    virtual void init();

    /** \brief Solve the SDQP */
    virtual void evaluate();

    /// Underlying SDP solver
    SdpSolver sdpsolver_;

    /// Cholesky Decomposition
    LinearSolver cholesky_;

    /// Mapping
    Function mapping_;
  };

} // namespace casadi

/// \endcond
#endif // CASADI_SDP_SDQP_INTERNAL_HPP
