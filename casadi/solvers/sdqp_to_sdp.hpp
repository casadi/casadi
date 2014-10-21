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


#ifndef CASADI_SDQP_TO_SDP_HPP
#define CASADI_SDQP_TO_SDP_HPP

#include "casadi/core/function/sdqp_solver_internal.hpp"
#include "casadi/core/function/sdp_solver_internal.hpp"
#include "casadi/core/function/linear_solver.hpp"

#include <casadi/solvers/casadi_sdqpsolver_sdp_export.h>

/** \defgroup plugin_SdqpSolver_sdp

    Solve an SQDP using an SdpSolver
   *  Note: this implementation relies on Cholesky decomposition:
   *        <tt>Chol(H) = L ->  H = LL'</tt> with L lower triangular
   *   This requires Pi, H to be positive definite. Positive semi-definite is not sufficient.
   *    Notably, H==0  will not work.
   *
   *  A better implementation would rely on matrix square root,
   *  but we need singular value decomposition to implement that.
   *
*/

/** \pluginsection{SdqpSolver,sdp} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{SdqpSolver,sdp}
   *
   *
   @copydoc SdqpSolver_doc
   @copydoc plugin_SdqpSolver_sdp

   \author Joris Gillis
   \date 2013
  */
  class CASADI_SDQPSOLVER_SDP_EXPORT SdqpToSdp : public SdqpSolverInternal,
                                                 public Adaptor<SdqpToSdp, SdpSolverInternal> {
  public:
    /// Solve with
    SdpSolver solver_;

    /** \brief Constructor */
    explicit SdqpToSdp(const std::vector<Sparsity> &st);

    /** \brief Clone */
    virtual SdqpToSdp* clone() const { return new SdqpToSdp(*this);}

    /** \brief  Create a new SDQP Solver */
    static SdqpSolverInternal* creator(const SDQPStructure& st)
    { return new SdqpToSdp(st);}

    /// Deep copy data members
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief Destructor */
    virtual ~SdqpToSdp();

    /** \brief Initialize */
    virtual void init();

    /** \brief Solve the SDQP */
    virtual void evaluate();

    /// Cholesky Decomposition
    LinearSolver cholesky_;

    /// Mapping
    Function mapping_;

    /// A documentation string
    static const std::string meta_doc;

  };

} // namespace casadi

/// \endcond
#endif // CASADI_SDQP_TO_SDP_HPP
