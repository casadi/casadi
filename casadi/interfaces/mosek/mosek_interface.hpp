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


#ifndef CASADI_MOSEK_INTERFACE_HPP
#define CASADI_MOSEK_INTERFACE_HPP

#include "casadi/core/function/sdp_solver_internal.hpp"
#include <casadi/interfaces/mosek/casadi_sdpsolver_mosek_export.h>
#include <mosek5.h>

/** \defgroup plugin_SdpSolver_mosek
      Interface to the SDP solver MOSEK
      Warning: The solver MOSEK is not good at handling linear equalities.
      There are several options if you notice difficulties:
      * play around with the parameter "_penalty"
      * leave a gap manually
      * switch to another SDP Solver
*/

/** \pluginsection{SdpSolver,mosek} */

/// \cond INTERNAL
namespace casadi {
  /** \brief \pluginbrief{SdpSolver,mosek}


   \author Joris Gillis
   \date 2013
   *
   @copydoc SdpSolver_doc
   @copydoc plugin_SdpSolver_mosek
   * */
  class CASADI_SDPSOLVER_MOSEK_EXPORT MosekInterface : public SdpSolverInternal {
  public:

    /** \brief Constructor */
    explicit MosekInterface(const std::vector<Sparsity> &st);

    /** \brief Clone */
    virtual MosekInterface* clone() const;

    /** \brief  Create a new SDP Solver */
    static SdpSolverInternal* creator(const SDPStructure& st)
    { return new MosekInterface(st);}

    /** \brief Destructor */
    virtual ~MosekInterface();

    /** \brief Initialize */
    virtual void init();

    /** \brief Solve the SDP */
    virtual void evaluate();

  };

} // namespace casadi
/// \endcond

#endif // CASADI_MOSEK_INTERFACE_HPP
