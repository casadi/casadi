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

#ifndef OCP_SOLVER_INTERNAL_HPP
#define OCP_SOLVER_INTERNAL_HPP

#include <vector>
#include "ocp_solver.hpp"
#include "function_internal.hpp"

/// \cond INTERNAL

namespace casadi{

/** \brief  Internal node class for OCPSolver
  \author Joel Andersson
  \date 2010
*/
class CASADI_SYMBOLIC_EXPORT OCPSolverInternal : public FunctionInternal{
  friend class OCPSolver;
  public:

    /** \brief Constructor
    *
    *
    * \param ffcn Continuous time dynamics
    * \param mfcn Mayer term
    * \param cfcn Path constraints
    * \param rfcn Initial value constraints
    *
    * The signatures (number and order of inputs/outputs) of these functions are not restricted at this stage.
    *
    * Only ffcn has a general requirement for input interface: { DAE_T, DAE_Y, DAE_P, DAE_YDOT } from casadi::DAEInput
    *
    * For example:
    *
    * When using the ACADO interface, all functions should have the same input interface: casadi::ACADO_FCN_Input \n
    * When using MultipleShooting, mfcn_ is a single input -> single output mapping
    *
    */
    explicit OCPSolverInternal(const Function& ffcn, const Function& mfcn, const Function& cfcn, const Function& rfcn);

    /// Destructor
    virtual ~OCPSolverInternal();

    /// Initialize
    virtual void init();

    /// Discrete time dynamics
    Function ffcn_;

    /// Mayer term
    Function mfcn_;

    /// Path constraints
    Function cfcn_;

    /// Initial value constraints
    Function rfcn_;

    /// Number of grid points
    int nk_;

    /// Number of states
    int nx_;

    /// Number of parameters
    int np_;

    /// Number of controls
    int nu_;

    /// Number of point constraints
    int nh_;

    /// Number of point coupling constraints
    int ng_;

    /// Final time
    double tf_;

};



} // namespace casadi

/// \endcond
#endif // OCP_SOLVER_INTERNAL_HPP

