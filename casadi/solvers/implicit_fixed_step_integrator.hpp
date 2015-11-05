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


#ifndef CASADI_IMPLICIT_FIXED_STEP_IVPSOL_HPP
#define CASADI_IMPLICIT_FIXED_STEP_IVPSOL_HPP

#include "fixed_step_integrator.hpp"
#include <casadi/solvers/casadi_ivpsols_export.h>

/// \cond INTERNAL
namespace casadi {

  class CASADI_IVPSOLS_EXPORT ImplicitFixedStepIvpsol
      : public FixedStepIvpsol {
  public:

    /// Constructor
    explicit ImplicitFixedStepIvpsol(const std::string& name, const XProblem& dae);

    /// Destructor
    virtual ~ImplicitFixedStepIvpsol();

    /// Initialize stage
    virtual void init();

    /// Get explicit dynamics
    virtual Function& getExplicit() { return implicit_solver_;}

    /// Get explicit dynamics (backward problem)
    virtual Function& getExplicitB() { return backward_implicit_solver_;}

    // Implicit function solver
    Function implicit_solver_, backward_implicit_solver_;

  };

} // namespace casadi
/// \endcond
#endif // CASADI_IMPLICIT_FIXED_STEP_IVPSOL_HPP
