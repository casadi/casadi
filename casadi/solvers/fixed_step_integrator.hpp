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


#ifndef CASADI_FIXED_STEP_INTEGRATOR_HPP
#define CASADI_FIXED_STEP_INTEGRATOR_HPP

#include "casadi/core/function/integrator_internal.hpp"
#include "casadi/core/function/mx_function.hpp"
#include <casadi/solvers/casadi_integrators_export.h>

/// \cond INTERNAL
namespace casadi {

  class CASADI_INTEGRATORS_EXPORT FixedStepIntegrator
    : public IntegratorInternal {
  public:

    /// Constructor
    explicit FixedStepIntegrator(const Function& f, const Function& g);

    /// Deep copy data members
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /// Clone
    virtual FixedStepIntegrator* clone() const = 0;

    /// Create a new integrator
    virtual FixedStepIntegrator* create(const Function& f, const Function& g) const = 0;

    /// Destructor
    virtual ~FixedStepIntegrator();

    /// Initialize stage
    virtual void init();

    /// Setup F and G
    virtual void setupFG() = 0;

    ///  Integrate until a specified time point
    virtual void integrate(double t_out);

    /// Integrate backward in time until a specified time point
    virtual void integrateB(double t_out);

    /// Reset the forward problem and bring the time back to t0
    virtual void reset();

    /// Reset the backward problem and take time to tf
    virtual void resetB();

    /// Get initial guess for the algebraic variable
    virtual void calculateInitialConditions();

    /// Get initial guess for the algebraic variable (backward problem)
    virtual void calculateInitialConditionsB();

    /// Get explicit dynamics
    virtual Function& getExplicit() { return F_;}

    /// Get explicit dynamics (backward problem)
    virtual Function& getExplicitB() { return G_;}

    // Discrete time dynamics
    Function F_, G_;

    // Number of finite elements
    int nk_;

    // Discrete time
    int k_;

    // Time step size
    double h_;

    /// Number of algebraic variables for the discrete time integration
    int nZ_, nRZ_;

    /// Algebraic variables for the discrete time integration
    DMatrix Z_, RZ_;

    // Tape
    std::vector<std::vector<double> > x_tape_, Z_tape_;

  };

} // namespace casadi
/// \endcond
#endif // CASADI_FIXED_STEP_INTEGRATOR_HPP
