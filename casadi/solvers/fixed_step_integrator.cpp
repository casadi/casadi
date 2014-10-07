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


#include "fixed_step_integrator.hpp"
#include "casadi/core/std_vector_tools.hpp"
#include "casadi/core/matrix/sparsity_tools.hpp"
#include "casadi/core/matrix/matrix_tools.hpp"
#include "casadi/core/sx/sx_tools.hpp"
#include "casadi/core/function/sx_function.hpp"
#include "casadi/core/mx/mx_tools.hpp"

using namespace std;
namespace casadi {

  FixedStepIntegrator::FixedStepIntegrator(const Function& f,
                                                           const Function& g)
      : IntegratorInternal(f, g) {
    addOption("number_of_finite_elements",     OT_INTEGER,  20, "Number of finite elements");
  }

  void FixedStepIntegrator::deepCopyMembers(
      std::map<SharedObjectNode*, SharedObject>& already_copied) {
    IntegratorInternal::deepCopyMembers(already_copied);
    F_ = deepcopy(F_, already_copied);
    G_ = deepcopy(G_, already_copied);
  }

  FixedStepIntegrator::~FixedStepIntegrator() {
  }

  void FixedStepIntegrator::init() {
    // Call the base class init
    IntegratorInternal::init();

    // Number of finite elements and time steps
    nk_ = getOption("number_of_finite_elements");
    casadi_assert(nk_>0);
    h_ = (tf_ - t0_)/nk_;

    // Setup discrete time dynamics
    setupFG();

    // Get discrete time dimensions
    Z_ = F_.input(DAE_Z);
    nZ_ = Z_.size();
    RZ_ = G_.isNull() ? DMatrix() : G_.input(RDAE_RZ);
    nRZ_ =  RZ_.size();

    // Allocate tape if backward states are present
    if (nrx_>0) {
      x_tape_.resize(nk_+1, vector<double>(nx_));
      Z_tape_.resize(nk_, vector<double>(nZ_));
    }
  }

  void FixedStepIntegrator::integrate(double t_out) {
    // Get discrete time sought
    int k_out = std::ceil((t_out-t0_)/h_);
    k_out = std::min(k_out, nk_); //  make sure that rounding errors does not result in k_out>nk_
    casadi_assert(k_out>=0);

    // Explicit discrete time dynamics
    Function& F = getExplicit();

    // Take time steps until end time has been reached
    while (k_<k_out) {
      // Take step
      F.input(DAE_T).set(t_);
      F.input(DAE_X).set(output(INTEGRATOR_XF));
      F.input(DAE_Z).set(Z_);
      F.input(DAE_P).set(input(INTEGRATOR_P));
      F.evaluate();
      F.output(DAE_ODE).get(output(INTEGRATOR_XF));
      F.output(DAE_ALG).get(Z_);
      transform(F.output(DAE_QUAD).begin(),
                F.output(DAE_QUAD).end(),
                output(INTEGRATOR_QF).begin(),
                output(INTEGRATOR_QF).begin(),
                std::plus<double>());

      // Tape
      if (nrx_>0) {
        output(INTEGRATOR_XF).get(x_tape_.at(k_+1));
        Z_.get(Z_tape_.at(k_));
      }

      // Advance time
      k_++;
      t_ = t0_ + k_*h_;
    }
  }

  void FixedStepIntegrator::integrateB(double t_out) {
    // Get discrete time sought
    int k_out = std::floor((t_out-t0_)/h_);
    k_out = std::max(k_out, 0); //  make sure that rounding errors does not result in k_out>nk_
    casadi_assert(k_out<=nk_);

    // Explicit discrete time dynamics
    Function& G = getExplicitB();

    // Take time steps until end time has been reached
    while (k_>k_out) {
      // Advance time
      k_--;
      t_ = t0_ + k_*h_;

      // Take step
      G.input(RDAE_T).set(t_);
      G.input(RDAE_X).set(x_tape_.at(k_));
      G.input(RDAE_Z).set(Z_tape_.at(k_));
      G.input(RDAE_P).set(input(INTEGRATOR_P));
      G.input(RDAE_RX).set(output(INTEGRATOR_RXF));
      G.input(RDAE_RZ).set(RZ_);
      G.input(RDAE_RP).set(input(INTEGRATOR_RP));
      G.evaluate();
      G.output(RDAE_ODE).get(output(INTEGRATOR_RXF));
      G.output(RDAE_ALG).get(RZ_);
      transform(G.output(RDAE_QUAD).begin(),
                G.output(RDAE_QUAD).end(),
                output(INTEGRATOR_RQF).begin(),
                output(INTEGRATOR_RQF).begin(),
                std::plus<double>());
    }
  }

  void FixedStepIntegrator::reset() {
    // Reset the base classes
    IntegratorInternal::reset();

    // Bring discrete time to the beginning
    k_ = 0;

    // Get consistent initial conditions
    calculateInitialConditions();

    // Add the first element in the tape
    if (nrx_>0) {
      output(INTEGRATOR_XF).get(x_tape_.at(0));
    }
  }

  void FixedStepIntegrator::resetB() {
    // Reset the base classes
    IntegratorInternal::resetB();

    // Bring discrete time to the end
    k_ = nk_;

    // Get consistent initial conditions
    calculateInitialConditionsB();
  }

  void FixedStepIntegrator::calculateInitialConditions() {
    Z_.set(numeric_limits<double>::quiet_NaN());
  }

  void FixedStepIntegrator::calculateInitialConditionsB() {
    RZ_.set(numeric_limits<double>::quiet_NaN());
  }

} // namespace casadi
