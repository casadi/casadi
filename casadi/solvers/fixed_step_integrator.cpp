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

using namespace std;
namespace casadi {

  FixedStepIvpsol::FixedStepIvpsol(const std::string& name, const XProblem& dae)
      : Ivpsol(name, dae) {

    addOption("number_of_finite_elements",     OT_INTEGER,  20, "Number of finite elements");
  }

  FixedStepIvpsol::~FixedStepIvpsol() {
  }

  void FixedStepIvpsol::init() {
    // Call the base class init
    Ivpsol::init();

    // Number of finite elements and time steps
    nk_ = option("number_of_finite_elements");
    casadi_assert(nk_>0);
    h_ = (grid_.back() - grid_.front())/nk_;

    // Setup discrete time dynamics
    setupFG();

    // Get discrete time dimensions
    Z_ = F_.input(DAE_Z);
    nZ_ = Z_.nnz();
    RZ_ = G_.isNull() ? DMatrix() : G_.input(RDAE_RZ);
    nRZ_ =  RZ_.nnz();

    // Allocate tape if backward states are present
    if (nrx_>0) {
      x_tape_.resize(nk_+1, vector<double>(nx_));
      Z_tape_.resize(nk_, vector<double>(nZ_));
    }
  }

  void FixedStepIvpsol::advance(int k) {
    // Get discrete time sought
    int k_out = std::ceil((grid_.at(k) - grid_.front())/h_);
    k_out = std::min(k_out, nk_); //  make sure that rounding errors does not result in k_out>nk_
    casadi_assert(k_out>=0);

    // Explicit discrete time dynamics
    Function& F = getExplicit();

    // Take time steps until end time has been reached
    while (k_<k_out) {
      // Take step
      F.input(DAE_T).set(t_);
      F.input(DAE_X).set(output(IVPSOL_XF));
      F.input(DAE_Z).set(Z_);
      F.input(DAE_P).set(input(IVPSOL_P));
      F.evaluate();
      F.output(DAE_ODE).get(output(IVPSOL_XF));
      F.output(DAE_ALG).get(Z_);
      copy(Z_->end()-nz_, Z_->end(), output(IVPSOL_ZF)->begin());
      transform(F.output(DAE_QUAD)->begin(),
                F.output(DAE_QUAD)->end(),
                output(IVPSOL_QF)->begin(),
                output(IVPSOL_QF)->begin(),
                std::plus<double>());

      // Tape
      if (nrx_>0) {
        output(IVPSOL_XF).getNZ(x_tape_.at(k_+1));
        Z_.getNZ(Z_tape_.at(k_));
      }

      // Advance time
      k_++;
      t_ = grid_.front() + k_*h_;
    }
  }

  void FixedStepIvpsol::retreat(int k) {
    // Get discrete time sought
    int k_out = std::floor((grid_.at(k) - grid_.front())/h_);
    k_out = std::max(k_out, 0); //  make sure that rounding errors does not result in k_out>nk_
    casadi_assert(k_out<=nk_);

    // Explicit discrete time dynamics
    Function& G = getExplicitB();

    // Take time steps until end time has been reached
    while (k_>k_out) {
      // Advance time
      k_--;
      t_ = grid_.front() + k_*h_;

      // Take step
      G.input(RDAE_T).set(t_);
      G.input(RDAE_X).setNZ(x_tape_.at(k_));
      G.input(RDAE_Z).setNZ(Z_tape_.at(k_));
      G.input(RDAE_P).set(input(IVPSOL_P));
      G.input(RDAE_RX).set(output(IVPSOL_RXF));
      G.input(RDAE_RZ).set(RZ_);
      G.input(RDAE_RP).set(input(IVPSOL_RP));
      G.evaluate();
      G.output(RDAE_ODE).get(output(IVPSOL_RXF));
      G.output(RDAE_ALG).get(RZ_);
      copy(RZ_->end()-nrz_, RZ_->end(), output(IVPSOL_RZF)->begin());
      transform(G.output(RDAE_QUAD)->begin(),
                G.output(RDAE_QUAD)->end(),
                output(IVPSOL_RQF)->begin(),
                output(IVPSOL_RQF)->begin(),
                std::plus<double>());
    }
  }

  void FixedStepIvpsol::reset(const double** arg, double** res, int* iw, double* w) {
    // Reset the base classes
    Ivpsol::reset(arg, res, iw, w);

    // Bring discrete time to the beginning
    k_ = 0;

    // Get consistent initial conditions
    calculateInitialConditions();

    // Add the first element in the tape
    if (nrx_>0) {
      output(IVPSOL_XF).getNZ(x_tape_.at(0));
    }
  }

  void FixedStepIvpsol::resetB() {
    // Reset the base classes
    Ivpsol::resetB();

    // Bring discrete time to the end
    k_ = nk_;

    // Get consistent initial conditions
    calculateInitialConditionsB();
  }

  void FixedStepIvpsol::calculateInitialConditions() {
    Z_.set(numeric_limits<double>::quiet_NaN());
  }

  void FixedStepIvpsol::calculateInitialConditionsB() {
    RZ_.set(numeric_limits<double>::quiet_NaN());
  }

} // namespace casadi
