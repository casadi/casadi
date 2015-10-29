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


#include "simulator_internal.hpp"
#include "integrator.hpp"
#include "../std_vector_tools.hpp"

using namespace std;
namespace casadi {


  SimulatorInternal::SimulatorInternal(const std::string& name, const Function& integrator,
                                       const DMatrix& grid)
    : FunctionInternal(name), integrator_(integrator), grid_(grid.data()) {

    casadi_assert_message(grid.iscolumn(),
                          "Simulator::Simulator: grid must be a column vector, but got "
                          << grid.dim());
    casadi_assert_message(grid.isdense(),
                          "Simulator::Simulator: grid must be dense, but got "
                          << grid.dim());
    addOption("monitor",      OT_STRINGVECTOR, GenericType(),  "", "initial|step", true);
  }

  SimulatorInternal::~SimulatorInternal() {
  }

  void SimulatorInternal::init() {
    if (!grid_.empty()) {
      casadi_assert_message(isNonDecreasing(grid_),
                            "The supplied time grid must be non-decreasing.");

      // Create new integrator object
      auto dae = make_pair(dynamic_cast<Integrator*>(integrator_.get())->f_,
                           dynamic_cast<Integrator*>(integrator_.get())->g_);
      string solver = dynamic_cast<Integrator*>(integrator_.get())->plugin_name();
      Function I = Function::integrator(integrator_.name(), solver,
                                        dae, integrator_.dictionary());

      // Let the integration time start from the first point of the time grid.
      I.setOption("t0", grid_[0]);
      // Let the integration time stop at the last point of the time grid.
      I.setOption("tf", grid_[grid_.size()-1]);
      I.init();

      // Overwrite
      integrator_ = I;
    }

    ischeme_ = Function::integrator_in();
    oscheme_ = Function::integrator_out();

    // Allocate inputs
    ibuf_.resize(get_n_in());
    for (int i=0; i<ibuf_.size(); ++i) {
      input(i) = DMatrix(get_sparsity_in(i));
    }

    // Allocate outputs
    obuf_.resize(get_n_out());
    for (int i=0; i<obuf_.size(); ++i) {
      output(i) = DMatrix(get_sparsity_out(i));
    }

    // Call base class method
    FunctionInternal::init();
  }

  void SimulatorInternal::evaluate() {

    // Pass the parameters and initial state
    integrator_.setInput(input(INTEGRATOR_X0), INTEGRATOR_X0);
    integrator_.setInput(input(INTEGRATOR_Z0), INTEGRATOR_Z0);
    integrator_.setInput(input(INTEGRATOR_P), INTEGRATOR_P);

    if (monitored("initial")) {
      userOut() << "SimulatorInternal::evaluate: initial condition:" << std::endl;
      userOut() << " x0     = "  << input(INTEGRATOR_X0) << std::endl;
      userOut() << " z0     = "  << input(INTEGRATOR_Z0) << std::endl;
      userOut() << " p      = "   << input(INTEGRATOR_P) << std::endl;
    }

    // Reset the integrator_
    dynamic_cast<Integrator*>(integrator_.get())->reset();

    // Advance solution in time
    for (int k=0; k<grid_.size(); ++k) {

      if (monitored("step")) {
        userOut() << "SimulatorInternal::evaluate: integrating up to: " <<  grid_[k] << std::endl;
        userOut() << " x0       = "  << integrator_.input(INTEGRATOR_X0) << std::endl;
        userOut() << " z0       = "  << integrator_.input(INTEGRATOR_Z0) << std::endl;
        userOut() << " p        = "   << integrator_.input(INTEGRATOR_P) << std::endl;
      }

      // Integrate to the output time
      dynamic_cast<Integrator*>(integrator_.get())->integrate(grid_[k]);

      if (monitored("step")) {
        userOut() << " xf  = "  << integrator_.output(INTEGRATOR_XF) << std::endl;
        userOut() << " zf  = "  << integrator_.output(INTEGRATOR_ZF) << std::endl;
      }

      // Save the outputs of the function
      for (int i=0; i<n_out(); ++i) {
        const Matrix<double> &res = integrator_.output(i);
        copy(res.begin(), res.end(), output(i).begin() + k*res.nnz());
      }
    }
  }

} // namespace casadi
