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
#include "integrator_internal.hpp"
#include "../std_vector_tools.hpp"
#include "sx_function.hpp"
#include "../sx/sx_tools.hpp"

INPUTSCHEME(IntegratorInput)

using namespace std;
namespace casadi {


  SimulatorInternal::SimulatorInternal(const Integrator& integrator,
                                       const Function& output_fcn,
                                       const vector<double>& grid) :
      integrator_(integrator), output_fcn_(output_fcn), grid_(grid) {
    setOption("name", "unnamed simulator");
    addOption("monitor",      OT_STRINGVECTOR, GenericType(),  "", "initial|step", true);

    input_.scheme = SCHEME_IntegratorInput;
  }

  SimulatorInternal::~SimulatorInternal() {
  }

  void SimulatorInternal::deepCopyMembers(
      std::map<SharedObjectNode*, SharedObject>& already_copied) {
    FunctionInternal::deepCopyMembers(already_copied);
    integrator_ = deepcopy(integrator_, already_copied);
    output_fcn_ = deepcopy(output_fcn_, already_copied);
  }

  void SimulatorInternal::init() {
    // Let the integration time start from the first point of the time grid.
    if (!grid_.empty()) integrator_.setOption("t0", grid_[0]);
    // Let the integration time stop at the last point of the time grid.
    if (!grid_.empty()) integrator_.setOption("tf", grid_[grid_.size()-1]);

    casadi_assert_message(isNonDecreasing(grid_), "The supplied time grid must be non-decreasing.");

    // Initialize the integrator
    integrator_.init();

    // Generate an output function if there is none (returns the whole state)
    if (output_fcn_.isNull()) {
      SX t = SX::sym("t");
      SX x = SX::sym("x", integrator_.input(INTEGRATOR_X0).sparsity());
      SX z = SX::sym("z", integrator_.input(INTEGRATOR_Z0).sparsity());
      SX p = SX::sym("p", integrator_.input(INTEGRATOR_P).sparsity());

      vector<SX> arg(DAE_NUM_IN);
      arg[DAE_T] = t;
      arg[DAE_X] = x;
      arg[DAE_Z] = z;
      arg[DAE_P] = p;

      vector<SX> out(INTEGRATOR_NUM_OUT);
      out[INTEGRATOR_XF] = x;
      out[INTEGRATOR_ZF] = z;

      // Create the output function
      output_fcn_ = SXFunction(arg, out);

      output_.scheme = SCHEME_IntegratorOutput;
    }

    // Initialize the output function
    output_fcn_.init();

    // Allocate inputs
    setNumInputs(INTEGRATOR_NUM_IN);
    for (int i=0; i<INTEGRATOR_NUM_IN; ++i) {
      input(i) = integrator_.input(i);
    }

    // Allocate outputs
    setNumOutputs(output_fcn_->getNumOutputs());
    for (int i=0; i<getNumOutputs(); ++i) {
      output(i) = Matrix<double>::zeros(output_fcn_.output(i).numel(), grid_.size());
      if (!output_fcn_.output(i).isEmpty()) {
        casadi_assert_message(output_fcn_.output(i).isVector(),
                              "SimulatorInternal::init: Output function output #" << i
                              << " has shape " << output_fcn_.output(i).dimString()
                              << ", while a column-matrix shape is expected.");
      }
    }

    casadi_assert_message(output_fcn_.input(DAE_T).numel() <=1,
                          "SimulatorInternal::init: output_fcn DAE_T argument must be "
                          "scalar or empty, but got " << output_fcn_.input(DAE_T).dimString());

    casadi_assert_message(
        output_fcn_.input(DAE_P).isEmpty() ||
        integrator_.input(INTEGRATOR_P).sparsity() == output_fcn_.input(DAE_P).sparsity(),
        "SimulatorInternal::init: output_fcn DAE_P argument must be empty or"
        << " have dimension " << integrator_.input(INTEGRATOR_P).dimString()
        << ", but got " << output_fcn_.input(DAE_P).dimString());

    casadi_assert_message(
        output_fcn_.input(DAE_X).isEmpty() ||
        integrator_.input(INTEGRATOR_X0).sparsity() == output_fcn_.input(DAE_X).sparsity(),
        "SimulatorInternal::init: output_fcn DAE_X argument must be empty or have dimension "
        << integrator_.input(INTEGRATOR_X0).dimString()
        << ", but got " << output_fcn_.input(DAE_X).dimString());

    // Call base class method
    FunctionInternal::init();

    // Output iterators
    output_its_.resize(getNumOutputs());
  }

  void SimulatorInternal::evaluate() {

    // Pass the parameters and initial state
    integrator_.setInput(input(INTEGRATOR_X0), INTEGRATOR_X0);
    integrator_.setInput(input(INTEGRATOR_Z0), INTEGRATOR_Z0);
    integrator_.setInput(input(INTEGRATOR_P), INTEGRATOR_P);

    if (monitored("initial")) {
      std::cout << "SimulatorInternal::evaluate: initial condition:" << std::endl;
      std::cout << " x0     = "  << input(INTEGRATOR_X0) << std::endl;
      std::cout << " z0     = "  << input(INTEGRATOR_Z0) << std::endl;
      std::cout << " p      = "   << input(INTEGRATOR_P) << std::endl;
    }

    // Reset the integrator_
    integrator_.reset();

    // Iterators to output data structures
    for (int i=0; i<output_its_.size(); ++i) output_its_[i] = output(i).begin();

    // Advance solution in time
    for (int k=0; k<grid_.size(); ++k) {

      if (monitored("step")) {
        std::cout << "SimulatorInternal::evaluate: integrating up to: " <<  grid_[k] << std::endl;
        std::cout << " x0       = "  << integrator_.input(INTEGRATOR_X0) << std::endl;
        std::cout << " z0       = "  << integrator_.input(INTEGRATOR_Z0) << std::endl;
        std::cout << " p        = "   << integrator_.input(INTEGRATOR_P) << std::endl;
      }

      // Integrate to the output time
      integrator_.integrate(grid_[k]);

      if (monitored("step")) {
        std::cout << " xf  = "  << integrator_.output(INTEGRATOR_XF) << std::endl;
        std::cout << " zf  = "  << integrator_.output(INTEGRATOR_ZF) << std::endl;
      }

      // Pass integrator output to the output function
      if (output_fcn_.input(DAE_T).size()!=0)
        output_fcn_.setInput(grid_[k], DAE_T);
      if (output_fcn_.input(DAE_X).size()!=0)
        output_fcn_.setInput(integrator_.output(INTEGRATOR_XF), DAE_X);
      if (output_fcn_.input(DAE_Z).size()!=0)
        output_fcn_.setInput(integrator_.output(INTEGRATOR_ZF), DAE_Z);
      if (output_fcn_.input(DAE_P).size()!=0)
        output_fcn_.setInput(input(INTEGRATOR_P), DAE_P);

      // Evaluate output function
      output_fcn_.evaluate();

      // Save the outputs of the function
      for (int i=0; i<getNumOutputs(); ++i) {
        const Matrix<double> &res = output_fcn_.output(i);
        copy(res.begin(), res.end(), output_its_.at(i));
        output_its_.at(i) += res.size();
      }
    }

    // Consistency check
    for (int i=0; i<output_its_.size(); ++i) {
      casadi_assert(output_its_[i] == output(i).end());
    }
  }

} // namespace casadi
