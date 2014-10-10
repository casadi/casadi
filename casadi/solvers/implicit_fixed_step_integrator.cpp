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


#include "implicit_fixed_step_integrator.hpp"
#include "casadi/core/std_vector_tools.hpp"
#include "casadi/core/matrix/sparsity_tools.hpp"
#include "casadi/core/matrix/matrix_tools.hpp"
#include "casadi/core/sx/sx_tools.hpp"
#include "casadi/core/function/sx_function.hpp"
#include "casadi/core/mx/mx_tools.hpp"

using namespace std;
namespace casadi {

  ImplicitFixedStepIntegrator::ImplicitFixedStepIntegrator(const Function& f,
                                                                           const Function& g)
      : FixedStepIntegrator(f, g) {
    addOption("implicit_solver",               OT_STRING,  GenericType(),
              "An implicit function solver");
    addOption("implicit_solver_options",       OT_DICTIONARY, GenericType(),
              "Options to be passed to the NLP Solver");
  }

  void ImplicitFixedStepIntegrator::deepCopyMembers(
      std::map<SharedObjectNode*, SharedObject>& already_copied) {
    FixedStepIntegrator::deepCopyMembers(already_copied);
    implicit_solver_ = deepcopy(implicit_solver_, already_copied);
    backward_implicit_solver_ = deepcopy(backward_implicit_solver_, already_copied);
  }

  ImplicitFixedStepIntegrator::~ImplicitFixedStepIntegrator() {
  }

  void ImplicitFixedStepIntegrator::init() {
    // Call the base class init
    FixedStepIntegrator::init();

    // Get the NLP creator function
    std::string implicit_function_name = getOption("implicit_solver");

    // Allocate an NLP solver
    implicit_solver_ = ImplicitFunction(implicit_function_name, F_, Function(), LinearSolver());
    implicit_solver_.setOption("name", string(getOption("name")) + "_implicit_solver");
    implicit_solver_.setOption("implicit_input", DAE_Z);
    implicit_solver_.setOption("implicit_output", DAE_ALG);

    // Pass options
    if (hasSetOption("implicit_solver_options")) {
      const Dictionary& implicit_solver_options = getOption("implicit_solver_options");
      implicit_solver_.setOption(implicit_solver_options);
    }

    // Initialize the solver
    implicit_solver_.init();

    // Allocate a root-finding solver for the backward problem
    if (nRZ_>0) {

      // Get the NLP creator function
      std::string backward_implicit_function_name = getOption("implicit_solver");

      // Allocate an NLP solver
      backward_implicit_solver_ = ImplicitFunction(backward_implicit_function_name,
                                                   G_, Function(), LinearSolver());
      backward_implicit_solver_.setOption("name",
                                          string(getOption("name")) + "_backward_implicit_solver");
      backward_implicit_solver_.setOption("implicit_input", RDAE_RZ);
      backward_implicit_solver_.setOption("implicit_output", RDAE_ALG);

      // Pass options
      if (hasSetOption("implicit_solver_options")) {
        const Dictionary& backward_implicit_solver_options = getOption("implicit_solver_options");
        backward_implicit_solver_.setOption(backward_implicit_solver_options);
      }

      // Initialize the solver
      backward_implicit_solver_.init();
    }
  }

} // namespace casadi
