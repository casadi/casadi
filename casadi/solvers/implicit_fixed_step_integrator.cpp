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

using namespace std;
namespace casadi {

  ImplicitFixedStepIntegrator::
  ImplicitFixedStepIntegrator(const std::string& name, const XProblem& dae)
    : FixedStepIntegrator(name, dae) {
    addOption("implicit_solver",               OT_STRING,  GenericType(),
              "An implicit function solver");
    addOption("implicit_solver_options",       OT_DICT, GenericType(),
              "Options to be passed to the NLP Solver");
  }

  ImplicitFixedStepIntegrator::~ImplicitFixedStepIntegrator() {
  }

  void ImplicitFixedStepIntegrator::init() {
    // Call the base class init
    FixedStepIntegrator::init();

    // Get the NLP creator function
    std::string implicit_function_name = option("implicit_solver");

    // Options
    Dict implicit_solver_options;
    if (hasSetOption("implicit_solver_options")) {
      implicit_solver_options = option("implicit_solver_options");
    }
    implicit_solver_options["implicit_input"] = DAE_Z;
    implicit_solver_options["implicit_output"] = DAE_ALG;

    // Allocate a solver
    implicit_solver_ =
      F_.rootfinder(name_ + "_implicit_solver", implicit_function_name,
                    implicit_solver_options);

    // Allocate a root-finding solver for the backward problem
    if (nRZ_>0) {

      // Get the NLP creator function
      std::string backward_implicit_function_name = option("implicit_solver");

      // Options
      Dict backward_implicit_solver_options;
      if (hasSetOption("implicit_solver_options")) {
        backward_implicit_solver_options = option("implicit_solver_options");
      }
      backward_implicit_solver_options["implicit_input"] = RDAE_RZ;
      backward_implicit_solver_options["implicit_output"] = RDAE_ALG;

      // Allocate an NLP solver
      backward_implicit_solver_ =
        G_.rootfinder(name_+ "_backward_implicit_solver", backward_implicit_function_name,
                      backward_implicit_solver_options);
    }
  }

} // namespace casadi
