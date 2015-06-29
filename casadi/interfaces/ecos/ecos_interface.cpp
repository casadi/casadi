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


#include "ecos_socp_interface.hpp"

#include "casadi/core/std_vector_tools.hpp"
#include "casadi/core/matrix/matrix_tools.hpp"
#include "casadi/core/mx/mx_tools.hpp"
#include "casadi/core/function/mx_function.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_SOCPSOLVER_ECOS_EXPORT
  casadi_register_socpsolver_ecos(SocpSolverInternal::Plugin* plugin) {
    plugin->creator = EcosInterface::creator;
    plugin->name = "ecos";
    plugin->doc = EcosInterface::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_SOCPSOLVER_ECOS_EXPORT casadi_load_socpsolver_ecos() {
    SocpSolverInternal::registerPlugin(casadi_register_socpsolver_ecos);
  }

  EcosInterface* EcosInterface::clone() const {
    // Return a deep copy
    EcosInterface* node = new EcosInterface(st_);
    if (!node->is_init_)
      node->init();
    return node;
  }

  EcosInterface::EcosInterface(const std::vector<Sparsity> &st) : SocpSolverInternal(st) {
    // Define options
  }

  EcosInterface::~EcosInterface() {
    if (is_init_) {
      // Destructor
    }
  }

  void EcosInterface::init() {
    // Initialize the base classes
    SocpSolverInternal::init();
    log("EcosInterface::init", "Enter");

   
  }

  void EcosInterface::evaluate() {
    if (inputs_check_) checkInputs();
    if (print_problem_) printProblem();

    // ECOS expects the SOCP in conic form, obtain this form by formulating dual SOCP
    convertToDualSocp();

  }


} // namespace casadi
