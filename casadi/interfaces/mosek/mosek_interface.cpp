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


#include "mosek_interface.hpp"

#include "casadi/core/std_vector_tools.hpp"
#include "casadi/core/matrix/matrix_tools.hpp"
#include "casadi/core/mx/mx_tools.hpp"
#include "casadi/core/function/mx_function.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_SDPSOLVER_MOSEK_EXPORT
  casadi_register_sdpsolver_mosek(SdpSolverInternal::Plugin* plugin) {
    plugin->creator = MosekInterface::creator;
    plugin->name = "mosek";
    plugin->doc = MosekInterface::meta_doc.c_str();
    plugin->version = 22;
    return 0;
  }

  extern "C"
  void CASADI_SDPSOLVER_MOSEK_EXPORT casadi_load_sdpsolver_mosek() {
    SdpSolverInternal::registerPlugin(casadi_register_sdpsolver_mosek);
  }

  MosekInterface* MosekInterface::clone() const {
    // Return a deep copy
    MosekInterface* node = new MosekInterface(st_);
    if (!node->is_init_)
      node->init();
    return node;
  }

  MosekInterface::MosekInterface(const std::vector<Sparsity> &st) : SdpSolverInternal(st) {
 

  }

  MosekInterface::~MosekInterface() {

  }

  const char* MosekInterface::terminationReason(int flag) {

  }

  const char* MosekInterface::solutionType(int flag) {

  }

  void MosekInterface::init() {
    // Initialize the base classes
    SdpSolverInternal::init();
    log("MosekInterface::init", "Enter");


  }

  void MosekInterface::evaluate() {
    if (inputs_check_) checkInputs();
    if (print_problem_) printProblem();


  }

} // namespace casadi
