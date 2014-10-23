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


#include "dple_to_dle.hpp"
#include <cassert>
#include "../core/std_vector_tools.hpp"
#include "../core/matrix/matrix_tools.hpp"
#include "../core/mx/mx_tools.hpp"
#include "../core/sx/sx_tools.hpp"
#include "../core/function/mx_function.hpp"
#include "../core/function/sx_function.hpp"

#include <numeric>

INPUTSCHEME(DLEInput)

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_DLESOLVER_DPLE_EXPORT
  casadi_register_dlesolver_dple(DleInternal::Plugin* plugin) {
    plugin->creator = DpleToDle::creator;
    plugin->name = "dple";
    plugin->doc = DpleToDle::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_DLESOLVER_DPLE_EXPORT casadi_load_dlesolver_dple() {
    DleInternal::registerPlugin(casadi_register_dlesolver_dple);
  }

  DpleToDle::DpleToDle(
         const DleStructure& st, const std::vector<int> &Hs) : DleInternal(st, Hs) {

    // set default options
    setOption("name", "unnamed_dple_to_dle"); // name of the function

    addOption("dple_solver",            OT_STRING, GenericType(),
              "User-defined DPLE solver class.");
    addOption("dple_solver_options",    OT_DICTIONARY,   GenericType(),
              "Options to be passed to the DPLE solver.");

  }

  DpleToDle::~DpleToDle() {

  }

  void DpleToDle::init() {

    DleInternal::init();

    // Create an dplesolver instance
    std::string dplesolver_name = getOption("dple_solver");
    dplesolver_ = DpleSolver(dplesolver_name, dpleStruct("a",
      std::vector<Sparsity>(1,A_), "v",std::vector<Sparsity>(1,V_),"c",std::vector<Sparsity>(1,C_),"h",std::vector<Sparsity>(1,H_))
    );

    if (hasSetOption("dple_solver_options")) {
      dplesolver_.setOption(getOption("dple_solver_options"));
    }

    // Initialize the NLP solver
    dplesolver_.init();

  }



  void DpleToDle::evaluate() {
    for (int i=0;i<getNumInputs();++i) {
      std::copy(input(i).begin(), input(i).end(), dplesolver_.input(i).begin());
    }
    dplesolver_.evaluate();
    for (int i=0;i<getNumOutputs();++i) {
      std::copy(dplesolver_.output(i).begin(), dplesolver_.output(i).end(), output(i).begin());
    }
  }

  Function DpleToDle::getDerivative(int nfwd, int nadj) {
    return dplesolver_.derivative(nfwd, nadj);
  }


  void DpleToDle::deepCopyMembers(
      std::map<SharedObjectNode*, SharedObject>& already_copied) {
    DleInternal::deepCopyMembers(already_copied);
  }

  DpleToDle* DpleToDle::clone() const {
    // Return a deep copy
    DpleToDle* node = new DpleToDle(st_, Hs_);
    node->setOption(dictionary());
    return node;
  }


} // namespace casadi


