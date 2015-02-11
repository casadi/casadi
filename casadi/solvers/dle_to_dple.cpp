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


#include "dle_to_dple.hpp"
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
    plugin->creator = DleToDple::creator;
    plugin->name = "dple";
    plugin->doc = DleToDple::meta_doc.c_str();
    plugin->version = 22;
    return 0;
  }

  extern "C"
  void CASADI_DLESOLVER_DPLE_EXPORT casadi_load_dlesolver_dple() {
    DleInternal::registerPlugin(casadi_register_dlesolver_dple);
  }

  DleToDple::DleToDple(
         const DleStructure& st) : DleInternal(st) {

    // set default options
    setOption("name", "unnamed_dle_to_dple"); // name of the function

    Adaptor::addOptions();
  }

  DleToDple::~DleToDple() {

  }

  void DleToDple::init() {
    // Initialize the base classes
    DleInternal::init();

    // Create a DpleSolver instance
    solver_ = DpleSolver(getOption(solvername()),
                         dpleStruct("a", std::vector<Sparsity>(1, A_),
                                    "v", std::vector<Sparsity>(1, V_)));
    solver_.setOption(getOption(optionsname()));
    solver_.init();
  }

  void DleToDple::evaluate() {
    for (int i=0;i<getNumInputs();++i) {
      std::copy(input(i).begin(), input(i).end(), solver_.input(i).begin());
    }
    solver_.evaluate();
    for (int i=0;i<getNumOutputs();++i) {
      std::copy(solver_.output(i).begin(), solver_.output(i).end(), output(i).begin());
    }
  }

  Function DleToDple::getDerivativeFwd(int nfwd) {
    return solver_.derivativeFwd(nfwd);
  }

  Function DleToDple::getDerivativeAdj(int nadj) {
    return solver_.derivativeAdj(nadj);
  }

  void DleToDple::deepCopyMembers(
      std::map<SharedObjectNode*, SharedObject>& already_copied) {
    DleInternal::deepCopyMembers(already_copied);
  }

  DleToDple* DleToDple::clone() const {
    // Return a deep copy
    DleToDple* node = new DleToDple(st_);
    node->setOption(dictionary());
    return node;
  }


} // namespace casadi


