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


#include "fixed_smith_dle_internal.hpp"
#include <cassert>
#include "../core/std_vector_tools.hpp"
#include "../core/matrix/matrix_tools.hpp"
#include "../core/mx/mx_tools.hpp"
#include "../core/sx/sx_tools.hpp"
#include "../core/function/mx_function.hpp"
#include "../core/function/sx_function.hpp"
#include <iomanip>

#include <numeric>

INPUTSCHEME(DLEInput)

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_DLESOLVER_FIXED_SMITH_EXPORT
  casadi_register_dlesolver_fixed_smith(DleInternal::Plugin* plugin) {
    plugin->creator = FixedSmithDleInternal::creator;
    plugin->name = "fixed_smith";
    plugin->doc = FixedSmithDleInternal::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_DLESOLVER_FIXED_SMITH_EXPORT casadi_load_dlesolver_fixed_smith() {
    DleInternal::registerPlugin(casadi_register_dlesolver_fixed_smith);
  }

  FixedSmithDleInternal::FixedSmithDleInternal(
      const DleStructure& st) :
      DleInternal(st) {

    // set default options
    setOption("name", "unnamed_fixed_smith_indef_dle_solver"); // name of the function

    addOption("iter", OT_INTEGER, 100,   "Number of Smith iterations");

    addOption("freq_doubling", OT_BOOLEAN, false,   "Use frequency doubling");

  }

  FixedSmithDleInternal::~FixedSmithDleInternal() {

  }

  void FixedSmithDleInternal::init() {
    iter_  = getOption("iter");
    freq_doubling_ = getOption("freq_doubling");

    DleInternal::init();

    casadi_assert_message(!pos_def_,
      "pos_def option set to True: Solver only handles the indefinite case.");

    MX As = MX::sym("A", A_);
    MX Vs = MX::sym("V", V_);

    MX V = (Vs+Vs.T())/2;

    MX P = V;
    MX A = As;

    for (int i=0;i<iter_;++i) {
      P = mul(A, mul(P, A.T())) + V;
      if (freq_doubling_) {
        V = mul(A, mul(V, A.T())) + V;
        A = mul(A, A);
      }
    }

    f_ = MXFunction(dleIn("a", As, "v", Vs), dleOut("p", P));
    f_.init();

    Wrapper::checkDimensions();

  }

  void FixedSmithDleInternal::evaluate() {
    Wrapper::evaluate();
  }

  Function FixedSmithDleInternal::getDerivative(int nfwd, int nadj) {
    return f_.derivative(nfwd, nadj);
  }

  void FixedSmithDleInternal::deepCopyMembers(
      std::map<SharedObjectNode*, SharedObject>& already_copied) {
    DleInternal::deepCopyMembers(already_copied);
  }

  FixedSmithDleInternal* FixedSmithDleInternal::clone() const {
    // Return a deep copy
    FixedSmithDleInternal* node = new FixedSmithDleInternal(st_);
    node->setOption(dictionary());
    return node;
  }


} // namespace casadi


