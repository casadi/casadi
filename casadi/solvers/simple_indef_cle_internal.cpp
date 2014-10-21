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


#include "simple_indef_cle_internal.hpp"
#include <cassert>
#include "../core/std_vector_tools.hpp"
#include "../core/matrix/matrix_tools.hpp"
#include "../core/mx/mx_tools.hpp"
#include "../core/sx/sx_tools.hpp"
#include "../core/function/mx_function.hpp"
#include "../core/function/sx_function.hpp"

#include <numeric>

INPUTSCHEME(CLEInput)
OUTPUTSCHEME(CLEOutput)

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_CLESOLVER_SIMPLE_EXPORT
  casadi_register_clesolver_simple(CleInternal::Plugin* plugin) {
    plugin->creator = SimpleIndefCleInternal::creator;
    plugin->name = "simple";
    plugin->doc = SimpleIndefCleInternal::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_CLESOLVER_SIMPLE_EXPORT casadi_load_clesolver_simple() {
    CleInternal::registerPlugin(casadi_register_clesolver_simple);
  }

  SimpleIndefCleInternal::SimpleIndefCleInternal(
      const CleStructure& st) : CleInternal(st) {

    // set default options
    setOption("name", "unnamed_simple_indef_cle_solver"); // name of the function

    addOption("linear_solver",            OT_STRING, GenericType(),
              "User-defined linear solver class. Needed for sensitivities.");
    addOption("linear_solver_options",    OT_DICTIONARY,   GenericType(),
              "Options to be passed to the linear solver.");

  }

  SimpleIndefCleInternal::~SimpleIndefCleInternal() {

  }

  void SimpleIndefCleInternal::init() {

    CleInternal::init();

    casadi_assert_message(!pos_def_,
      "pos_def option set to True: Solver only handles the indefinite case.");

    n_ = A_.size1();

    MX As = MX::sym("A", A_);
    MX Vs = MX::sym("V", V_);

    MX Vss = (Vs+Vs.T())/2;

    DMatrix e = DMatrix::eye(n_);

    MX A_total = -kron(As, e)-kron(e, As);
    MX Pf = solve(A_total, vec(Vss), getOption("linear_solver"));

    std::vector<MX> v_in;
    v_in.push_back(As);
    v_in.push_back(Vs);
    f_ = MXFunction(v_in, reshape(Pf, n_, n_));
    f_.setInputScheme(SCHEME_CLEInput);
    f_.setOutputScheme(SCHEME_CLEOutput);
    f_.init();
  }



  void SimpleIndefCleInternal::evaluate() {
    for (int i=0;i<getNumInputs();++i) {
      std::copy(input(i).begin(), input(i).end(), f_.input(i).begin());
    }
    f_.evaluate();
    for (int i=0;i<getNumOutputs();++i) {
      std::copy(f_.output(i).begin(), f_.output(i).end(), output(i).begin());
    }
  }

  Function SimpleIndefCleInternal::getDerivative(int nfwd, int nadj) {
    return f_.derivative(nfwd, nadj);
  }


  void SimpleIndefCleInternal::deepCopyMembers(
      std::map<SharedObjectNode*, SharedObject>& already_copied) {
    CleInternal::deepCopyMembers(already_copied);
  }

  SimpleIndefCleInternal* SimpleIndefCleInternal::clone() const {
    // Return a deep copy
    SimpleIndefCleInternal* node = new SimpleIndefCleInternal(st_);
    node->setOption(dictionary());
    return node;
  }


} // namespace casadi


