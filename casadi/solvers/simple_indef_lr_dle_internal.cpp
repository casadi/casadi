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


#include "simple_indef_dle_internal.hpp"
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
  int CASADI_DLESOLVER_SIMPLE_EXPORT
  casadi_register_dlesolver_simple(DleInternal::Plugin* plugin) {
    plugin->creator = SimpleIndefDleInternal::creator;
    plugin->name = "simple";
    plugin->doc = SimpleIndefDleInternal::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_DLESOLVER_SIMPLE_EXPORT casadi_load_dlesolver_simple() {
    DleInternal::registerPlugin(casadi_register_dlesolver_simple);
  }

  SimpleIndefDleInternal::SimpleIndefDleInternal(
      const DleStructure& st, const std::vector<int> &Hs) : DleInternal(st,Hs) {

    // set default options
    setOption("name", "unnamed_simple_indef_dle_solver"); // name of the function

    addOption("compressed_solve",         OT_BOOLEAN, true,
              "When a system with sparse rhs arises, compress to"
              "a smaller system with dense rhs.");
    addOption("linear_solver",            OT_STRING, GenericType(),
              "User-defined linear solver class. Needed for sensitivities.");
    addOption("linear_solver_options",    OT_DICTIONARY,   GenericType(),
              "Options to be passed to the linear solver.");

  }

  SimpleIndefDleInternal::~SimpleIndefDleInternal() {

  }

  void SimpleIndefDleInternal::init() {

    DleInternal::init();

    casadi_assert_message(!pos_def_,
      "pos_def option set to True: Solver only handles the indefinite case.");


    n_ = A_.size1();

    MX As = MX::sym("A", A_);
    MX Vs = MX::sym("V", V_);
    MX Cs = MX::sym("C", C_);
    MX Hs = MX::sym("H", H_);

    MX Vss = (Vs+Vs.T())/2;
    if (with_C_) Vss = mul(mul(Cs, Vss), Cs.T());

    MX A_total = DMatrix::eye(n_*n_) - kron(As,As);

    // Should be treated by solve node
    MX Pf = solve(A_total, vec(Vss), getOption("linear_solver"));

    std::vector<MX> v_in;
    v_in.push_back(As);
    v_in.push_back(Vs);
    v_in.push_back(Cs);

    MX P = reshape(Pf,n_,n_);

    std::vector<MX> HPH;

    if (with_H_) {
      std::vector<MX> H = horzsplit(Hs,Hi_);

      for (int k=0;k<H.size();++k) {
        HPH.push_back(mul(H[k].T(),mul(P,H[k])));
      }
    }

    std::vector<MX> dle_in(DLE_NUM_IN);
    dle_in[DLE_A] = As;
    dle_in[DLE_V] = Vs;
    if (with_C_) dle_in[DLE_C] = Cs;
    if (with_H_) dle_in[DLE_H] = Hs;

    f_ = MXFunction(dle_in,dleOut("p",with_H_? blkdiag(HPH) : P(output().sparsity())));

    f_.init();

    casadi_assert(getNumOutputs()==f_.getNumOutputs());
    for (int i=0;i<getNumInputs();++i) {
      casadi_assert_message(input(i).sparsity()==f_.input(i).sparsity(),
        "Sparsity mismatch for input " << i << ":" <<
        input(i).dimString() << " <-> " << f_.input(i).dimString() << ".");
    }
    for (int i=0;i<getNumOutputs();++i) {
      casadi_assert_message(output(i).sparsity()==f_.output(i).sparsity(),
        "Sparsity mismatch for output " << i << ":" <<
        output(i).dimString() << " <-> " << f_.output(i).dimString() << ".");
    }
  }



  void SimpleIndefDleInternal::evaluate() {
    std::cout << "eval" << std::endl;
    input(DLE_A).printDense();
    for (int i=0;i<getNumInputs();++i) {
      std::copy(input(i).begin(), input(i).end(), f_.input(i).begin());
    }
    f_.evaluate();
    for (int i=0;i<getNumOutputs();++i) {
      std::copy(f_.output(i).begin(), f_.output(i).end(), output(i).begin());
    }
  }

  Function SimpleIndefDleInternal::getDerivative(int nfwd, int nadj) {
    return f_.derivative(nfwd, nadj);
  }


  void SimpleIndefDleInternal::deepCopyMembers(
      std::map<SharedObjectNode*, SharedObject>& already_copied) {
    DleInternal::deepCopyMembers(already_copied);
  }

  SimpleIndefDleInternal* SimpleIndefDleInternal::clone() const {
    // Return a deep copy
    SimpleIndefDleInternal* node = new SimpleIndefDleInternal(st_, Hs_);
    node->setOption(dictionary());
    return node;
  }


} // namespace casadi


