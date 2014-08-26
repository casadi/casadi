/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#include "condensing_indef_dple_internal.hpp"
#include <cassert>
#include "../core/std_vector_tools.hpp"
#include "../core/matrix/matrix_tools.hpp"
#include "../core/mx/mx_tools.hpp"
#include "../core/sx/sx_tools.hpp"
#include "../core/function/mx_function.hpp"
#include "../core/function/sx_function.hpp"

#include <numeric>

INPUTSCHEME(DPLEInput)
OUTPUTSCHEME(DPLEOutput)

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_DPLESOLVER_CONDENSING_EXPORT
  casadi_register_dplesolver_condensing(DpleInternal::Plugin* plugin) {
    plugin->creator = CondensingIndefDpleInternal::creator;
    plugin->name = "condensing";
    plugin->doc = CondensingIndefDpleInternal::meta_doc.c_str();
    plugin->version = 20;
    return 0;
  }

  extern "C"
  void CASADI_DPLESOLVER_CONDENSING_EXPORT casadi_load_dplesolver_condensing() {
    DpleInternal::registerPlugin(casadi_register_dplesolver_condensing);
  }

  CondensingIndefDpleInternal::CondensingIndefDpleInternal(
      const std::vector< Sparsity > & A,
      const std::vector< Sparsity > &V) : DpleInternal(A, V) {

    // set default options
    setOption("name", "unnamed_condensing_indef_dple_solver"); // name of the function

    addOption("dle_solver",            OT_STRING, GenericType(),
              "User-defined Dle solver class. Needed for sensitivities.");
    addOption("dle_solver_options",    OT_DICTIONARY,   GenericType(),
              "Options to be passed to the Dle solver.");

  }

  CondensingIndefDpleInternal::~CondensingIndefDpleInternal() {

  }

  void CondensingIndefDpleInternal::init() {

    DpleInternal::init();

    casadi_assert_message(!pos_def_,
      "pos_def option set to True: Solver only handles the indefinite case.");
    casadi_assert_message(const_dim_,
      "const_dim option set to False: Solver only handles the True case.");

    n_ = A_[0].size1();


    MX As = MX::sym("A", n_, K_*n_);
    MX Vs = MX::sym("V", n_, K_*n_);

    std::vector< MX > Vss = horzsplit(Vs, n_);
    std::vector< MX > Ass = horzsplit(As, n_);

    for (int k=0;k<K_;++k) {
      Vss[k] = (Vss[k]+Vss[k].T())/2;
    }
    
    MX R = MX::zeros(n_,n_);
    
    for (int k=0;k<K_;++k) {
      R = mul(mul(Ass[k],R),Ass[k].T()) + Vss[k];
    }
    
    std::vector< MX > Assr(K_);    
    std::reverse_copy (Ass.begin(), Ass.end(), Assr.begin());

    MX Ap = mul(Assr);
    
    // Create an dlesolver instance
    std::string dlesolver_name = getOption("dle_solver");
    dlesolver_ = DleSolver(dlesolver_name, Ap.sparsity(), R.sparsity());

    if (hasSetOption("dle_solver_options")) {
      dlesolver_.setOption(getOption("dle_solver_options"));
    }

    // Initialize the NLP solver
    dlesolver_.init();

    std::vector<MX> v_in;
    v_in.push_back(As);
    v_in.push_back(Vs);
    
    std::vector<MX> Pr = dlesolver_.call(dpleIn("a",Ap,"v",R));
    
    std::vector<MX> Ps(K_);
    Ps[0] = Pr[0];
    
    for (int k=0;k<K_-1;++k) {
      Ps[k+1] = mul(mul(Ass[k],Ps[k]),Ass[k].T()) + Vss[k];
    }

    f_ = MXFunction(v_in, horzcat(Ps));
    f_.setInputScheme(SCHEME_DPLEInput);
    f_.setOutputScheme(SCHEME_DPLEOutput);
    f_.init();
  }



  void CondensingIndefDpleInternal::evaluate() {
    for (int i=0;i<getNumInputs();++i) {
      std::copy(input(i).begin(), input(i).end(), f_.input(i).begin());
    }
    f_.evaluate();
    for (int i=0;i<getNumOutputs();++i) {
      std::copy(f_.output(i).begin(), f_.output(i).end(), output(i).begin());
    }
  }

  Function CondensingIndefDpleInternal::getDerivative(int nfwd, int nadj) {
    return f_.derivative(nfwd, nadj);
  }


  void CondensingIndefDpleInternal::deepCopyMembers(
      std::map<SharedObjectNode*, SharedObject>& already_copied) {
    DpleInternal::deepCopyMembers(already_copied);
  }

  CondensingIndefDpleInternal* CondensingIndefDpleInternal::clone() const {
    // Return a deep copy
    CondensingIndefDpleInternal* node = new CondensingIndefDpleInternal(A_, V_);
    node->setOption(dictionary());
    return node;
  }


} // namespace casadi


