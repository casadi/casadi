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
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_DPLESOLVER_CONDENSING_EXPORT casadi_load_dplesolver_condensing() {
    DpleInternal::registerPlugin(casadi_register_dplesolver_condensing);
  }

  CondensingIndefDpleInternal::CondensingIndefDpleInternal(
      const DpleStructure& st) : DpleInternal(st) {

    // set default options
    setOption("name", "unnamed_condensing_indef_dple_solver"); // name of the function

    Adaptor::addOptions();

  }

  CondensingIndefDpleInternal::~CondensingIndefDpleInternal() {

  }

  void CondensingIndefDpleInternal::init() {
    // Initialize the base classes
    DpleInternal::init();

    casadi_assert_message(!pos_def_,
      "pos_def option set to True: Solver only handles the indefinite case.");
    casadi_assert_message(const_dim_,
      "const_dim option set to False: Solver only handles the True case.");

    n_ = A_[0].size1();


    MX As = MX::sym("A", horzcat(A_));
    MX Vs = MX::sym("V", horzcat(V_));

    std::vector< MX > Vss = horzsplit(Vs, n_);
    std::vector< MX > Ass = horzsplit(As, n_);

    for (int k=0;k<K_;++k) {
      Vss[k] = (Vss[k]+Vss[k].T())/2;
    }

    MX R = MX::zeros(n_, n_);

    for (int k=0;k<K_;++k) {
      R = mul(mul(Ass[k], R), Ass[k].T()) + Vss[k];
    }

    std::vector< MX > Assr(K_);
    std::reverse_copy(Ass.begin(), Ass.end(), Assr.begin());

    MX Ap = mul(Assr);

    // Create an dlesolver instance
    solver_ = DleSolver(getOption(solvername()), dleStruct("a", Ap.sparsity(), "v", R.sparsity()));
    solver_.setOption(getOption(optionsname()));

    // Initialize the NLP solver
    solver_.init();

    std::vector<MX> Pr = solver_.call(dpleIn("a", Ap, "v", R));

    std::vector<MX> Ps(K_);
    Ps[0] = Pr[0];

    for (int k=0;k<K_-1;++k) {
      Ps[k+1] = mul(mul(Ass[k], Ps[k]), Ass[k].T()) + Vss[k];
    }

    f_ = MXFunction(dpleIn("a", As, "v", Vs), dpleOut("p", horzcat(Ps)));
    f_.init();

    Wrapper::checkDimensions();

  }



  void CondensingIndefDpleInternal::evaluate() {
    Wrapper::evaluate();
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
    CondensingIndefDpleInternal* node = new CondensingIndefDpleInternal(st_);
    node->setOption(dictionary());
    return node;
  }


} // namespace casadi


