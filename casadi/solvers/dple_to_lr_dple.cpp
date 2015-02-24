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


#include "dple_to_lr_dple.hpp"
#include <cassert>
#include "../core/std_vector_tools.hpp"
#include "../core/matrix/matrix_tools.hpp"
#include "../core/mx/mx_tools.hpp"
#include "../core/sx/sx_tools.hpp"
#include "../core/function/mx_function.hpp"
#include "../core/function/sx_function.hpp"

#include <numeric>

INPUTSCHEME(LR_DPLEInput)

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_DPLESOLVER_LRDPLE_EXPORT
      casadi_register_dplesolver_lrdple(DpleInternal::Plugin* plugin) {
    plugin->creator = DpleToLrDple::creator;
    plugin->name = "lrdple";
    plugin->doc = DpleToLrDple::meta_doc.c_str();
    plugin->version = 22;
    return 0;
  }

  extern "C"
  void CASADI_DPLESOLVER_LRDPLE_EXPORT casadi_load_dplesolver_lrdple() {
    DpleInternal::registerPlugin(casadi_register_dplesolver_lrdple);
  }

  DpleToLrDple::DpleToLrDple(
         const DpleStructure& st) :
          DpleInternal(st) {

    // set default options
    setOption("name", "unnamed_dple_to_lr_dple"); // name of the function

    Adaptor::addOptions();

  }

  DpleToLrDple::~DpleToLrDple() {

  }

  void DpleToLrDple::init() {
    // Initialize the base classes
    DpleInternal::init();

    MX A = MX::sym("A", input(DPLE_A).sparsity());
    MX V = MX::sym("V", input(DPLE_V).sparsity());

    int n = A.size1();

    DMatrix C = DMatrix::eye(n);
    DMatrix H = DMatrix::eye(n);

    int K = A_.size();

    // Create an LrDpleSolver instance
    solver_ = LrDpleSolver(getOption(solvername()),
                           lrdpleStruct("a", st_[DPLE_A],
                                        "v", st_[DPLE_V],
                                        "c", std::vector<Sparsity>(K, C.sparsity()),
                                        "h", std::vector<Sparsity>(K, H.sparsity())));
    solver_.setOption("Hs", std::vector< std::vector<int> >(K, std::vector<int>(1, n)));
    solver_.setOption(getOption(optionsname()));
    solver_.init();

    std::vector<MX> Pr = solver_(
      lrdpleIn("a", A, "v", V, "c", MX(repmat(C, 1, K)), "h", MX(repmat(H, 1, K))));

    f_ = MXFunction(dpleIn("a", A, "v", V),
                    dpleOut("p", Pr[DPLE_P]));
    f_.init();

    Wrapper::checkDimensions();
  }



  void DpleToLrDple::evaluate() {
    Wrapper::evaluate();
  }

  Function DpleToLrDple::getDerForward(int nfwd) {
    return f_.derForward(nfwd);
  }

  Function DpleToLrDple::getDerReverse(int nadj) {
    return f_.derReverse(nadj);
  }

  void DpleToLrDple::deepCopyMembers(
      std::map<SharedObjectNode*, SharedObject>& already_copied) {
    DpleInternal::deepCopyMembers(already_copied);
  }

  DpleToLrDple* DpleToLrDple::clone() const {
    // Return a deep copy
    DpleToLrDple* node = new DpleToLrDple(st_);
    node->setOption(dictionary());
    return node;
  }


} // namespace casadi


