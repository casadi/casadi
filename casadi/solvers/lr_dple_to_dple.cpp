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


#include "lr_dple_to_dple.hpp"
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
    plugin->creator = LrDpleToDple::creator;
    plugin->name = "lrdple";
    plugin->doc = LrDpleToDple::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_DPLESOLVER_LRDPLE_EXPORT casadi_load_dplesolver_lrdple() {
    DpleInternal::registerPlugin(casadi_register_dplesolver_lrdple);
  }

  LrDpleToDple::LrDpleToDple(
         const DpleStructure& st) :
          DpleInternal(st) {

    // set default options
    setOption("name", "unnamed_lr_dple_to_dple"); // name of the function

    Adaptor::addOptions();

  }

  LrDpleToDple::~LrDpleToDple() {

  }

  void LrDpleToDple::init() {
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
                                        "h", std::vector<Sparsity>(K, H.sparsity())),
                           std::vector< std::vector<int> >(K, std::vector<int>(1, n)));
    solver_.setOption(getOption(optionsname()));
    solver_.init();

    std::vector<MX> Pr = solver_.call(
      lrdpleIn("a", A, "v", V, "c", MX(repmat(C, 1, K)), "h", MX(repmat(H, 1, K))));

    f_ = MXFunction(dpleIn("a", A, "v", V),
                    dpleOut("p", Pr[DPLE_P]));
    f_.init();

    Wrapper::checkDimensions();
  }



  void LrDpleToDple::evaluate() {
    Wrapper::evaluate();
  }

  Function LrDpleToDple::getDerivative(int nfwd, int nadj) {
    return f_.derivative(nfwd, nadj);
  }


  void LrDpleToDple::deepCopyMembers(
      std::map<SharedObjectNode*, SharedObject>& already_copied) {
    DpleInternal::deepCopyMembers(already_copied);
  }

  LrDpleToDple* LrDpleToDple::clone() const {
    // Return a deep copy
    LrDpleToDple* node = new LrDpleToDple(st_);
    node->setOption(dictionary());
    return node;
  }


} // namespace casadi


