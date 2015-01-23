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


#include "lr_dle_to_dle.hpp"
#include <cassert>
#include "../core/std_vector_tools.hpp"
#include "../core/matrix/matrix_tools.hpp"
#include "../core/mx/mx_tools.hpp"
#include "../core/sx/sx_tools.hpp"
#include "../core/function/mx_function.hpp"
#include "../core/function/sx_function.hpp"

#include <numeric>

INPUTSCHEME(LR_DLEInput)

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_LRDLESOLVER_DLE_EXPORT
  casadi_register_lrdlesolver_dle(LrDleInternal::Plugin* plugin) {
    plugin->creator = LrDleToDle::creator;
    plugin->name = "dle";
    plugin->doc = LrDleToDle::meta_doc.c_str();
    plugin->version = 22;
    plugin->adaptorHasPlugin = DleSolver::hasPlugin;
    return 0;
  }

  extern "C"
  void CASADI_LRDLESOLVER_DLE_EXPORT casadi_load_lrdlesolver_dle() {
    LrDleInternal::registerPlugin(casadi_register_lrdlesolver_dle);
  }

  LrDleToDle::LrDleToDle(
         const LrDleStructure& st) : LrDleInternal(st) {

    // set default options
    setOption("name", "unnamed_lr_dle_to_dle"); // name of the function

    Adaptor::addOptions();
  }

  LrDleToDle::~LrDleToDle() {

  }

  void LrDleToDle::init() {
    // Initialize the base classes
    LrDleInternal::init();

    MX H = MX::sym("H", H_);
    MX A = MX::sym("A", A_);
    MX C = MX::sym("C", C_);
    MX V = MX::sym("V", V_);

    MX Vs = (V+V.T())/2;

    MX CVC = mul(C, mul(V, C.T()));

    // Create an DleSolver instance
    solver_ = DleSolver(getOption(solvername()),
                        dleStruct("a", A_, "v", CVC.sparsity()));
    if (hasSetOption(optionsname())) solver_.setOption(getOption(optionsname()));
    solver_.init();

    std::vector<MX> Pr = solver_.call(dleIn("a", A, "v", CVC));
    MX P = Pr[DLE_P];

    std::vector<MX> HPH(Hs_.size(), 0);
    std::vector<MX> Hs = horzsplit(H, Hi_);
    MX out = 0;

    for (int k=0;k<Hs.size();++k) {
      HPH[k] = mul(Hs[k].T(), mul(P, Hs[k]));
    }

    f_ = MXFunction(lrdpleIn("a", A, "v", V, "c", C, "h", H),
                    lrdleOut("y", diagcat(HPH)));
    f_.init();

    Wrapper::checkDimensions();

  }

  void LrDleToDle::evaluate() {
    Wrapper::evaluate();
  }

  Function LrDleToDle::getDerivative(int nfwd, int nadj) {
    return f_.derivative(nfwd, nadj);
  }


  void LrDleToDle::deepCopyMembers(
      std::map<SharedObjectNode*, SharedObject>& already_copied) {
    LrDleInternal::deepCopyMembers(already_copied);
  }

  LrDleToDle* LrDleToDle::clone() const {
    // Return a deep copy
    LrDleToDle* node = new LrDleToDle(st_);
    node->setOption(dictionary());
    return node;
  }


} // namespace casadi


