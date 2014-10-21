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


#include "fixed_smith_lr_dle_internal.hpp"
#include <cassert>
#include "../core/std_vector_tools.hpp"
#include "../core/matrix/matrix_tools.hpp"
#include "../core/mx/mx_tools.hpp"
#include "../core/sx/sx_tools.hpp"
#include "../core/function/mx_function.hpp"
#include "../core/function/sx_function.hpp"
#include <iomanip>

#include <numeric>

INPUTSCHEME(LR_DLEInput)

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_LRDLESOLVER_FIXED_SMITH_EXPORT
  casadi_register_lrdlesolver_fixed_smith(LrDleInternal::Plugin* plugin) {
    plugin->creator = FixedSmithLrDleInternal::creator;
    plugin->name = "fixed_smith";
    plugin->doc = FixedSmithLrDleInternal::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_LRDLESOLVER_FIXED_SMITH_EXPORT casadi_load_lrdlesolver_fixed_smith() {
    LrDleInternal::registerPlugin(casadi_register_lrdlesolver_fixed_smith);
  }

  FixedSmithLrDleInternal::FixedSmithLrDleInternal(
      const LrDleStructure& st, const std::vector<int> &Hs) :
      LrDleInternal(st, Hs) {

    // set default options
    setOption("name", "unnamed_fixed_smith_indef_dle_solver"); // name of the function

    addOption("iter", OT_INTEGER, 100,   "Number of Smith iterations");


  }

  FixedSmithLrDleInternal::~FixedSmithLrDleInternal() {

  }

  void FixedSmithLrDleInternal::init() {
    iter_  = getOption("iter");

    LrDleInternal::init();

    casadi_assert_message(!pos_def_,
      "pos_def option set to True: Solver only handles the indefinite case.");

    MX H = MX::sym("H", H_);
    MX A = MX::sym("A", A_);
    MX C = MX::sym("C", C_);
    MX V = MX::sym("V", V_);

    MX Vs = (V+V.T())/2;

    MX D = with_C_ ? C : DMatrix::eye(A_.size1());


    std::vector<MX> HPH(Hs_.size(), 0);
    std::vector<MX> Hs = with_H_? horzsplit(H, Hi_) : std::vector<MX>();
    MX out = 0;

    for (int i=0;i<iter_;++i) {
      if (with_H_) {
        for (int k=0;k<Hs.size();++k) {
          HPH[k]+= quad_form(mul(D.T(), Hs[k]), Vs);
        }
      } else {
        out += mul(D, mul(Vs, D.T()));
      }
      D = mul(A, D);
    }

    std::vector<MX> dle_in(LR_DLE_NUM_IN);
    dle_in[LR_DLE_A] = A;
    dle_in[LR_DLE_V] = V;
    if (with_C_) dle_in[LR_DLE_C] = C;
    if (with_H_) dle_in[LR_DLE_H] = H;

    f_ = MXFunction(dle_in, lrdleOut("y", with_H_? blkdiag(HPH): out));
    f_.init();

    Wrapper::checkDimensions();

  }

  void FixedSmithLrDleInternal::evaluate() {
    Wrapper::evaluate();
  }

  Function FixedSmithLrDleInternal::getDerivative(int nfwd, int nadj) {
    return f_.derivative(nfwd, nadj);
  }

  void FixedSmithLrDleInternal::deepCopyMembers(
      std::map<SharedObjectNode*, SharedObject>& already_copied) {
    LrDleInternal::deepCopyMembers(already_copied);
  }

  FixedSmithLrDleInternal* FixedSmithLrDleInternal::clone() const {
    // Return a deep copy
    FixedSmithLrDleInternal* node = new FixedSmithLrDleInternal(st_, Hs_);
    node->setOption(dictionary());
    return node;
  }


} // namespace casadi


