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
  int CASADI_LRDPLESOLVER_DPLE_EXPORT
  casadi_register_lrdplesolver_dple(LrDpleInternal::Plugin* plugin) {
    plugin->creator = LrDpleToDple::creator;
    plugin->name = "dple";
    plugin->doc = LrDpleToDple::meta_doc.c_str();
    plugin->version = 22;
    plugin->adaptorHasPlugin = DpleSolver::hasPlugin;
    return 0;
  }

  extern "C"
  void CASADI_LRDPLESOLVER_DPLE_EXPORT casadi_load_lrdplesolver_dple() {
    LrDpleInternal::registerPlugin(casadi_register_lrdplesolver_dple);
  }

  LrDpleToDple::LrDpleToDple(
         const LrDpleStructure& st) :
          LrDpleInternal(st) {

    // set default options
    setOption("name", "unnamed_lr_dple_to_dple"); // name of the function

    Adaptor::addOptions();

  }

  LrDpleToDple::~LrDpleToDple() {

  }

  void LrDpleToDple::init() {
    // Initialize the base classes
    LrDpleInternal::init();

    MX As = MX::sym("As", input(LR_DPLE_A).sparsity());
    MX Vs = MX::sym("Vs", input(LR_DPLE_V).sparsity());
    MX Cs = MX::sym("Cs", input(LR_DPLE_C).sparsity());
    MX Hs = MX::sym("Hs", input(LR_DPLE_H).sparsity());

    int n_ = A_[0].size1();

    // Chop-up the arguments
    std::vector<MX> As_ = horzsplit(As, n_);
    std::vector<MX> Vs_ = horzsplit(Vs, V_[0].size2());
    std::vector<MX> Cs_ = horzsplit(Cs, V_[0].size2());
    std::vector<MX> Hss_ = horzsplit(Hs, Hsi_);

    std::vector<MX> V_(Vs_.size());

    for (int k=0;k<V_.size();++k) {
      V_[k] = mul(Cs_[k], mul(Vs_[k], Cs_[k].T()));
    }

    std::vector<Sparsity> Vsp(Vs_.size());
    for (int k=0;k<V_.size();++k) {
      Vsp[k] = V_[k].sparsity();
    }

    // Create an dplesolver instance
    solver_ = DpleSolver(getOption(solvername()),
                         dpleStruct("a", A_, "v", Vsp));
    if (hasSetOption(optionsname())) solver_.setOption(getOption(optionsname()));
    solver_.init();

    std::vector<MX> Pr = solver_.call(dpleIn("a", horzcat(As_), "v", horzcat(V_)));
    std::vector<MX> Ps_ = horzsplit(Pr[DPLE_P], n_);

    std::vector<MX> HPH(K_);

    for (int k=0;k<K_;++k) {
      std::vector<MX> hph = horzsplit(Hss_[k], cumsum0(Hs_[k]));

      for (int kk=0;kk<hph.size();++kk) {
        hph[kk] = mul(hph[kk].T(), mul(Ps_[k], hph[kk]));
      }
      HPH[k] = diagcat(hph);
    }


    f_ = MXFunction(lrdpleIn("a", As, "v", Vs, "c", Cs, "h", Hs),
                    lrdpleOut("y", horzcat(HPH)));
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
    LrDpleInternal::deepCopyMembers(already_copied);
  }

  LrDpleToDple* LrDpleToDple::clone() const {
    // Return a deep copy
    LrDpleToDple* node = new LrDpleToDple(st_);
    node->setOption(dictionary());
    return node;
  }

} // namespace casadi


