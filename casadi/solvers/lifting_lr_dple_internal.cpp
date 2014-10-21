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


#include "lifting_lr_dple_internal.hpp"
#include <cassert>
#include "../core/std_vector_tools.hpp"
#include "../core/matrix/matrix_tools.hpp"
#include "../core/mx/mx_tools.hpp"
#include "../core/sx/sx_tools.hpp"
#include "../core/function/mx_function.hpp"
#include "../core/function/sx_function.hpp"

#include <numeric>

INPUTSCHEME(LR_DPLEInput)
OUTPUTSCHEME(LR_DPLEOutput)

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_LRDPLESOLVER_LIFTING_EXPORT
  casadi_register_lrdplesolver_lifting(LrDpleInternal::Plugin* plugin) {
    plugin->creator = LiftingLrDpleInternal::creator;
    plugin->name = "lifting";
    plugin->doc = LiftingLrDpleInternal::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_LRDPLESOLVER_LIFTING_EXPORT casadi_load_lrdplesolver_lifting() {
    LrDpleInternal::registerPlugin(casadi_register_lrdplesolver_lifting);
  }

  LiftingLrDpleInternal::LiftingLrDpleInternal(
      const LrDpleStructure & st,
      const std::vector< std::vector<int> > &Hs) : LrDpleInternal(st, Hs) {

    // set default options
    setOption("name", "unnamed_lifting_lr_dple_solver"); // name of the function

    addOption("form",    OT_STRING,    "A",
              "The form of the lifting", "A:0|B:1");

    Adaptor::addOptions();
  }

  LiftingLrDpleInternal::~LiftingLrDpleInternal() {

  }

  void LiftingLrDpleInternal::init() {

    form_ = getOptionEnumValue("form");

    // Initialize the base classes
    LrDpleInternal::init();

    casadi_assert_message(!pos_def_,
      "pos_def option set to True: Solver only handles the indefinite case.");
    casadi_assert_message(const_dim_,
      "const_dim option set to False: Solver only handles the True case.");

    // We will construct an MXFunction to facilitate the calculation of derivatives

    MX As = MX::sym("As", input(LR_DLE_A).sparsity());
    MX Vs = MX::sym("Vs", input(LR_DLE_V).sparsity());
    MX Cs = MX::sym("Cs", input(LR_DLE_C).sparsity());
    MX Hs = MX::sym("Hs", input(LR_DLE_H).sparsity());

    n_ = A_[0].size1();

    // Chop-up the arguments
    std::vector<MX> As_ = horzsplit(As, n_);
    std::vector<MX> Vs_ = horzsplit(Vs, V_[0].size2());
    std::vector<MX> Cs_ = horzsplit(Cs, V_[0].size2());
    std::vector<MX> Hs_;
    if (with_H_) {
      Hs_ = horzsplit(Hs, Hsi_);
    }

    MX A;
    if (K_==1) {
      A = As;
    } else {
      if (form_==0) {
        MX AL = blkdiag(vector_slice(As_, range(As_.size()-1)));

        MX AL2 = horzcat(AL, Sparsity::sparse(AL.size1(), As_[0].size2()));
        MX AT = horzcat(Sparsity::sparse(As_[0].size1(), AL.size2()), As_.back());
        A = vertcat(AT, AL2);
      } else {
        MX AL = blkdiag(reverse(vector_slice(As_, range(As_.size()-1))));

        MX AL2 = horzcat(Sparsity::sparse(AL.size1(), As_[0].size2()), AL);
        MX AT = horzcat(As_.back(), Sparsity::sparse(As_[0].size1(), AL.size2()));
        A = vertcat(AL2, AT);
      }
    }

    MX V;
    MX C;

    MX H;

    if (form_==0) {
      V = blkdiag(Vs_.back(), blkdiag(vector_slice(Vs_, range(Vs_.size()-1))));
      if (with_C_) {
        C = blkdiag(Cs_.back(), blkdiag(vector_slice(Cs_, range(Cs_.size()-1))));
      }
    } else {
      V = blkdiag(blkdiag(reverse(vector_slice(Vs_, range(Vs_.size()-1)))), Vs_.back());
      if (with_C_) {
        C = blkdiag(blkdiag(reverse(vector_slice(Cs_, range(Cs_.size()-1)))), Cs_.back());
      }
    }

    if (with_H_) {
      H = blkdiag(form_==0? Hs_ : reverse(Hs_));
    }

    // Create an LrDleSolver instance
    solver_ = LrDleSolver(getOption(solvername()),
                          lrdleStruct("a", A.sparsity(),
                                      "v", V.sparsity(),
                                      "c", C.sparsity(),
                                      "h", H.sparsity()),
                          Hss_);
    if (hasSetOption(optionsname())) solver_.setOption(getOption(optionsname()));
    solver_.init();

    std::vector<MX> v_in(LR_DPLE_NUM_IN);
    v_in[LR_DLE_A] = As;
    v_in[LR_DLE_V] = Vs;
    if (with_C_) {
      v_in[LR_DLE_C] = Cs;
    }
    if (with_H_) {
      v_in[LR_DLE_H] = Hs;
    }

    std::vector<MX> Pr = solver_.call(lrdpleIn("a", A, "v", V, "c", C, "h", H));

    MX Pf = Pr[0];

    std::vector<MX> Ps = with_H_ ? diagsplit(Pf, Hsi_) : diagsplit(Pf, n_);

    if (form_==1) {
      Ps = reverse(Ps);
    }

    f_ = MXFunction(v_in, dpleOut("p", horzcat(Ps)));
    f_.setInputScheme(SCHEME_LR_DPLEInput);
    f_.setOutputScheme(SCHEME_LR_DPLEOutput);
    f_.init();

    Wrapper::checkDimensions();

  }



  void LiftingLrDpleInternal::evaluate() {
    Wrapper::evaluate();
  }

  Function LiftingLrDpleInternal::getDerivative(int nfwd, int nadj) {
    return f_.derivative(nfwd, nadj);
  }


  void LiftingLrDpleInternal::deepCopyMembers(
      std::map<SharedObjectNode*, SharedObject>& already_copied) {
    LrDpleInternal::deepCopyMembers(already_copied);
  }

  LiftingLrDpleInternal* LiftingLrDpleInternal::clone() const {
    // Return a deep copy
    LiftingLrDpleInternal* node = new LiftingLrDpleInternal(st_, Hs_);
    node->setOption(dictionary());
    return node;
  }


} // namespace casadi


