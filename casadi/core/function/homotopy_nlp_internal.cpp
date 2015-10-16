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


#include "homotopy_nlp_internal.hpp"
#include "mx_function.hpp"

INPUTSCHEME(NlpSolverInput)
OUTPUTSCHEME(NlpSolverOutput)

using namespace std;
namespace casadi {

  HomotopyNLPInternal::HomotopyNLPInternal(const std::string& name,
                                           const Function& hnlp)
  : FunctionInternal(name), hnlp_(hnlp) {

    addOption("expand",             OT_BOOLEAN,  false,
              "Expand the NLP function in terms of scalar operations, i.e. MX->SX");

    // Enable string notation for IO
    ischeme_ = IOScheme(SCHEME_NlpSolverInput);
    oscheme_ = IOScheme(SCHEME_NlpSolverOutput);
  }

  HomotopyNLPInternal::~HomotopyNLPInternal() {
  }

  Sparsity HomotopyNLPInternal::get_sparsity_in(int ind) const {
    switch (static_cast<NlpSolverInput>(ind)) {
    case NLP_SOLVER_X0:
    case NLP_SOLVER_LBX:
    case NLP_SOLVER_UBX:
    case NLP_SOLVER_LAM_X0:
      return get_sparsity_out(NLP_SOLVER_X);
    case NLP_SOLVER_LBG:
    case NLP_SOLVER_UBG:
    case NLP_SOLVER_LAM_G0:
      return get_sparsity_out(NLP_SOLVER_G);
    case NLP_SOLVER_P:
      return get_sparsity_out(NLP_SOLVER_P);
    case NLP_SOLVER_NUM_IN:
      return Sparsity();
    }
  }

  Sparsity HomotopyNLPInternal::get_sparsity_out(int ind) const {
    switch (static_cast<NlpSolverOutput>(ind)) {
    case NLP_SOLVER_F:
      return Sparsity::scalar();
    case NLP_SOLVER_X:
    case NLP_SOLVER_LAM_X:
      return hnlp_.sparsity_in(HNL_X);
    case NLP_SOLVER_LAM_G:
    case NLP_SOLVER_G:
      return hnlp_.sparsity_out(NL_G);
    case NLP_SOLVER_LAM_P:
      return hnlp_.sparsity_in(HNL_P);
    case NLP_SOLVER_NUM_OUT:
      return Sparsity();
    }
  }

  void HomotopyNLPInternal::init() {
    // Initialize the Homotopy NLP
    hnlp_.init();

    casadi_assert_message(hnlp_.n_in()==HNL_NUM_IN,
                          "The HNLP function must have exactly three input");
    casadi_assert_message(hnlp_.n_out()==NL_NUM_OUT,
                          "The HNLP function must have exactly two outputs");

    // Sparsity patterns
    const Sparsity& x_sparsity = hnlp_.input(HNL_X).sparsity();
    const Sparsity& p_sparsity = hnlp_.input(HNL_P).sparsity();
    const Sparsity& g_sparsity = hnlp_.output(NL_G).sparsity();

    // Get dimensions
    nx_ = x_sparsity.nnz();
    np_ = p_sparsity.nnz();
    ng_ = g_sparsity.nnz();

    // Allocate space for inputs
    ibuf_.resize(NLP_SOLVER_NUM_IN);
    input(NLP_SOLVER_X0)       =  DMatrix::zeros(x_sparsity);
    input(NLP_SOLVER_LBX)      = -DMatrix::inf(x_sparsity);
    input(NLP_SOLVER_UBX)      =  DMatrix::inf(x_sparsity);
    input(NLP_SOLVER_LBG)      = -DMatrix::inf(g_sparsity);
    input(NLP_SOLVER_UBG)      =  DMatrix::inf(g_sparsity);
    input(NLP_SOLVER_LAM_X0)   =  DMatrix::zeros(x_sparsity);
    input(NLP_SOLVER_LAM_G0)   =  DMatrix::zeros(g_sparsity);
    input(NLP_SOLVER_P)        =  DMatrix::zeros(p_sparsity);

    // Allocate space for outputs
    obuf_.resize(NLP_SOLVER_NUM_OUT);
    output(NLP_SOLVER_X)       = DMatrix::zeros(x_sparsity);
    output(NLP_SOLVER_F)       = DMatrix::zeros(1);
    output(NLP_SOLVER_LAM_X)   = DMatrix::zeros(x_sparsity);
    output(NLP_SOLVER_LAM_G)   = DMatrix::zeros(g_sparsity);
    output(NLP_SOLVER_LAM_P)   = DMatrix::zeros(p_sparsity);
    output(NLP_SOLVER_G)       = DMatrix::zeros(g_sparsity);

    // Call the initialization method of the base class
    FunctionInternal::init();

    // Find out if we are to expand the NLP in terms of scalar operations
    bool expand = getOption("expand");
    if (expand) {
      log("Expanding NLP in scalar operations");

      // Cast to MXFunction
      MXFunction hnlp_mx = shared_cast<MXFunction>(hnlp_);
      if (hnlp_mx.isNull()) {
        casadi_warning("Cannot expand NLP as it is not an MXFunction");
      } else {
        hnlp_ = SX::fun(hnlp_mx.name(), hnlp_mx);
        hnlp_.copyOptions(hnlp_mx, true);
        hnlp_.init();
      }
    }
  }

  std::map<std::string, HomotopyNLPInternal::Plugin> HomotopyNLPInternal::solvers_;

  const std::string HomotopyNLPInternal::infix_ = "homotopynlpsolver";


} // namespace casadi
