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


#include "qp_to_nlp.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_QPSOL_NLPSOL_EXPORT
  casadi_register_qpsol_nlpsol(Qpsol::Plugin* plugin) {
    plugin->creator = QpToNlp::creator;
    plugin->name = "nlpsol";
    plugin->doc = QpToNlp::meta_doc.c_str();
    plugin->version = 23;
    return 0;
  }

  extern "C"
  void CASADI_QPSOL_NLPSOL_EXPORT casadi_load_qpsol_nlpsol() {
    Qpsol::registerPlugin(casadi_register_qpsol_nlpsol);
  }

  QpToNlp::QpToNlp(const std::string& name, const std::map<std::string, Sparsity> &st)
    : Qpsol(name, st) {

    addOption("nlpsol", OT_STRING,
              "Name of solver.");
    addOption("nlpsol_options", OT_DICT,
              "Options to be passed to solver.");
  }

  QpToNlp::~QpToNlp() {
  }

  void QpToNlp::init(const Dict& opts) {
    // Initialize the base classes
    Qpsol::init(opts);

    // Default options
    string nlpsol_plugin;
    Dict nlpsol_options;

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="nlpsol") {
        nlpsol_plugin = op.second.to_string();
      } else if (op.first=="nlpsol_options") {
        nlpsol_options = op.second;
      }
    }

    // Create a symbolic matrix for the decision variables
    SX X = SX::sym("X", n_, 1);

    // Parameters to the problem
    SX H = SX::sym("H", sparsity_in(QPSOL_H));
    SX G = SX::sym("G", sparsity_in(QPSOL_G));
    SX A = SX::sym("A", sparsity_in(QPSOL_A));

    // Put parameters in a vector
    std::vector<SX> par;
    par.push_back(H.data());
    par.push_back(G.data());
    par.push_back(A.data());

    // The nlp looks exactly like a mathematical description of the NLP
    SXDict nlp = {{"x", X}, {"p", vertcat(par)},
                  {"f", mtimes(G.T(), X) + 0.5*mtimes(mtimes(X.T(), H), X)},
                  {"g", mtimes(A, X)}};

    // Create an Nlpsol instance
    casadi_assert_message(!nlpsol_plugin.empty(), "'nlpsol' option has not been set");
    solver_ = nlpsol("nlpsol", nlpsol_plugin, nlp, nlpsol_options);
    alloc(solver_);

    // Allocate storage for NLP solver  parameters
    alloc_w(solver_.nnz_in(NLPSOL_P), true);
  }

  void QpToNlp::
  eval(Memory& mem, const double** arg, double** res, int* iw, double* w) const {
    // Inputs
    const double *h_, *g_, *a_, *lba_, *uba_, *lbx_, *ubx_, *x0_;
    // Outputs
    double *x_, *f_, *lam_a_, *lam_x_;

    // Get input pointers
    h_ = arg[QPSOL_H];
    g_ = arg[QPSOL_G];
    a_ = arg[QPSOL_A];
    lba_ = arg[QPSOL_LBA];
    uba_ = arg[QPSOL_UBA];
    lbx_ = arg[QPSOL_LBX];
    ubx_ = arg[QPSOL_UBX];
    x0_ = arg[QPSOL_X0];

    // Get output pointers
    x_ = res[QPSOL_X];
    f_ = res[QPSOL_COST];
    lam_a_ = res[QPSOL_LAM_A];
    lam_x_ = res[QPSOL_LAM_X];

    // Buffers for calling the NLP solver
    const double** arg1 = arg + n_in();
    double** res1 = res + n_out();
    fill_n(arg1, static_cast<int>(NLPSOL_NUM_IN), nullptr);
    fill_n(res1, static_cast<int>(NLPSOL_NUM_OUT), nullptr);

    // NLP inputs
    arg1[NLPSOL_X0] = x0_;
    arg1[NLPSOL_LBG] = lba_;
    arg1[NLPSOL_UBG] = uba_;
    arg1[NLPSOL_LBX] = lbx_;
    arg1[NLPSOL_UBX] = ubx_;

    // NLP parameters
    arg1[NLPSOL_P] = w;

    // Quadratic term
    int nh = nnz_in(QPSOL_H);
    if (h_) {
      copy_n(h_, nh, w);
    } else {
      fill_n(w, nh, 0);
    }
    w += nh;

    // Linear objective term
    int ng = nnz_in(QPSOL_G);
    if (g_) {
      copy_n(g_, ng, w);
    } else {
      fill_n(w, ng, 0);
    }
    w += ng;

    // Linear constraints
    int na = nnz_in(QPSOL_A);
    if (a_) {
      copy_n(a_, na, w);
    } else {
      fill_n(w, na, 0);
    }
    w += na;

    // Solution
    res1[NLPSOL_X] = x_;
    res1[NLPSOL_F] = f_;
    res1[NLPSOL_LAM_X] = lam_x_;
    res1[NLPSOL_LAM_G] = lam_a_;

    // Solve the NLP
    solver_(arg1, res1, iw, w, 0);
  }

} // namespace casadi
