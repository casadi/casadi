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


#include "implicit_to_nlp.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_ROOTFINDER_NLPSOL_EXPORT
  casadi_register_rootfinder_nlpsol(Rootfinder::Plugin* plugin) {
    plugin->creator = ImplicitToNlp::creator;
    plugin->name = "nlpsol";
    plugin->doc = ImplicitToNlp::meta_doc.c_str();
    plugin->version = 23;
    plugin->adaptorHasPlugin = Function::has_nlpsol;
    return 0;
  }

  extern "C"
  void CASADI_ROOTFINDER_NLPSOL_EXPORT casadi_load_rootfinder_nlpsol() {
    Rootfinder::registerPlugin(casadi_register_rootfinder_nlpsol);
  }

  ImplicitToNlp::ImplicitToNlp(const std::string& name, const Function& f)
    : Rootfinder(name, f) {

    addOption("nlpsol", OT_STRING, GenericType(), "Name of solver.");
    addOption("nlpsol_options", OT_DICT,  Dict(), "Options to be passed to solver.");
  }

  ImplicitToNlp::~ImplicitToNlp() {
  }

  void ImplicitToNlp::init() {
    // Call the base class initializer
    Rootfinder::init();

    // Free variable in the NLP
    MX u = MX::sym("u", input(iin_).sparsity());

    // So that we can pass it on to createParent
    std::vector<MX> inputs;
    for (int i=0; i<n_in(); ++i) {
      if (i!=iin_) {
        stringstream ss;
        ss << "p" << i;
        inputs.push_back(MX::sym(ss.str(), input(i).sparsity()));
      }
    }
    MX p = veccat(inputs);

    // Dummy NLP objective
    MX nlp_f = 0;

    // NLP constraints
    std::vector< MX > args_call(n_in());
    args_call[iin_] = u;
    for (int i=0, i2=0; i<n_in(); ++i)
      if (i!=iin_) args_call[i] = inputs[i2++];
    MX nlp_g = f_(args_call).at(iout_);

    // We're going to use two-argument objective and constraints to allow the use of parameters
    MXDict nlp = {{"x", u}, {"p", p}, {"f", nlp_f}, {"g", nlp_g}};

    Dict options;
    if (hasSetOption("nlpsol_options")) options = option("nlpsol_options");
    // Create an Nlpsol instance
    solver_ = nlpsol("nlpsol", option("nlpsol"), nlp, options);
    alloc(solver_);

    // Allocate storage for variable bounds
    alloc_w(n_, true); // lbx
    alloc_w(n_, true); // ubx

    // Allocate storage for NLP solver parameters
    alloc_w(f_.nnz_in() - nnz_in(iin_), true);

    // Allocate storage for NLP primal solution
    alloc_w(n_, true);
  }

  void ImplicitToNlp::eval(const double** arg, double** res, int* iw, double* w, void* mem) {
    // Buffers for calling the NLP solver
    const double** arg1 = arg + n_in();
    double** res1 = res + n_out();
    fill(arg1, arg1+NLPSOL_NUM_IN, nullptr);
    fill(res1, res1+NLPSOL_NUM_OUT, nullptr);

    // Initial guess
    arg1[NLPSOL_X] = arg[iin_];

    // Nonlinear bounds
    arg1[NLPSOL_LBG] = 0;
    arg1[NLPSOL_UBG] = 0;

    // Variable bounds
    double* lbx = w; w += n_;
    fill_n(lbx, n_, -std::numeric_limits<double>::infinity());
    arg1[NLPSOL_LBX] = lbx;
    double* ubx = w; w += n_;
    fill_n(ubx, n_,  std::numeric_limits<double>::infinity());
    arg1[NLPSOL_UBX] = ubx;
    for (int k=0; k<u_c_.size(); ++k) {
      if (u_c_[k] > 0) lbx[k] = 0;
      if (u_c_[k] < 0) ubx[k] = 0;
    }

    // NLP parameters
    arg1[NLPSOL_P] = w;
    for (int i=0; i<n_in(); ++i) {
      if (i!=iin_) {
        int n = f_.nnz_in(i);
        if (arg[i]) {
          copy_n(arg[i], n, w);
        } else {
          fill_n(w, n, 0.);
        }
        w += n;
      }
    }

    // Primal solution
    double* x = w; w += n_;
    res1[NLPSOL_X] = x;

    // Solve the NLP
    solver_(arg1, res1, iw, w, 0);
    stats_["solver_stats"] = solver_.getStats();

    // Get the implicit variable
    if (res[iout_]) {
      copy_n(x, n_, res[iout_]);
    }

    // Check if any auxilary outputs to evaluate
    bool has_aux = false;
    for (int i=0; i<n_out(); ++i) {
      if (i!=iout_ && res[i]) {
        has_aux = true;
        break;
      }
    }

    // Evaluate auxilary outputs, if necessary
    if (has_aux) {
      copy_n(arg, n_in(), arg1);
      arg1[iin_] = x;
      copy_n(res, n_out(), res1);
      res1[iout_] = 0;
      f_(arg1, res1, iw, w, 0);
    }
  }

} // namespace casadi
