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
  int CASADI_NLSOL_NLPSOL_EXPORT
  casadi_register_nlsol_nlpsol(Nlsol::Plugin* plugin) {
    plugin->creator = QpToImplicit::creator;
    plugin->name = "nlpsol";
    plugin->doc = QpToImplicit::meta_doc.c_str();
    plugin->version = 23;
    plugin->adaptorHasPlugin = Function::has_nlpsol;
    return 0;
  }

  extern "C"
  void CASADI_NLSOL_NLPSOL_EXPORT casadi_load_nlsol_nlpsol() {
    Nlsol::registerPlugin(casadi_register_nlsol_nlpsol);
  }

  QpToImplicit::QpToImplicit(const std::string& name, const Function& f)
    : Nlsol(name, f) {

    Adaptor<QpToImplicit, Nlpsol>::addOptions();
  }

  QpToImplicit::~QpToImplicit() {
  }

  void QpToImplicit::evalD(void* mem, const double** arg, double** res, int* iw, double* w) {
    // Copy to buffers
    for (int i=0; i<n_in(); ++i) {
      if (arg[i]) {
        setInputNZ(arg[i], i);
      } else {
        setInputNZ(0, i);
      }
    }
    setOutput(input(iin_), iout_);

    // Equality nonlinear constraints
    solver_.input(NLPSOL_LBG).set(0.);
    solver_.input(NLPSOL_UBG).set(0.);

    // Simple constraints
    vector<double>& lbx = solver_.input(NLPSOL_LBX).data();
    vector<double>& ubx = solver_.input(NLPSOL_UBX).data();
    for (int k=0; k<u_c_.size(); ++k) {
      lbx[k] = u_c_[k] <= 0 ? -std::numeric_limits<double>::infinity() : 0;
      ubx[k] = u_c_[k] >= 0 ?  std::numeric_limits<double>::infinity() : 0;
    }

    // Pass initial guess
    solver_.input(NLPSOL_X0).set(output(iout_));

    // Add auxiliary inputs
    auto nlp_p = solver_.input(NLPSOL_P)->begin();
    for (int i=0; i<n_in(); ++i) {
      if (i!=iin_) {
        std::copy(input(i)->begin(), input(i)->end(), nlp_p);
        nlp_p += input(i).nnz();
      }
    }

    // Solve the NLP
    solver_.evaluate();
    stats_["solver_stats"] = solver_.getStats();

    // Get the implicit variable
    output(iout_).set(solver_.output(NLPSOL_X));

    // Evaluate auxilary outputs, if necessary
    if (n_out()>0) {
      f_.setInput(output(iout_), iin_);
      for (int i=0; i<n_in(); ++i)
        if (i!=iin_) f_.setInput(input(i), i);
      f_.evaluate();
      for (int i=0; i<n_out(); ++i) {
        if (i!=iout_) f_.getOutput(output(i), i);
      }
    }

    // Get from buffers
    for (int i=0; i<n_out(); ++i) {
      if (res[i]) {
        getOutputNZ(res[i], i);
      }
    }
  }

  void QpToImplicit::init() {
    // Call the base class initializer
    Nlsol::init();

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
    if (hasSetOption(optionsname())) options = option(optionsname());
    // Create an Nlpsol instance
    solver_ = Function::nlpsol("nlpsol", option(solvername()), nlp, options);
  }

} // namespace casadi
