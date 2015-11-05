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
    plugin->adaptorHasPlugin = Function::has_nlpsol;
    return 0;
  }

  extern "C"
  void CASADI_QPSOL_NLPSOL_EXPORT casadi_load_qpsol_nlpsol() {
    Qpsol::registerPlugin(casadi_register_qpsol_nlpsol);
  }

  QpToNlp::QpToNlp(const std::string& name, const std::map<std::string, Sparsity> &st)
    : Qpsol(name, st) {

    Adaptor<QpToNlp, Nlpsol>::addOptions();
  }

  QpToNlp::~QpToNlp() {
  }

  void QpToNlp::evaluate() {
    if (inputs_check_) checkInputs();

    int k = 0;

    // Pass inputs of QP to NLP form

    std::copy(input(QPSOL_H).data().begin(),
              input(QPSOL_H).data().end(),
              solver_.input(NLPSOL_P).data().begin()+k); k+= input(QPSOL_H).nnz();
    std::copy(input(QPSOL_G).data().begin(),
              input(QPSOL_G).data().end(),
              solver_.input(NLPSOL_P).data().begin()+k); k+= input(QPSOL_G).nnz();
    std::copy(input(QPSOL_A).data().begin(),
              input(QPSOL_A).data().end(),
              solver_.input(NLPSOL_P).data().begin()+k);


    solver_.input(NLPSOL_LBX).set(input(QPSOL_LBX));
    solver_.input(NLPSOL_UBX).set(input(QPSOL_UBX));

    solver_.input(NLPSOL_LBG).set(input(QPSOL_LBA));
    solver_.input(NLPSOL_UBG).set(input(QPSOL_UBA));

    // Delegate computation to NLP Solver
    solver_.evaluate();

    // Pass the stats
    stats_["nlpsol_stats"] = solver_.getStats();

    // Read the outputs from Ipopt
    output(QPSOL_X).set(solver_.output(NLPSOL_X));
    output(QPSOL_COST).set(solver_.output(NLPSOL_F));
    output(QPSOL_LAM_A).set(solver_.output(NLPSOL_LAM_G));
    output(QPSOL_LAM_X).set(solver_.output(NLPSOL_LAM_X));
  }

  void QpToNlp::init() {
    // Initialize the base classes
    Qpsol::init();

    // Create a symbolic matrix for the decision variables
    SX X = SX::sym("X", n_, 1);

    // Parameters to the problem
    SX H = SX::sym("H", input(QPSOL_H).sparsity());
    SX G = SX::sym("G", input(QPSOL_G).sparsity());
    SX A = SX::sym("A", input(QPSOL_A).sparsity());

    // Put parameters in a vector
    std::vector< SX > par;
    par.push_back(H.data());
    par.push_back(G.data());
    par.push_back(A.data());



    // The nlp looks exactly like a mathematical description of the NLP
    SXDict nlp = {{"x", X}, {"p", vertcat(par)},
                  {"f", mul(G.T(), X) + 0.5*mul(mul(X.T(), H), X)}, {"g", mul(A, X)}};

    Dict options;
    if (hasSetOption(optionsname())) options = option(optionsname());
    options = OptionsFunctionality::addOptionRecipe(options, "qp");

    // Create an Nlpsol instance
    solver_ = Function::nlpsol("nlpsol", option(solvername()), nlp, options);
  }

} // namespace casadi
