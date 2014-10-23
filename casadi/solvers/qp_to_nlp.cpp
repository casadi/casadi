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

#include "casadi/core/sx/sx_tools.hpp"
#include "casadi/core/function/sx_function.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_QPSOLVER_NLP_EXPORT
  casadi_register_qpsolver_nlp(QpSolverInternal::Plugin* plugin) {
    plugin->creator = QpToNlp::creator;
    plugin->name = "nlp";
    plugin->doc = QpToNlp::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_QPSOLVER_NLP_EXPORT casadi_load_qpsolver_nlp() {
    QpSolverInternal::registerPlugin(casadi_register_qpsolver_nlp);
  }

  QpToNlp* QpToNlp::clone() const {
    // Return a deep copy
    QpToNlp* node = new QpToNlp(st_);
    if (!node->is_init_)
      node->init();
    return node;
  }

  QpToNlp::QpToNlp(const std::vector<Sparsity> &st) : QpSolverInternal(st) {
    Adaptor::addOptions();
  }

  QpToNlp::~QpToNlp() {
  }

  void QpToNlp::evaluate() {
    if (inputs_check_) checkInputs();

    int k = 0;

    // Pass inputs of QP to NLP form

    std::copy(input(QP_SOLVER_H).data().begin(),
              input(QP_SOLVER_H).data().end(),
              solver_.input(NLP_SOLVER_P).data().begin()+k); k+= input(QP_SOLVER_H).size();
    std::copy(input(QP_SOLVER_G).data().begin(),
              input(QP_SOLVER_G).data().end(),
              solver_.input(NLP_SOLVER_P).data().begin()+k); k+= input(QP_SOLVER_G).size();
    std::copy(input(QP_SOLVER_A).data().begin(),
              input(QP_SOLVER_A).data().end(),
              solver_.input(NLP_SOLVER_P).data().begin()+k);


    solver_.input(NLP_SOLVER_LBX).set(input(QP_SOLVER_LBX));
    solver_.input(NLP_SOLVER_UBX).set(input(QP_SOLVER_UBX));

    solver_.input(NLP_SOLVER_LBG).set(input(QP_SOLVER_LBA));
    solver_.input(NLP_SOLVER_UBG).set(input(QP_SOLVER_UBA));

    // Delegate computation to NLP Solver
    solver_.evaluate();

    // Pass the stats
    stats_["nlp_solver_stats"] = solver_.getStats();

    // Read the outputs from Ipopt
    output(QP_SOLVER_X).set(solver_.output(NLP_SOLVER_X));
    output(QP_SOLVER_COST).set(solver_.output(NLP_SOLVER_F));
    output(QP_SOLVER_LAM_A).set(solver_.output(NLP_SOLVER_LAM_G));
    output(QP_SOLVER_LAM_X).set(solver_.output(NLP_SOLVER_LAM_X));
  }

  void QpToNlp::init() {
    // Initialize the base classes
    QpSolverInternal::init();

    // Create a symbolic matrix for the decision variables
    SX X = SX::sym("X", n_, 1);

    // Parameters to the problem
    SX H = SX::sym("H", input(QP_SOLVER_H).sparsity());
    SX G = SX::sym("G", input(QP_SOLVER_G).sparsity());
    SX A = SX::sym("A", input(QP_SOLVER_A).sparsity());

    // Put parameters in a vector
    std::vector< SX > par;
    par.push_back(H.data());
    par.push_back(G.data());
    par.push_back(A.data());

    // The nlp looks exactly like a mathematical description of the NLP
    SXFunction QP_SOLVER_nlp(nlpIn("x", X, "p", vertcat(par)),
                             nlpOut("f", mul(G.T(), X) + 0.5*mul(mul(X.T(), H), X),
                                    "g", mul(A, X)));

    // Create an NlpSolver instance
    solver_ = NlpSolver(getOption(solvername()), QP_SOLVER_nlp);
    solver_.setQPOptions();
    if (hasSetOption(optionsname())) solver_.setOption(getOption(optionsname()));
    solver_.init();
  }

} // namespace casadi
