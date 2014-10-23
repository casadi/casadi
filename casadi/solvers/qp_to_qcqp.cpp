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


#include "qp_to_qcqp.hpp"

#include "casadi/core/sx/sx_tools.hpp"
#include "casadi/core/function/sx_function.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_QPSOLVER_QCQP_EXPORT
  casadi_register_qpsolver_qcqp(QpSolverInternal::Plugin* plugin) {
    plugin->creator = QpToQcqp::creator;
    plugin->name = "qcqp";
    plugin->doc = QpToQcqp::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_QPSOLVER_QCQP_EXPORT casadi_load_qpsolver_qcqp() {
    QpSolverInternal::registerPlugin(casadi_register_qpsolver_qcqp);
  }

  QpToQcqp* QpToQcqp::clone() const {
    // Return a deep copy
    QpToQcqp* node = new QpToQcqp(st_);
    if (!node->is_init_)
      node->init();
    return node;
  }

  QpToQcqp::QpToQcqp(const std::vector<Sparsity> &st) : QpSolverInternal(st) {
    Adaptor::addOptions();
  }

  QpToQcqp::~QpToQcqp() {
  }

  void QpToQcqp::evaluate() {

    // Pass inputs of QP to QCQP form
    solver_.input(QCQP_SOLVER_A).set(input(QP_SOLVER_A));
    solver_.input(QCQP_SOLVER_G).set(input(QP_SOLVER_G));
    solver_.input(QCQP_SOLVER_H).set(input(QP_SOLVER_H));

    solver_.input(QCQP_SOLVER_LBX).set(input(QP_SOLVER_LBX));
    solver_.input(QCQP_SOLVER_UBX).set(input(QP_SOLVER_UBX));

    solver_.input(QCQP_SOLVER_LBA).set(input(QP_SOLVER_LBA));
    solver_.input(QCQP_SOLVER_UBA).set(input(QP_SOLVER_UBA));

    // Delegate computation to QCQP Solver
    solver_.evaluate();

    // Pass the stats
    stats_["qcqp_solver_stats"] = solver_.getStats();

    // Read the outputs from Ipopt
    output(QCQP_SOLVER_X).set(solver_.output(QP_SOLVER_X));
    output(QCQP_SOLVER_COST).set(solver_.output(QP_SOLVER_COST));
    output(QCQP_SOLVER_LAM_A).set(solver_.output(QP_SOLVER_LAM_A));
    output(QCQP_SOLVER_LAM_X).set(solver_.output(QP_SOLVER_LAM_X));
  }

  void QpToQcqp::init() {
    // Initialize the base classes
    QpSolverInternal::init();

    // Create an QcqpSolver instance
    solver_ = QcqpSolver(getOption(solvername()),
                         qcqpStruct("h", input(QP_SOLVER_H).sparsity(),
                                    "p", Sparsity::sparse(n_, 0),
                                    "a", input(QP_SOLVER_A).sparsity()));
    solver_.setQPOptions();
    if (hasSetOption(optionsname())) solver_.setOption(getOption(optionsname()));
    solver_.init();
  }

} // namespace casadi
