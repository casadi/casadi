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


#include "lp_to_qp.hpp"

#include "casadi/core/sx/sx_tools.hpp"
#include "casadi/core/function/sx_function.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_LPSOLVER_QP_EXPORT
  casadi_register_lpsolver_qp(LpSolverInternal::Plugin* plugin) {
    plugin->creator = LpToQp::creator;
    plugin->name = "qp";
    plugin->doc = LpToQp::meta_doc.c_str();;
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_LPSOLVER_QP_EXPORT casadi_load_lpsolver_qp() {
    LpSolverInternal::registerPlugin(casadi_register_lpsolver_qp);
  }

  LpToQp* LpToQp::clone() const {
    // Return a deep copy
    LpToQp* node = new LpToQp(st_);
    if (!node->is_init_)
      node->init();
    return node;
  }

  LpToQp::LpToQp(const std::vector<Sparsity> &st) : LpSolverInternal(st) {
    Adaptor::addOptions();
  }

  LpToQp::~LpToQp() {
  }

  void LpToQp::evaluate() {

    // Pass inputs of LP to QP form
    solver_.input(QP_SOLVER_A).set(input(LP_SOLVER_A));
    solver_.input(QP_SOLVER_G).set(input(LP_SOLVER_C));

    solver_.input(QP_SOLVER_LBX).set(input(LP_SOLVER_LBX));
    solver_.input(QP_SOLVER_UBX).set(input(LP_SOLVER_UBX));

    solver_.input(QP_SOLVER_LBA).set(input(LP_SOLVER_LBA));
    solver_.input(QP_SOLVER_UBA).set(input(LP_SOLVER_UBA));

    // Delegate computation to NLP Solver
    solver_.evaluate();

    // Pass the stats
    stats_["qp_solver_stats"] = solver_.getStats();

    // Read the outputs from Ipopt
    output(QP_SOLVER_X).set(solver_.output(LP_SOLVER_X));
    output(QP_SOLVER_COST).set(solver_.output(LP_SOLVER_COST));
    output(QP_SOLVER_LAM_A).set(solver_.output(LP_SOLVER_LAM_A));
    output(QP_SOLVER_LAM_X).set(solver_.output(LP_SOLVER_LAM_X));
  }

  void LpToQp::init() {
    // Initialize the base classes
    LpSolverInternal::init();

    // Create a QpSolver instance
    solver_ = QpSolver(getOption(solvername()),
                       qpStruct("h", Sparsity::sparse(n_, n_),
                                "a", input(LP_SOLVER_A).sparsity()));
    if (hasSetOption(optionsname())) solver_.setOption(getOption(optionsname()));
    solver_.init();
  }

} // namespace casadi

