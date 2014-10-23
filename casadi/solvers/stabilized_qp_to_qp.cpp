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


#include "stabilized_qp_to_qp.hpp"

#include "casadi/core/sx/sx_tools.hpp"
#include "casadi/core/function/sx_function.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_STABILIZEDQPSOLVER_QP_EXPORT
  casadi_register_stabilizedqpsolver_qp(StabilizedQpSolverInternal::Plugin* plugin) {
    plugin->creator = StabilizedQpToQp::creator;
    plugin->name = "qp";
    plugin->doc = StabilizedQpToQp::meta_doc.c_str();
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_STABILIZEDQPSOLVER_QP_EXPORT casadi_load_stabilizedqpsolver_qp() {
    StabilizedQpSolverInternal::registerPlugin(casadi_register_stabilizedqpsolver_qp);
  }

  StabilizedQpToQp::StabilizedQpToQp(const std::vector<Sparsity> &st)
      : StabilizedQpSolverInternal(st) {
    addOption("qp_solver",         OT_STRING,   GenericType(),
              "The QP solver used to solve the stabilized QPs.");
    addOption("qp_solver_options", OT_DICTIONARY, GenericType(),
              "Options to be passed to the QP solver instance");
  }

  StabilizedQpToQp::~StabilizedQpToQp() {
  }

  void StabilizedQpToQp::deepCopyMembers(
      std::map<SharedObjectNode*, SharedObject>& already_copied) {
    StabilizedQpSolverInternal::deepCopyMembers(already_copied);
    qp_solver_ = deepcopy(qp_solver_, already_copied);
  }

  void StabilizedQpToQp::init() {
    // Initialize the base classes
    StabilizedQpSolverInternal::init();

    // Form augmented QP
    Sparsity H_sparsity_qp = blkdiag(st_[QP_STRUCT_H], Sparsity::diag(nc_));
    Sparsity A_sparsity_qp = horzcat(st_[QP_STRUCT_A], Sparsity::diag(nc_));
    std::string qp_solver_name = getOption("qp_solver");
    qp_solver_ = QpSolver(qp_solver_name,
                          qpStruct("h", H_sparsity_qp, "a", A_sparsity_qp));

    // Pass options if provided
    if (hasSetOption("qp_solver_options")) {
      Dictionary qp_solver_options = getOption("qp_solver_options");
      qp_solver_.setOption(qp_solver_options);
    }

    // Initialize the QP solver
    qp_solver_.init();
  }

  void StabilizedQpToQp::evaluate() {
    double muR = input(STABILIZED_QP_SOLVER_MUR).at(0);
    std::vector<double>& muE = input(STABILIZED_QP_SOLVER_MUE).data();
    std::vector<double>& mu = input(STABILIZED_QP_SOLVER_MU).data();

    // Construct stabilized H
    DMatrix& H_qp = qp_solver_.input(QP_SOLVER_H);
    std::copy(input(STABILIZED_QP_SOLVER_H).begin(),
              input(STABILIZED_QP_SOLVER_H).end(),
              H_qp.begin());
    std::fill(H_qp.begin()+input(STABILIZED_QP_SOLVER_H).size(), H_qp.end(), muR);

    // Linear constraints
    if (nc_>0) {
      // Upper and lower bounds
      qp_solver_.setInput(input(STABILIZED_QP_SOLVER_LBA), QP_SOLVER_LBA);
      qp_solver_.setInput(input(STABILIZED_QP_SOLVER_UBA), QP_SOLVER_UBA);

      // Matrix term in the linear constraint
      DMatrix& A_qp = qp_solver_.input(QP_SOLVER_A);
      DMatrix& A = input(STABILIZED_QP_SOLVER_A);
      std::copy(A.begin(), A.end(), A_qp.begin());
      std::fill(A_qp.begin()+A.size(), A_qp.end(), -muR);

      // Add constant to linear inequality
      for (int i=0; i<mu.size(); ++i) {
        double extra = muR*(mu[i]-muE[i]);
        qp_solver_.input(QP_SOLVER_LBA).at(i) += extra;
        qp_solver_.input(QP_SOLVER_UBA).at(i) += extra;
      }
    }

    // Bounds on x
    std::copy(input(STABILIZED_QP_SOLVER_LBX).begin(),
              input(STABILIZED_QP_SOLVER_LBX).end(),
              qp_solver_.input(QP_SOLVER_LBX).begin());
    std::copy(input(STABILIZED_QP_SOLVER_UBX).begin(),
              input(STABILIZED_QP_SOLVER_UBX).end(),
              qp_solver_.input(QP_SOLVER_UBX).begin());
    std::fill(qp_solver_.input(QP_SOLVER_LBX).begin()+n_,
              qp_solver_.input(QP_SOLVER_LBX).end(),
              -numeric_limits<double>::infinity());
    std::fill(qp_solver_.input(QP_SOLVER_UBX).begin()+n_,
              qp_solver_.input(QP_SOLVER_UBX).end(),
              numeric_limits<double>::infinity());

    // Gradient term in the objective
    DMatrix &g = input(STABILIZED_QP_SOLVER_G);
    std::copy(g.begin(), g.end(), qp_solver_.input(QP_SOLVER_G).begin());
    for (int i=0;i<nc_;++i) {
      qp_solver_.input(QP_SOLVER_G).at(g.size()+i) = muR * mu[i];
    }

    // Hot-starting if possible
    std::copy(input(STABILIZED_QP_SOLVER_X0).begin(),
              input(STABILIZED_QP_SOLVER_X0).end(),
              qp_solver_.input(QP_SOLVER_X0).begin());

    // Solve the QP
    qp_solver_.evaluate();

    // Pass the stats
    stats_["qp_solver_stats"] = qp_solver_.getStats();

    // Get the optimal solution
    std::copy(qp_solver_.output(QP_SOLVER_X).begin(),
              qp_solver_.output(QP_SOLVER_X).begin()+n_,
              output(QP_SOLVER_X).begin());
    std::copy(qp_solver_.output(QP_SOLVER_LAM_X).begin(),
              qp_solver_.output(QP_SOLVER_LAM_X).begin()+n_,
              output(QP_SOLVER_LAM_X).begin());
    std::copy(qp_solver_.output(QP_SOLVER_LAM_A).begin(),
              qp_solver_.output(QP_SOLVER_LAM_A).begin()+nc_,
              output(QP_SOLVER_LAM_A).begin());
  }

  void StabilizedQpToQp::generateNativeCode(std::ostream &file) const {
    qp_solver_.generateNativeCode(file);
  }

} // namespace casadi
