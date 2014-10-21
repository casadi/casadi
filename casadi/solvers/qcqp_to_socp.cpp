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


#include "qcqp_to_socp.hpp"

#include "casadi/core/sx/sx_tools.hpp"
#include "casadi/core/function/sx_function.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_QCQPSOLVER_SOCP_EXPORT
  casadi_register_qcqpsolver_socp(QcqpSolverInternal::Plugin* plugin) {
    plugin->creator = QcqpToSocp::creator;
    plugin->name = "socp";
    plugin->doc = QcqpToSocp::meta_doc.c_str();;
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_QCQPSOLVER_SOCP_EXPORT casadi_load_qcqpsolver_socp() {
    QcqpSolverInternal::registerPlugin(casadi_register_qcqpsolver_socp);
  }

  QcqpToSocp* QcqpToSocp::clone() const {
    // Return a deep copy
    QcqpToSocp* node = new QcqpToSocp(st_);
    if (!node->is_init_)
      node->init();
    return node;
  }

  QcqpToSocp::QcqpToSocp(const std::vector<Sparsity> &st) : QcqpSolverInternal(st) {
    Adaptor::addOptions();
  }

  QcqpToSocp::~QcqpToSocp() {
  }

  void QcqpToSocp::evaluate() {
    if (inputs_check_) checkInputs();

    // Pass inputs of QCQP to SOCP form
    std::copy(input(QCQP_SOLVER_A).begin(),
              input(QCQP_SOLVER_A).end(),
              solver_.input(SOCP_SOLVER_A).begin());

    // Transform QCQP_SOLVER_P to SOCP_SOLVER_G
    // G = chol(H/2)
    int qcqp_p_offset = 0;
    cholesky_[0].input(0).setArray(&input(QCQP_SOLVER_H).data().front()+qcqp_p_offset);
    for (int i=0;i<nq_;++i) {
      cholesky_[i+1].input(0).setArray(&input(QCQP_SOLVER_P).data().front()+qcqp_p_offset);
      qcqp_p_offset+= cholesky_[i+1].input(0).size();
    }

    for (int i=0;i<nq_+1;++i) {
      for (int k=0;k<cholesky_[i].input(0).size();++k) {
        cholesky_[i].input(0).at(k)*=0.5;
      }
    }

    int socp_g_offset = 0;
    for (int i=0;i<nq_+1;++i) {
      cholesky_[i].prepare();
      DMatrix G = cholesky_[i].getFactorization(true);
      std::copy(G.begin(), G.end(), solver_.input(SOCP_SOLVER_G).begin()+socp_g_offset);
      socp_g_offset += G.size()+1;
    }

    // Transform QCQP_SOLVER_Q to SOCP_SOLVER_H (needs SOCP_SOLVER_G)
    //   2h'G = Q   -> 2G'h = Q  ->  solve for h
    double *x = &solver_.input(SOCP_SOLVER_H).data().front();
    std::copy(input(QCQP_SOLVER_G).begin(), input(QCQP_SOLVER_G).begin()+n_, x);
    cholesky_[0].solveL(x, 1, true);
    int x_offset = n_+1;
    for (int i=0;i<nq_;++i) {
      std::copy(input(QCQP_SOLVER_Q).begin()+i*(n_+1),
                input(QCQP_SOLVER_Q).begin()+(i+1)*(n_+1)-1,
                x+x_offset);
      cholesky_[i+1].solveL(x+x_offset, 1, true);
      x_offset += n_+1;
    }

    for (int k=0;k<solver_.input(SOCP_SOLVER_H).size();++k) {
      solver_.input(SOCP_SOLVER_H).at(k)*= 0.5;
    }

    // Transform QCQP_SOLVER_R to SOCP_SOLVER_F (needs SOCP_SOLVER_H)
    // f = sqrt(h'h-r)

    for (int i=0;i<nq_;++i) {
      solver_.input(SOCP_SOLVER_F)[i+1] = -input(QCQP_SOLVER_R)[i];
    }
    solver_.input(SOCP_SOLVER_F)[0] = 0;

    for (int i=0;i<nq_+1;++i) {
      for (int k=0;k<n_+1;++k) {
        double v = solver_.input(SOCP_SOLVER_H).at(i*(n_+1)+k);
        solver_.input(SOCP_SOLVER_F)[i]+=v*v;
      }
      solver_.input(SOCP_SOLVER_F).at(i) = sqrt(solver_.input(SOCP_SOLVER_F).at(i));
    }

    solver_.input(SOCP_SOLVER_E)[n_] = 0.5/solver_.input(SOCP_SOLVER_F)[0];

    // Fix the first qc arising from epigraph reformulation: we must make use of e here.
    solver_.input(SOCP_SOLVER_G)[cholesky_[0].getFactorization().size()] =
      solver_.input(SOCP_SOLVER_E)[n_];

    /// Objective of the epigraph form
    solver_.input(SOCP_SOLVER_C)[n_] = 1;

    std::copy(input(QCQP_SOLVER_LBX).begin(),
              input(QCQP_SOLVER_LBX).end(),
              solver_.input(SOCP_SOLVER_LBX).begin());
    std::copy(input(QCQP_SOLVER_UBX).begin(),
              input(QCQP_SOLVER_UBX).end(),
              solver_.input(SOCP_SOLVER_UBX).begin());

    std::copy(input(QCQP_SOLVER_LBA).begin(),
              input(QCQP_SOLVER_LBA).end(),
              solver_.input(SOCP_SOLVER_LBA).begin());
    std::copy(input(QCQP_SOLVER_UBA).begin(),
              input(QCQP_SOLVER_UBA).end(),
              solver_.input(SOCP_SOLVER_UBA).begin());

    // Delegate computation to SOCP Solver
    solver_.evaluate();

    // Pass the stats
    stats_["socp_solver_stats"] = solver_.getStats();

    // Read the outputs from SOCP Solver
    output(SOCP_SOLVER_COST).set(solver_.output(QCQP_SOLVER_COST));
    output(SOCP_SOLVER_LAM_A).set(solver_.output(QCQP_SOLVER_LAM_A));
    std::copy(solver_.output(QCQP_SOLVER_X).begin(),
              solver_.output(QCQP_SOLVER_X).begin()+n_,
              output(SOCP_SOLVER_X).begin());
    std::copy(solver_.output(QCQP_SOLVER_LAM_X).begin(),
              solver_.output(QCQP_SOLVER_LAM_X).begin()+n_,
              output(SOCP_SOLVER_LAM_X).begin());
  }

  void QcqpToSocp::init() {
    // Initialize the base classes
    QcqpSolverInternal::init();

    // Collection of sparsities that will make up SOCP_SOLVER_G
    std::vector<Sparsity> socp_g;

    // Allocate Cholesky solvers
    cholesky_.push_back(LinearSolver("csparsecholesky", st_[QCQP_STRUCT_H]));
    for (int i=0;i<nq_;++i) {
      cholesky_.push_back(
        LinearSolver("csparsecholesky",
          DMatrix(st_[QCQP_STRUCT_P])(range(i*n_, (i+1)*n_), ALL).sparsity()));
    }

    for (int i=0;i<nq_+1;++i) {
      // Initialize Cholesky solve
      cholesky_[i].init();

      // Harvest Cholsesky sparsity patterns
      // Note that we add extra scalar to make room for the epigraph-reformulation variable
      socp_g.push_back(blkdiag(
        cholesky_[i].getFactorizationSparsity(false), Sparsity::dense(1, 1)));
    }

    // Create an SocpSolver instance
    solver_ = SocpSolver(getOption(solvername()),
                         socpStruct("g", horzcat(socp_g),
                                    "a", horzcat(input(QCQP_SOLVER_A).sparsity(),
                                                 Sparsity::sparse(nc_, 1))));
    //solver_.setQCQPOptions();
    if (hasSetOption(optionsname())) solver_.setOption(getOption(optionsname()));

    std::vector<int> ni(nq_+1);
    for (int i=0;i<nq_+1;++i) {
      ni[i] = n_+1;
    }
    solver_.setOption("ni", ni);

    // Initialize the SocpSolver
    solver_.init();
  }

} // namespace casadi
