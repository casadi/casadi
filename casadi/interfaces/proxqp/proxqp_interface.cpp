/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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


#include "proxqp_interface.hpp"
#include "casadi/core/casadi_misc.hpp"


using namespace proxsuite::proxqp;

namespace casadi {

  extern "C"
  int CASADI_CONIC_PROXQP_EXPORT
  casadi_register_conic_proxqp(Conic::Plugin* plugin) {
    plugin->creator = ProxqpInterface::creator;
    plugin->name = "proxqp";
    plugin->doc = ProxqpInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &ProxqpInterface::options_;
    plugin->deserialize = &ProxqpInterface::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_PROXQP_EXPORT casadi_load_conic_proxqp() {
    Conic::registerPlugin(casadi_register_conic_proxqp);
  }

  ProxqpInterface::ProxqpInterface(const std::string& name,
                                   const std::map<std::string, Sparsity>& st)
    : Conic(name, st), sparse_backend(true), max_iter(0.0) {
  }

  ProxqpInterface::~ProxqpInterface() {
    clear_mem();
  }

  const Options ProxqpInterface::options_
  = {{&Conic::options_},
     {{"proxqp",
       {OT_DICT,
        "const proxqp options."}},
      {"warm_start_primal",
       {OT_BOOL,
        "Use x input to warmstart [Default: true]."}},
      {"warm_start_dual",
       {OT_BOOL,
        "Use y and z input to warmstart [Default: true]."}}
     }
  };

  void ProxqpInterface::init(const Dict& opts) {
    // Initialize the base classes
    Conic::init(opts);

    warm_start_primal_ = true;
    warm_start_dual_ = true;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="warm_start_primal") {
        warm_start_primal_ = op.second;
      } else if (op.first=="warm_start_dual") {
        warm_start_dual_ = op.second;
      } else if (op.first=="proxqp") {
        const Dict& opts = op.second;
        for (auto&& op : opts) {
          if (op.first=="default_rho") {
            settings_.default_rho = op.second;
          } else if (op.first=="default_mu_eq") {
            settings_.default_mu_eq = op.second;
          } else if (op.first=="default_mu_in") {
            settings_.default_mu_in = op.second;
          } else if (op.first=="eps_abs") {
            settings_.eps_abs = op.second;
          } else if (op.first=="eps_rel") {
            settings_.eps_rel = op.second;
          } else if (op.first=="max_iter") {
            settings_.max_iter = static_cast<double>(op.second);
            max_iter = static_cast<double>(op.second);
          } else if (op.first=="verbose") {
            settings_.verbose = op.second;
          } else if (op.first=="backend") {
            if (op.second == "sparse") {
              sparse_backend = true;
            } else if (op.second == "dense") {
              sparse_backend = false;
            } else {
              casadi_error("[Backend option] Please specify either sparse or dense");
            }
          } else {
            casadi_error("[ProxQP settings] User-specified option "
                         "'" + str(op.first) + "' not recognized.");
          }
        }
      }
    }

    // Allocate memory for problem
    nA_ = nnz_in(CONIC_A);
    nH_ = nnz_in(CONIC_H);

    // Allocate work vectors
    alloc_w(nx_, true); // g
    alloc_w(nx_, true); // lbx
    alloc_w(nx_, true); // ubx
    alloc_w(na_, true); // lba
    alloc_w(na_, true); // uba
    alloc_w(nH_, true); // H
    alloc_w(nA_, true); // A

  }

  int ProxqpInterface::init_mem(void* mem) const {
    if (Conic::init_mem(mem)) return 1;
    auto m = static_cast<ProxqpMemory*>(mem);

    m->tripletList.reserve(H_.nnz());
    m->tripletListEq.reserve(na_);

    m->g_vector.resize(nx_);
    m->uba_vector.resize(na_);
    m->lba_vector.resize(na_);
    m->ubx_vector.resize(na_);
    m->lbx_vector.resize(na_);
    m->ub_vector.resize(na_);
    m->lb_vector.resize(na_);
    m->b_vector.resize(na_);

    m->add_stat("preprocessing");
    m->add_stat("solver");
    m->add_stat("postprocessing");
    return 0;
  }

  inline const char* return_status_string(casadi_int status) {
    return "Unknown";
  }

  int ProxqpInterface::
  solve(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    typedef Eigen::Triplet<double> T;

    auto m = static_cast<ProxqpMemory*>(mem);
    m->fstats.at("preprocessing").tic();

    // Get problem data
    double* g=w; w += nx_;
    casadi_copy(arg[CONIC_G], nx_, g);
    double* lbx=w; w += nx_;
    casadi_copy(arg[CONIC_LBX], nx_, lbx);
    double* ubx=w; w += nx_;
    casadi_copy(arg[CONIC_UBX], nx_, ubx);
    double* lba=w; w += na_;
    casadi_copy(arg[CONIC_LBA], na_, lba);
    double* uba=w; w += na_;
    casadi_copy(arg[CONIC_UBA], na_, uba);
    double* H=w; w += nnz_in(CONIC_H);
    casadi_copy(arg[CONIC_H], nnz_in(CONIC_H), H);
    double* A=w; w += nnz_in(CONIC_A);
    casadi_copy(arg[CONIC_A], nnz_in(CONIC_A), A);

    m->g_vector = Eigen::Map<const Eigen::VectorXd>(g, nx_);
    m->uba_vector = Eigen::Map<Eigen::VectorXd>(uba, na_);
    m->lba_vector = Eigen::Map<Eigen::VectorXd>(lba, na_);
    m->ubx_vector = Eigen::Map<Eigen::VectorXd>(ubx, nx_);
    m->lbx_vector = Eigen::Map<Eigen::VectorXd>(lbx, nx_);

    // Use lhs_equals_rhs_constraint to split double-sided bounds into one-sided
    // bound for equality constraints and double-sided for inequality constraints
    const Eigen::Array<bool, Eigen::Dynamic, 1>
    lhs_equals_rhs_constraint = (m->uba_vector.array() == m->lba_vector.array()).eval();
    std::vector<unsigned int> number_of_prev_equality(lhs_equals_rhs_constraint.size(), 0);
    std::vector<unsigned int> number_of_prev_inequality(lhs_equals_rhs_constraint.size(), 0);
    std::vector<double> tmp_eq_vector;
    std::vector<double> tmp_ineq_lb_vector;
    std::vector<double> tmp_ineq_ub_vector;
    {

      // number_of_prev_equality and number_of_prev
      // inequality are two vectors that contains the number of
      // equality and inequality that can be found before the current index
      // number_of_prev_equality[i] = number of equality that can be found before index i
      // number_of_prev_inequality[i] = number of inequality that can be found before index i
      // For instance:
      //     equality and inequality   [i, e, e, e, i, i, i, e]
      //     lhs_equals_rgs_contraint  [f, t, t, t, f, f, f, t]
      //     number_of_prev_equality   [0, 0, 1, 3, 3, 3, 3, 3]
      //     number_of_prev_inequality [0, 1, 1, 1, 1, 2, 3, 4]
      for (std::size_t k=1; k<lhs_equals_rhs_constraint.size(); ++k) {
        if (lhs_equals_rhs_constraint[k-1]) {
          number_of_prev_equality[k] = number_of_prev_equality[k-1] + 1;
          number_of_prev_inequality[k] = number_of_prev_inequality[k-1];
        } else {
          number_of_prev_inequality[k] = number_of_prev_inequality[k-1] + 1;
          number_of_prev_equality[k] = number_of_prev_equality[k-1];
        }
      }

      for (std::size_t k=0; k<lhs_equals_rhs_constraint.size(); ++k) {
        if (lhs_equals_rhs_constraint[k]) {
          tmp_eq_vector.push_back(m->lba_vector[k]);
        } else {
          tmp_ineq_lb_vector.push_back(m->lba_vector[k]);
          tmp_ineq_ub_vector.push_back(m->uba_vector[k]);
        }
      }

      m->b_vector.resize(tmp_eq_vector.size());
      if (tmp_eq_vector.size() > 0) {
        m->b_vector = Eigen::Map<Eigen::VectorXd>(
          get_ptr(tmp_eq_vector), tmp_eq_vector.size());
      }

      m->lba_vector.resize(tmp_ineq_lb_vector.size());
      if (tmp_ineq_lb_vector.size() > 0) {
        m->lba_vector = Eigen::Map<Eigen::VectorXd>(
          get_ptr(tmp_ineq_lb_vector), tmp_ineq_lb_vector.size());
      }

      m->uba_vector.resize(tmp_ineq_ub_vector.size());
      if (tmp_ineq_ub_vector.size() > 0) {
        m->uba_vector = Eigen::Map<Eigen::VectorXd>(
          get_ptr(tmp_ineq_ub_vector), tmp_ineq_ub_vector.size());
      }
    }
    std::size_t n_eq = m->b_vector.size();
    std::size_t n_ineq = m->lba_vector.size() + m->lbx_vector.size();

    // Convert H_ from casadi::Sparsity to Eigen::SparseMatrix
    H_.get_triplet(m->row, m->col);
    for (int k=0; k<H_.nnz(); ++k) {
      m->tripletList.push_back(T(
        static_cast<double>(m->row[k]),
        static_cast<double>(m->col[k]),
        static_cast<double>(H[k])));
    }
    Eigen::SparseMatrix<double> H_spa(H_.size1(), H_.size2());
    H_spa.setFromTriplets(m->tripletList.begin(), m->tripletList.end());
    m->tripletList.clear();

    // Convert A_ from casadi Sparsity to Eigen::SparseMatrix and split
    // in- and equality constraints into different matrices
    m->tripletList.reserve(A_.nnz());
    A_.get_triplet(m->row, m->col);

    for (int k=0; k<A_.nnz(); ++k) {
      // Detect equality constraint
      if (lhs_equals_rhs_constraint[m->row[k]]) {
        // Equality constraint the row[k] is decreased
        // by the number of previous inequality constraints
        m->tripletListEq.push_back(T(
          static_cast<double>(m->row[k] - number_of_prev_inequality[m->row[k]]),
          static_cast<double>(m->col[k]),
          static_cast<double>(A[k])));
      } else {
        // Inequality constraint the row[k] is decreased
        // by the number of previous equality constraints
        m->tripletList.push_back(T(
          static_cast<double>(m->row[k] - number_of_prev_equality[m->row[k]]),
          static_cast<double>(m->col[k]),
          static_cast<double>(A[k])));
      }
    }

    // Handle constraints on decision variable x in inequality constraint matrix C
    uint32_t n_constraints_x = 0;
    if (m->ubx_vector.size() > 0 || m->lbx_vector.size() > 0) {
      for (int k=0; k<nx_; ++k) {
        m->tripletList.push_back(T(
          static_cast<double>(m->lba_vector.size() + k),
          static_cast<double>(k),
          static_cast<double>(1.0)));
        ++n_constraints_x;
      }
    }

    Eigen::SparseMatrix<double> A_spa(n_eq, nx_);
    A_spa.setFromTriplets(m->tripletListEq.begin(), m->tripletListEq.end());
    m->tripletListEq.clear();

    Eigen::SparseMatrix<double> C_spa(n_ineq, nx_);
    C_spa.setFromTriplets(m->tripletList.begin(), m->tripletList.end());
    m->tripletList.clear();

    // Get stacked lower and upper inequality bounds
    m->ub_vector.resize(n_ineq);
    m->lb_vector.resize(n_ineq);
    if (m->uba_vector.size() > 0) {
      m->ub_vector << m->uba_vector, m->ubx_vector;
    } else {
      m->ub_vector << m->ubx_vector;
    }

    if (m->lba_vector.size() > 0) {
      m->lb_vector << m->lba_vector, m->lbx_vector;
    } else {
      m->lb_vector << m->lbx_vector;
    }
    m->fstats.at("preprocessing").toc();

    // Solve Problem
    m->fstats.at("solver").tic();

    if (sparse_backend) {
      m->sparse_solver = proxsuite::proxqp::sparse::QP<double, long long> (nx_, n_eq, n_ineq);
      m->sparse_solver.init(H_spa, m->g_vector,
                            A_spa, m->b_vector,
                            C_spa, m->lb_vector, m->ub_vector);
      m->sparse_solver.settings = settings_;

      m->sparse_solver.solve();

      m->results_x = std::make_unique<Eigen::VectorXd>(m->sparse_solver.results.x);
      m->results_y = std::make_unique<Eigen::VectorXd>(m->sparse_solver.results.y);
      m->results_z = std::make_unique<Eigen::VectorXd>(m->sparse_solver.results.z);
      m->objValue = m->sparse_solver.results.info.objValue;
      m->status = m->sparse_solver.results.info.status;
    } else {
      m->dense_solver = proxsuite::proxqp::dense::QP<double> (nx_, n_eq, n_ineq);
      m->dense_solver.init(Eigen::MatrixXd(H_spa), m->g_vector,
        Eigen::MatrixXd(A_spa), m->b_vector,
        Eigen::MatrixXd(C_spa), m->lb_vector, m->ub_vector);
      m->dense_solver.settings = settings_;

      m->dense_solver.solve();

      m->results_x = std::make_unique<Eigen::VectorXd>(m->dense_solver.results.x);
      m->results_y = std::make_unique<Eigen::VectorXd>(m->dense_solver.results.y);
      m->results_z = std::make_unique<Eigen::VectorXd>(m->dense_solver.results.z);
      m->objValue = m->dense_solver.results.info.objValue;
      m->status = m->dense_solver.results.info.status;
    }
    m->fstats.at("solver").toc();

    // Post-processing to retrieve the results
    m->fstats.at("postprocessing").tic();
    casadi_copy(m->results_x->data(), nx_, res[CONIC_X]);

    // Copy back the multipliers.
    // Note, casadi has LAM_X (multipliers for constraints on variable x) and
    // LAM_A (multipliers for in- and equality constraints) while proxqp has results_y
    // (equality multipliers) and results_z (inequality multipliers).
    if (n_constraints_x > 0) {
      casadi_copy(m->results_z->tail(n_constraints_x).data(), n_constraints_x, res[CONIC_LAM_X]);
    }
    if (!n_eq) {
      uint32_t n_lam_a = m->results_z->size() - n_constraints_x;
      casadi_assert_dev(n_lam_a == na_);
      casadi_copy(m->results_z->head(n_lam_a).data(), n_lam_a, res[CONIC_LAM_A]);
    } else if (!n_ineq) {
      casadi_copy(m->results_z->data(), na_, res[CONIC_LAM_A]);
    } else {
      bool ineq_constraints_a = (lhs_equals_rhs_constraint == 0).any();
      if (!ineq_constraints_a) {
        casadi_copy(m->results_y->data(), na_, res[CONIC_LAM_A]);
      } else {
        Eigen::VectorXd lam_a(na_);

        for (int k=0; k<lhs_equals_rhs_constraint.size(); ++k) {
          if (lhs_equals_rhs_constraint[k]) {
            lam_a[k] = m->results_y->coeff(k);
          } else {
            lam_a[k] = m->results_z->coeff(k);
          }
        }
        casadi_copy(lam_a.data(), na_, res[CONIC_LAM_A]);
      }
    }

    if (res[CONIC_COST]) {
      *res[CONIC_COST] = m->objValue;
    }

    m->d_qp.success = m->status == QPSolverOutput::PROXQP_SOLVED;
    if (m->d_qp.success) {
      m->d_qp.unified_return_status = SOLVER_RET_SUCCESS;
    } else {
      if (m->status == QPSolverOutput::PROXQP_MAX_ITER_REACHED) {
        m->d_qp.unified_return_status = SOLVER_RET_LIMITED;
      } else { // primal or dual infeasibility
        m->d_qp.unified_return_status = SOLVER_RET_UNKNOWN;
      }
    }
    m->fstats.at("postprocessing").toc();

    return 0;
  }


  Dict ProxqpInterface::get_stats(void* mem) const {

    Dict stats = Conic::get_stats(mem);
    auto m = static_cast<ProxqpMemory*>(mem);

    std::string ret_status;
    if (m->status == QPSolverOutput::PROXQP_SOLVED) {
      ret_status = "PROXQP_SOLVED";
    } else if (m->status == QPSolverOutput::PROXQP_MAX_ITER_REACHED) {
      ret_status = "PROXQP_MAX_ITER_REACHED";
    } else if (m->status == QPSolverOutput::PROXQP_PRIMAL_INFEASIBLE) {
      ret_status = "PROXQP_PRIMAL_INFEASIBLE";
    } else if (m->status == QPSolverOutput::PROXQP_DUAL_INFEASIBLE) {
      ret_status = "PROXQP_DUAL_INFEASIBLE";
    }

    stats["return_status"] = ret_status;
    return stats;
  }

  ProxqpMemory::ProxqpMemory()
    : sparse_solver(1, 0, 0),
      dense_solver(1, 0, 0) {
  }

  ProxqpMemory::~ProxqpMemory() {
    this->sparse_solver.cleanup();
    this->dense_solver.cleanup();
  }

  ProxqpInterface::ProxqpInterface(DeserializingStream& s) : Conic(s) {
    s.version("ProxqpInterface", 1);
    s.unpack("ProxqpInterface::warm_start_primal", warm_start_primal_);
    s.unpack("ProxqpInterface::warm_start_dual", warm_start_dual_);

    s.unpack("ProxqpInterface::settings::default_rho", settings_.default_rho);
    s.unpack("ProxqpInterface::settings::default_mu_eq", settings_.default_mu_eq);
    s.unpack("ProxqpInterface::settings::default_mu_in", settings_.default_mu_in);
    s.unpack("ProxqpInterface::settings::eps_abs", settings_.eps_abs);
    s.unpack("ProxqpInterface::settings::eps_rel", settings_.eps_rel);
    s.unpack("ProxqpInterface::settings::max_iter", max_iter);
    settings_.max_iter = isize(max_iter);
    s.unpack("ProxqpInterface::settings::verbose", settings_.verbose);
    s.unpack("ProxqpInterface::settings::sparse_backend", sparse_backend);
  }

  void ProxqpInterface::serialize_body(SerializingStream &s) const {
    Conic::serialize_body(s);
    s.version("ProxqpInterface", 1);
    s.pack("ProxqpInterface::warm_start_primal", warm_start_primal_);
    s.pack("ProxqpInterface::warm_start_dual", warm_start_dual_);
    s.pack("ProxqpInterface::settings::default_rho", settings_.default_rho);
    s.pack("ProxqpInterface::settings::default_mu_eq", settings_.default_mu_eq);
    s.pack("ProxqpInterface::settings::default_mu_in", settings_.default_mu_in);
    s.pack("ProxqpInterface::settings::eps_abs", settings_.eps_abs);
    s.pack("ProxqpInterface::settings::eps_rel", settings_.eps_rel);
    s.pack("ProxqpInterface::settings::max_iter", static_cast<double>(settings_.max_iter));
    s.pack("ProxqpInterface::settings::verbose", settings_.verbose);
    s.pack("ProxqpInterface::settings::sparse_backend", sparse_backend);
  }

} // namespace casadi
