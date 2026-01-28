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


#include "piqp_interface.hpp"
#include "casadi/core/casadi_misc.hpp"


using namespace std;
namespace casadi {

  extern "C"
  int CASADI_CONIC_PIQP_EXPORT
  casadi_register_conic_piqp(Conic::Plugin* plugin) {
    plugin->creator = PiqpInterface::creator;
    plugin->name = "piqp";
    plugin->doc = PiqpInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &PiqpInterface::options_;
    plugin->deserialize = &PiqpInterface::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_PIQP_EXPORT casadi_load_conic_piqp() {
    Conic::registerPlugin(casadi_register_conic_piqp);
  }

  // Map a kkt_solver option string to the PIQP enum
  inline piqp::KKTSolver kkt_solver_from_string(const std::string& s) {
    if (s == "dense_cholesky") return piqp::KKTSolver::dense_cholesky;
    if (s == "sparse_ldlt") return piqp::KKTSolver::sparse_ldlt;
    if (s == "sparse_ldlt_eq_cond") return piqp::KKTSolver::sparse_ldlt_eq_cond;
    if (s == "sparse_ldlt_ineq_cond") return piqp::KKTSolver::sparse_ldlt_ineq_cond;
    if (s == "sparse_ldlt_cond") return piqp::KKTSolver::sparse_ldlt_cond;
    if (s == "sparse_multistage") return piqp::KKTSolver::sparse_multistage;
    casadi_error("Unknown kkt_solver '" + s + "'. Choose one of: dense_cholesky, "
      "sparse_ldlt, sparse_ldlt_eq_cond, sparse_ldlt_ineq_cond, sparse_ldlt_cond, "
      "sparse_multistage.");
    return piqp::KKTSolver::sparse_ldlt;
  }

  PiqpInterface::PiqpInterface(const std::string& name,
                                   const std::map<std::string, Sparsity>& st)
    : Conic(name, st) {
    // Default to the sparse backend (PIQP's struct default is dense_cholesky)
    settings_.kkt_solver = piqp::KKTSolver::sparse_ldlt;
  }

  PiqpInterface::~PiqpInterface() {
    clear_mem();
  }

  const Options PiqpInterface::options_
  = {{&Conic::options_},
     {{"piqp",
       {OT_DICT,
        "const Options to be passed to piqp."}},
     }
  };

  void PiqpInterface::init(const Dict& opts) {
    // Initialize the base classes
    Conic::init(opts);

    // Read options
    for (auto&& op : opts) {
      if (op.first=="piqp") {
        const Dict& opts = op.second;
        for (auto&& op : opts) {
            if (op.first == "rho_init") {
                settings_.rho_init = op.second;
            } else if (op.first == "delta_init") {
                settings_.delta_init = op.second;
            } else if (op.first == "eps_abs") {
                settings_.eps_abs = op.second;
            } else if (op.first == "eps_rel") {
                settings_.eps_rel = op.second;
            } else if (op.first == "check_duality_gap") {
                settings_.check_duality_gap = op.second;
            } else if (op.first == "eps_duality_gap_abs") {
                settings_.eps_duality_gap_abs = op.second;
            } else if (op.first == "eps_duality_gap_rel") {
                settings_.eps_duality_gap_rel = op.second;
            } else if (op.first == "reg_lower_limit") {
                settings_.reg_lower_limit = op.second;
            } else if (op.first == "reg_finetune_lower_limit") {
                settings_.reg_finetune_lower_limit = op.second;
            } else if (op.first == "reg_finetune_primal_update_threshold") {
                settings_.reg_finetune_primal_update_threshold = static_cast<piqp::isize>(
                        op.second.to_int());
            } else if (op.first == "reg_finetune_dual_update_threshold") {
                settings_.reg_finetune_dual_update_threshold = static_cast<piqp::isize>(
                    op.second.to_int());
            } else if (op.first == "max_iter") {
                settings_.max_iter = static_cast<piqp::isize>(op.second.to_int());
            } else if (op.first == "max_factor_retires") {
                settings_.max_factor_retires = static_cast<piqp::isize>(
                        op.second.to_int());
            } else if (op.first == "preconditioner_scale_cost") {
                settings_.preconditioner_scale_cost = op.second;
            } else if (op.first == "preconditioner_iter") {
                settings_.preconditioner_iter = static_cast<piqp::isize>(
                        op.second.to_int());
            } else if (op.first == "tau") {
                settings_.tau = op.second;
            } else if (op.first == "iterative_refinement_always_enabled") {
                settings_.iterative_refinement_always_enabled = op.second;
            } else if (op.first == "iterative_refinement_eps_abs") {
                settings_.iterative_refinement_eps_abs = op.second;
            } else if (op.first == "iterative_refinement_eps_rel") {
                settings_.iterative_refinement_eps_rel = op.second;
            } else if (op.first == "iterative_refinement_max_iter") {
                settings_.iterative_refinement_max_iter = static_cast<piqp::isize>(
                        op.second.to_int());
            } else if (op.first == "iterative_refinement_min_improvement_rate") {
                settings_.iterative_refinement_min_improvement_rate = op.second;
            } else if (op.first == "iterative_refinement_static_regularization_eps") {
                settings_.iterative_refinement_static_regularization_eps = op.second;
            } else if (op.first == "iterative_refinement_static_regularization_rel") {
                settings_.iterative_refinement_static_regularization_rel = op.second;
            } else if (op.first == "verbose") {
                settings_.verbose = op.second;
            } else if (op.first == "compute_timings") {
                settings_.compute_timings = op.second;
            } else if (op.first == "kkt_solver") {
                settings_.kkt_solver = kkt_solver_from_string(op.second.to_string());
            } else {
              casadi_error("Unrecognised PIQP option '" + op.first + "'.");
            }
        }
      }
    }

    nnzH_ = H_.nnz();
    nnzA_ = A_.nnz()+nx_;

    alloc_w(nx_, true); // g
    alloc_w(nx_, true); // lbx
    alloc_w(nx_, true); // ubx
    alloc_w(na_, true); // lba
    alloc_w(na_, true); // uba
    alloc_w(nnzH_, true); // H
    alloc_w(nnzA_, true); // A
  }

  void PiqpInterface::finalize() {
    // dense_cholesky uses the dense backend; all sparse_* solvers use the sparse one
    sparse_backend_ = settings_.kkt_solver != piqp::KKTSolver::dense_cholesky;
    Conic::finalize();
  }

  int PiqpInterface::init_mem(void* mem) const {
    if (Conic::init_mem(mem)) return 1;
    auto m = static_cast<PiqpMemory*>(mem);

    m->tripletList.reserve(2 * H_.nnz());
    m->tripletListEq.reserve(na_);

    m->g_vector.resize(nx_);
    m->uba_vector.resize(na_);
    m->lba_vector.resize(na_);
    m->ubx_vector.resize(nx_);
    m->lbx_vector.resize(nx_);
    m->eq_b_vector.resize(na_);

    m->add_stat("preprocessing");
    m->add_stat("solver");
    m->add_stat("postprocessing");
    return 0;
  }

  int PiqpInterface::
  solve(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    typedef Eigen::Triplet<double> TripletT;

    auto m = static_cast<PiqpMemory*>(mem);
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

    // Split constraints into equality (lba == uba) and inequality (lba != uba)
    // PIQP 0.6.0+ supports double-sided inequality constraints natively: h_l <= Gx <= h_u
    const Eigen::Array<bool, Eigen::Dynamic, 1>
    is_equality = (m->uba_vector.array() == m->lba_vector.array()).eval();

    // Count equalities and inequalities, track mapping from original to new indices
    std::vector<unsigned int> number_of_prev_equality(na_, 0);
    std::vector<unsigned int> number_of_prev_inequality(na_, 0);
    std::vector<double> tmp_eq_vector;
    std::vector<double> tmp_ineq_lb_vector;
    std::vector<double> tmp_ineq_ub_vector;

    for (std::size_t k = 1; k < static_cast<std::size_t>(na_); ++k) {
      if (is_equality[k-1]) {
        number_of_prev_equality[k] = number_of_prev_equality[k-1] + 1;
        number_of_prev_inequality[k] = number_of_prev_inequality[k-1];
      } else {
        number_of_prev_equality[k] = number_of_prev_equality[k-1];
        number_of_prev_inequality[k] = number_of_prev_inequality[k-1] + 1;
      }
    }

    for (casadi_int k = 0; k < na_; ++k) {
      if (is_equality[k]) {
        tmp_eq_vector.push_back(m->lba_vector[k]);
      } else {
        tmp_ineq_lb_vector.push_back(m->lba_vector[k]);
        tmp_ineq_ub_vector.push_back(m->uba_vector[k]);
      }
    }

    m->eq_b_vector.resize(tmp_eq_vector.size());
    if (tmp_eq_vector.size() > 0) {
      m->eq_b_vector = Eigen::Map<Eigen::VectorXd>(
        get_ptr(tmp_eq_vector), tmp_eq_vector.size());
    }

    // For inequalities, use double-sided bounds directly (PIQP 0.6.0+ feature)
    Eigen::VectorXd ineq_lb_vector(tmp_ineq_lb_vector.size());
    Eigen::VectorXd ineq_ub_vector(tmp_ineq_ub_vector.size());
    if (tmp_ineq_lb_vector.size() > 0) {
      ineq_lb_vector = Eigen::Map<Eigen::VectorXd>(
        get_ptr(tmp_ineq_lb_vector), tmp_ineq_lb_vector.size());
      ineq_ub_vector = Eigen::Map<Eigen::VectorXd>(
        get_ptr(tmp_ineq_ub_vector), tmp_ineq_ub_vector.size());
    }

    std::size_t n_eq = m->eq_b_vector.size();
    std::size_t n_ineq = tmp_ineq_lb_vector.size();

    // Convert H_ from casadi::Sparsity to Eigen::SparseMatrix (misuse tripletList)
    H_.get_triplet(m->row, m->col);
    for (int k=0; k<H_.nnz(); ++k) {
      m->tripletList.push_back(PiqpMemory::TripletT(
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
      if (is_equality[m->row[k]]) {
        m->tripletListEq.push_back(TripletT(
          static_cast<double>(number_of_prev_equality[m->row[k]]),
          static_cast<double>(m->col[k]),
          static_cast<double>(A[k])));
      } else {
        // Inequality constraint - add to G matrix (PIQP uses h_l <= Gx <= h_u)
        m->tripletList.push_back(TripletT(
          static_cast<double>(number_of_prev_inequality[m->row[k]]),
          static_cast<double>(m->col[k]),
          static_cast<double>(A[k])));
      }
    }

    // Build equality constraint matrix A
    Eigen::SparseMatrix<double> A_spa(n_eq, nx_);
    A_spa.setFromTriplets(m->tripletListEq.begin(), m->tripletListEq.end());
    m->tripletListEq.clear();

    // Build inequality constraint matrix G (for h_l <= Gx <= h_u)
    Eigen::SparseMatrix<double> G_spa(n_ineq, nx_);
    G_spa.setFromTriplets(m->tripletList.begin(), m->tripletList.end());
    m->tripletList.clear();

    m->fstats.at("preprocessing").toc();

    // Solve Problem using PIQP 0.6.0+ API with double-sided inequality constraints
    // Problem form: min 0.5*x'*P*x + c'*x s.t. Ax = b, h_l <= Gx <= h_u, x_l <= x <= x_u
    m->fstats.at("solver").tic();
    if (sparse_backend_) {
        piqp::SparseSolver<double> solver;
        solver.settings() = settings_;

        solver.setup(
            H_spa, m->g_vector,
            A_spa, m->eq_b_vector,
            G_spa, ineq_lb_vector, ineq_ub_vector,
            m->lbx_vector, m->ubx_vector);
        m->status = solver.solve();

        m->results_x = std::make_unique<Eigen::VectorXd>(solver.result().x);
        m->results_y = std::make_unique<Eigen::VectorXd>(solver.result().y);
        // Inequality duals: z_u for upper bounds, z_l for lower bounds
        // Combined dual for lba <= Ax <= uba is z_u - z_l
        m->results_z = std::make_unique<Eigen::VectorXd>(
            solver.result().z_u - solver.result().z_l);
        // Box constraint duals: z_bu for upper, z_bl for lower
        m->results_lam_x = std::make_unique<Eigen::VectorXd>(
            solver.result().z_bu - solver.result().z_bl);
        m->objValue = solver.result().info.primal_obj;
    } else {
        piqp::DenseSolver<double> solver;
        solver.settings() = settings_;
        solver.setup(
            Eigen::MatrixXd(H_spa), m->g_vector,
            Eigen::MatrixXd(A_spa), m->eq_b_vector,
            Eigen::MatrixXd(G_spa), ineq_lb_vector, ineq_ub_vector,
            m->lbx_vector, m->ubx_vector);
        m->status = solver.solve();

        m->results_x = std::make_unique<Eigen::VectorXd>(solver.result().x);
        m->results_y = std::make_unique<Eigen::VectorXd>(solver.result().y);
        m->results_z = std::make_unique<Eigen::VectorXd>(
            solver.result().z_u - solver.result().z_l);
        m->results_lam_x = std::make_unique<Eigen::VectorXd>(
            solver.result().z_bu - solver.result().z_bl);
        m->objValue = solver.result().info.primal_obj;
    }
    m->fstats.at("solver").toc();

    // Post-processing to retrieve the results
    m->fstats.at("postprocessing").tic();
    casadi_copy(m->results_x->data(), nx_, res[CONIC_X]);
    casadi_copy(m->results_lam_x->data(), nx_, res[CONIC_LAM_X]);

    // Copy back the multipliers.
    // CasADi has LAM_X (multipliers for box constraints on x) and
    // LAM_A (multipliers for in- and equality constraints Ax).
    // PIQP returns: results_y (equality multipliers), results_z (inequality multipliers)
    if (n_ineq + n_eq > 0) {
        Eigen::VectorXd lam_a(na_);

        for (casadi_int k = 0; k < na_; ++k) {
          if (is_equality[k]) {
            lam_a[k] = m->results_y->coeff(number_of_prev_equality[k]);
          } else {
            lam_a[k] = m->results_z->coeff(number_of_prev_inequality[k]);
          }
        }
        casadi_copy(lam_a.data(), na_, res[CONIC_LAM_A]);
    }

    if (res[CONIC_COST]) {
      *res[CONIC_COST] = m->objValue;
    }

    m->d_qp.success = m->status == piqp::Status::PIQP_SOLVED;
    if (m->d_qp.success) {
      m->d_qp.unified_return_status = SOLVER_RET_SUCCESS;
    } else {
      if (m->status == piqp::Status::PIQP_MAX_ITER_REACHED) {
        m->d_qp.unified_return_status = SOLVER_RET_LIMITED;
      } else { // primal or dual infeasibility
        m->d_qp.unified_return_status = SOLVER_RET_UNKNOWN;
      }
    }
    m->fstats.at("postprocessing").toc();

    return 0;
  }

  Dict PiqpInterface::get_stats(void* mem) const {
    Dict stats = Conic::get_stats(mem);
    auto m = static_cast<PiqpMemory*>(mem);

    stats["return_status"] = piqp::status_to_string(m->status);
    return stats;
  }

  PiqpMemory::PiqpMemory() {
  }

  PiqpMemory::~PiqpMemory() {
  }

  PiqpInterface::PiqpInterface(DeserializingStream& s) : Conic(s) {
    std::size_t tmp;
    s.version("PiqpInterface", 1);
    s.unpack("PiqpInterface::nnzH", nnzH_);
    s.unpack("PiqpInterface::nnzA", nnzA_);
    s.unpack("PiqpInterface::settings::rho_init", settings_.rho_init);
    s.unpack("PiqpInterface::settings::delta_init", settings_.delta_init);
    s.unpack("PiqpInterface::settings::eps_abs", settings_.eps_abs);
    s.unpack("PiqpInterface::settings::eps_rel", settings_.eps_rel);
    s.unpack("PiqpInterface::settings::check_duality_gap", settings_.check_duality_gap);
    s.unpack("PiqpInterface::settings::eps_duality_gap_abs", settings_.eps_duality_gap_abs);
    s.unpack("PiqpInterface::settings::eps_duality_gap_rel", settings_.eps_duality_gap_rel);
    s.unpack("PiqpInterface::settings::reg_lower_limit", settings_.reg_lower_limit);
    s.unpack("PiqpInterface::settings::reg_finetune_lower_limit",
      settings_.reg_finetune_lower_limit);
    s.unpack("PiqpInterface::settings::reg_finetune_primal_update_threshold", tmp);
    settings_.reg_finetune_primal_update_threshold = tmp;
    s.unpack("PiqpInterface::settings::reg_finetune_dual_update_threshold", tmp);
    settings_.reg_finetune_dual_update_threshold = tmp;
    s.unpack("PiqpInterface::settings::max_iter", tmp);
    settings_.max_iter = tmp;
    s.unpack("PiqpInterface::settings::max_factor_retires", tmp);
    settings_.max_factor_retires = tmp;
    s.unpack("PiqpInterface::settings::preconditioner_scale_cost",
      settings_.preconditioner_scale_cost);
    s.unpack("PiqpInterface::settings::preconditioner_iter", tmp);
    settings_.preconditioner_iter = tmp;
    s.unpack("PiqpInterface::settings::tau", settings_.tau);
    s.unpack("PiqpInterface::settings::iterative_refinement_always_enabled",
      settings_.iterative_refinement_always_enabled);
    s.unpack("PiqpInterface::settings::iterative_refinement_eps_abs",
      settings_.iterative_refinement_eps_abs);
    s.unpack("PiqpInterface::settings::iterative_refinement_eps_rel",
      settings_.iterative_refinement_eps_rel);
    s.unpack("PiqpInterface::settings::iterative_refinement_max_iter", tmp);
    settings_.iterative_refinement_max_iter = tmp;
    s.unpack("PiqpInterface::settings::iterative_refinement_min_improvement_rate",
      settings_.iterative_refinement_min_improvement_rate);
    s.unpack("PiqpInterface::settings::iterative_refinement_static_regularization_eps",
      settings_.iterative_refinement_static_regularization_eps);
    s.unpack("PiqpInterface::settings::iterative_refinement_static_regularization_rel",
      settings_.iterative_refinement_static_regularization_rel);
    s.unpack("PiqpInterface::settings::verbose", settings_.verbose);
    s.unpack("PiqpInterface::settings::compute_timings", settings_.compute_timings);
    std::string kkt_solver;
    s.unpack("PiqpInterface::settings::kkt_solver", kkt_solver);
    settings_.kkt_solver = kkt_solver_from_string(kkt_solver);
  }

  void PiqpInterface::serialize_body(SerializingStream &s) const {
    std::size_t tmp;

    Conic::serialize_body(s);
    s.version("PiqpInterface", 1);
    s.pack("PiqpInterface::nnzH", nnzH_);
    s.pack("PiqpInterface::nnzA", nnzA_);
    s.pack("PiqpInterface::settings::rho_init", settings_.rho_init);
    s.pack("PiqpInterface::settings::delta_init", settings_.delta_init);
    s.pack("PiqpInterface::settings::eps_abs", settings_.eps_abs);
    s.pack("PiqpInterface::settings::eps_rel", settings_.eps_rel);
    s.pack("PiqpInterface::settings::check_duality_gap", settings_.check_duality_gap);
    s.pack("PiqpInterface::settings::eps_duality_gap_abs", settings_.eps_duality_gap_abs);
    s.pack("PiqpInterface::settings::eps_duality_gap_rel", settings_.eps_duality_gap_rel);
    s.pack("PiqpInterface::settings::reg_lower_limit", settings_.reg_lower_limit);
    s.pack("PiqpInterface::settings::reg_finetune_lower_limit",
      settings_.reg_finetune_lower_limit);
    tmp = settings_.reg_finetune_primal_update_threshold;
    s.pack("PiqpInterface::settings::reg_finetune_primal_update_threshold", tmp);
    tmp = settings_.reg_finetune_dual_update_threshold;
    s.pack("PiqpInterface::settings::reg_finetune_dual_update_threshold", tmp);
    tmp = settings_.max_iter;
    s.pack("PiqpInterface::settings::max_iter", tmp);
    tmp = settings_.max_factor_retires;
    s.pack("PiqpInterface::settings::max_factor_retires", tmp);
    s.pack("PiqpInterface::settings::preconditioner_scale_cost",
      settings_.preconditioner_scale_cost);
    tmp = settings_.preconditioner_iter;
    s.pack("PiqpInterface::settings::preconditioner_iter", tmp);
    s.pack("PiqpInterface::settings::tau", settings_.tau);
    s.pack("PiqpInterface::settings::iterative_refinement_always_enabled",
      settings_.iterative_refinement_always_enabled);
    s.pack("PiqpInterface::settings::iterative_refinement_eps_abs",
      settings_.iterative_refinement_eps_abs);
    s.pack("PiqpInterface::settings::iterative_refinement_eps_rel",
      settings_.iterative_refinement_eps_rel);
    tmp = settings_.iterative_refinement_max_iter;
    s.pack("PiqpInterface::settings::iterative_refinement_max_iter", tmp);
    s.pack("PiqpInterface::settings::iterative_refinement_min_improvement_rate",
      settings_.iterative_refinement_min_improvement_rate);
    s.pack("PiqpInterface::settings::iterative_refinement_static_regularization_eps",
      settings_.iterative_refinement_static_regularization_eps);
    s.pack("PiqpInterface::settings::iterative_refinement_static_regularization_rel",
      settings_.iterative_refinement_static_regularization_rel);
    s.pack("PiqpInterface::settings::verbose", settings_.verbose);
    s.pack("PiqpInterface::settings::compute_timings", settings_.compute_timings);
    s.pack("PiqpInterface::settings::kkt_solver",
      std::string(piqp::kkt_solver_to_string(settings_.kkt_solver)));
  }

} // namespace casadi
