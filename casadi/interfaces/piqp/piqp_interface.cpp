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

  PiqpInterface::PiqpInterface(const std::string& name,
                                   const std::map<std::string, Sparsity>& st)
    : Conic(name, st), sparse_backend(true) {

    has_refcount_ = true;
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
            } else if (op.first=="backend") {
              if (op.second == "sparse") {
                sparse_backend = true;
              } else if (op.second == "dense") {
                sparse_backend = false;
              } else {
                casadi_error("[Backend option] Please specify either sparse or dense");
              }
            } else {
              casadi_error("Not recognised");
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

  int PiqpInterface::init_mem(void* mem) const {
    if (Conic::init_mem(mem)) return 1;
    auto m = static_cast<PiqpMemory*>(mem);

    m->tripletList.reserve(2 * H_.nnz());
    m->tripletListEq.reserve(na_);

    m->g_vector.resize(nx_);
    m->uba_vector.resize(na_);
    m->lba_vector.resize(na_);
    m->ubx_vector.resize(na_);
    m->lbx_vector.resize(na_);
    m->ineq_b_vector.resize(na_);
    m->eq_b_vector.resize(na_);

    m->add_stat("preprocessing");
    m->add_stat("solver");
    m->add_stat("postprocessing");
    return 0;
  }

  inline const char* return_status_string(casadi_int status) {
    return "Unknown";
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

    // Use lhs_equals_rhs_constraint to split double-sided bounds into one-sided
    // bound for equality constraints and double-sided for inequality constraints
    const Eigen::Array<bool, Eigen::Dynamic, 1>
    lhs_equals_rhs_constraint = (m->uba_vector.array() == m->lba_vector.array()).eval();
    const Eigen::Array<bool, Eigen::Dynamic, 1>
    lhs_is_inf = m->lba_vector.array().isInf();
    const Eigen::Array<bool, Eigen::Dynamic, 1>
    rhs_is_inf = m->uba_vector.array().isInf();
    std::vector<unsigned int> number_of_prev_equality(lhs_equals_rhs_constraint.size(), 0);
    std::vector<unsigned int> number_of_prev_lb_inequality(lhs_equals_rhs_constraint.size(), 0);
    std::vector<unsigned int> number_of_prev_ub_inequality(lhs_equals_rhs_constraint.size(), 0);
    std::vector<double> tmp_eq_vector;
    std::vector<double> tmp_ineq_lb_vector;
    std::vector<double> tmp_ineq_ub_vector;
    {

      // number_of_prev_equality and number_of_prev_inequality are two vectors that contains
      // the number of equality and inequality that can be found before the current index
      // number_of_prev_equality[i] = number of equality that can be found before index i
      // number_of_prev_inequality[i] = number of inequality that can be found before index i
      // For instance:
      //     equality and inequality   [i, e, e, e, i, i, i, e]
      //     lhs_equals_rgs_contraint  [f, t, t, t, f, f, f, t]
      //     number_of_prev_equality   [0, 0, 1, 3, 3, 3, 3, 3]
      //     number_of_prev_lb_inequality [0, 1, 1, 1, 1, 2, 3, 4]
      //     number_of_prev_ub_inequality [0, 1, 1, 1, 1, 2, 3, 4]
      for (std::size_t k=1; k<lhs_equals_rhs_constraint.size(); ++k) {
        if (lhs_equals_rhs_constraint[k-1]) {
          number_of_prev_equality[k] = number_of_prev_equality[k-1] + 1;
          number_of_prev_lb_inequality[k] = number_of_prev_lb_inequality[k-1];
          number_of_prev_ub_inequality[k] = number_of_prev_ub_inequality[k-1];
        } else {
          number_of_prev_equality[k] = number_of_prev_equality[k-1];
          if (!lhs_is_inf[k-1]) {
              number_of_prev_lb_inequality[k] = number_of_prev_lb_inequality[k-1] +1;
          } else {
              number_of_prev_lb_inequality[k] = number_of_prev_lb_inequality[k-1];
          }
          if (!rhs_is_inf[k-1]) {
              number_of_prev_ub_inequality[k] = number_of_prev_ub_inequality[k-1] +1;
          } else {
              number_of_prev_ub_inequality[k] = number_of_prev_ub_inequality[k-1];
          }
        }
      }

      for (std::size_t k=0; k<lhs_equals_rhs_constraint.size(); ++k) {
        if (lhs_equals_rhs_constraint[k]) {
          tmp_eq_vector.push_back(m->lba_vector[k]);
        } else {
          if (!lhs_is_inf[k]) {
              tmp_ineq_lb_vector.push_back(m->lba_vector[k]);
          }
          if (!rhs_is_inf[k]) {
              tmp_ineq_ub_vector.push_back(m->uba_vector[k]);
          }
        }
      }

      m->eq_b_vector.resize(tmp_eq_vector.size());
      if (tmp_eq_vector.size() > 0) {
        m->eq_b_vector = Eigen::Map<Eigen::VectorXd>(
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
    std::size_t n_eq = m->eq_b_vector.size();
    std::size_t n_ineq_ub = m->uba_vector.size();
    std::size_t n_ineq_lb = m->lba_vector.size();
    std::size_t n_ineq = n_ineq_lb + n_ineq_ub;

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
      if (lhs_equals_rhs_constraint[m->row[k]]) {
        m->tripletListEq.push_back(TripletT(
          static_cast<double>(number_of_prev_equality[m->row[k]]),
          static_cast<double>(m->col[k]),
          static_cast<double>(A[k])));
      } else {
        if (!rhs_is_inf[m->row[k]]) {
          m->tripletList.push_back(TripletT(
            static_cast<double>(number_of_prev_ub_inequality[m->row[k]]),
            static_cast<double>(m->col[k]),
            static_cast<double>(A[k])));
        }
        if (!lhs_is_inf[m->row[k]]) {
          m->tripletList.push_back(TripletT(
            static_cast<double>(n_ineq_ub + number_of_prev_lb_inequality[m->row[k]]),
            static_cast<double>(m->col[k]),
            static_cast<double>(-A[k]))); // Reverse sign!
        }
      }
    }

    // Handle constraints on decision variable x in inequality constraint matrix C
    Eigen::SparseMatrix<double> A_spa(n_eq, nx_);
    A_spa.setFromTriplets(m->tripletListEq.begin(), m->tripletListEq.end());
    m->tripletListEq.clear();

    Eigen::SparseMatrix<double> C_spa(n_ineq, nx_);
    C_spa.setFromTriplets(m->tripletList.begin(), m->tripletList.end());
    m->tripletList.clear();

    // Get stacked lower and upper inequality bounds
    m->ineq_b_vector.resize(n_ineq);
    m->ineq_b_vector << m->uba_vector, -m->lba_vector;

    m->fstats.at("preprocessing").toc();

    // Solve Problem
    m->fstats.at("solver").tic();
    bool sparse_backend = true;
    if (sparse_backend) {
        piqp::SparseSolver<double> solver;
        solver.settings() = settings_;

        solver.setup(
            H_spa, m->g_vector,
            A_spa, m->eq_b_vector,
            C_spa, m->ineq_b_vector,
            m->lbx_vector, m->ubx_vector);
        m->status = solver.solve();

        m->results_x = std::make_unique<Eigen::VectorXd>(solver.result().x);
        m->results_y = std::make_unique<Eigen::VectorXd>(solver.result().y);
        m->results_z = std::make_unique<Eigen::VectorXd>(solver.result().z);
        m->results_lam_x = std::make_unique<Eigen::VectorXd>(
            solver.result().z_ub -
            solver.result().z_lb);
        m->objValue = solver.result().info.primal_obj;
    } else {
        piqp::DenseSolver<double> solver;
        solver.settings() = settings_;
        solver.setup(
            Eigen::MatrixXd(H_spa), m->g_vector,
            Eigen::MatrixXd(A_spa), m->eq_b_vector,
            Eigen::MatrixXd(C_spa), m->ineq_b_vector,
            m->lbx_vector, m->ubx_vector);
        m->status = solver.solve();

        m->results_x = std::make_unique<Eigen::VectorXd>(solver.result().x);
        m->results_y = std::make_unique<Eigen::VectorXd>(solver.result().y);
        m->results_z = std::make_unique<Eigen::VectorXd>(solver.result().z);
        m->results_lam_x = std::make_unique<Eigen::VectorXd>(
            solver.result().z_ub -
            solver.result().z_lb);
        m->objValue = solver.result().info.primal_obj;
    }
    m->fstats.at("solver").toc();

    // Post-processing to retrieve the results
    m->fstats.at("postprocessing").tic();
    casadi_copy(m->results_x->data(), nx_, res[CONIC_X]);
    casadi_copy(m->results_lam_x->data(), nx_, res[CONIC_LAM_X]);

    // Copy back the multipliers.
    // Note, casadi has LAM_X (multipliers for constraints on variable x) and
    // LAM_A (multipliers for in- and equality constraints) while proxqp has results_y
    // (equality multipliers) and results_z (inequality multipliers).
    if (n_ineq + n_eq > 0) {
        Eigen::VectorXd lam_a(na_);

        for (int k=0; k<lhs_equals_rhs_constraint.size(); ++k) {
          if (lhs_equals_rhs_constraint[k]) {
            lam_a[k] = m->results_y->coeff(number_of_prev_equality[k]);
          } else {
            lam_a[k] = 0;
              if (!lhs_is_inf[m->row[k]]) {
                lam_a[k] += m->results_z->coeff(number_of_prev_ub_inequality[k]);
              }
              if (!lhs_is_inf[m->row[k]]) {
                lam_a[k] -= m->results_z->coeff(number_of_prev_ub_inequality[k]);
              }
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

  void PiqpInterface::codegen_free_mem(CodeGenerator& g) const {
    g << "piqp_cleanup(" + codegen_mem(g) + ");\n";
  }

  void PiqpInterface::codegen_init_mem(CodeGenerator& g) const {
    // Sparsity Asp = vertcat(Sparsity::diag(nx_), A_);
    // casadi_int dummy_size = max(nx_+na_, max(Asp.nnz(), H_.nnz()));

    // g.local("A", "piqp_csc");
    // g.local("dummy[" + str(dummy_size) + "]", "casadi_real");
    // g << g.clear("dummy", dummy_size) << "\n";

    // g.constant_copy("A_row", Asp.get_row());
    // g.constant_copy("A_colind", Asp.get_colind());
    // g.constant_copy("H_row", H_.get_row());
    // g.constant_copy("H_colind", H_.get_colind());

    // g.local("A", "piqp_csc");
    // g << "A.m = " << nx_ + na_ << ";\n";
    // g << "A.n = " << nx_ << ";\n";
    // g << "A.nz = " << nnzA_ << ";\n";
    // g << "A.nzmax = " << nnzA_ << ";\n";
    // g << "A.x = dummy;\n";
    // g << "A.i = A_row;\n";
    // g << "A.p = A_colind;\n";

    // g.local("H", "piqp_csc");
    // g << "H.m = " << nx_ << ";\n";
    // g << "H.n = " << nx_ << ";\n";
    // g << "H.nz = " << H_.nnz() << ";\n";
    // g << "H.nzmax = " << H_.nnz() << ";\n";
    // g << "H.x = dummy;\n";
    // g << "H.i = H_row;\n";
    // g << "H.p = H_colind;\n";

    // g.local("data", "PPIPData");
    // g << "data.n = " << nx_ << ";\n";
    // g << "data.m = " << nx_ + na_ << ";\n";
    // g << "data.P = &H;\n";
    // g << "data.q = dummy;\n";
    // g << "data.A = &A;\n";
    // g << "data.l = dummy;\n";
    // g << "data.u = dummy;\n";

    // g.local("settings", "piqp_settings");
    // g << "piqp_set_default_settings(&settings);\n";
    // g << "settings.rho = " << settings_.rho << ";\n";
    // g << "settings.sigma = " << settings_.sigma << ";\n";
    // g << "settings.scaling = " << settings_.scaling << ";\n";
    // g << "settings.adaptive_rho = " << settings_.adaptive_rho << ";\n";
    // g << "settings.adaptive_rho_interval = " << settings_.adaptive_rho_interval << ";\n";
    // g << "settings.adaptive_rho_tolerance = " << settings_.adaptive_rho_tolerance << ";\n";
    // //g << "settings.adaptive_rho_fraction = " << settings_.adaptive_rho_fraction << ";\n";
    // g << "settings.max_iter = " << settings_.max_iter << ";\n";
    // g << "settings.eps_abs = " << settings_.eps_abs << ";\n";
    // g << "settings.eps_rel = " << settings_.eps_rel << ";\n";
    // g << "settings.eps_prim_inf = " << settings_.eps_prim_inf << ";\n";
    // g << "settings.eps_dual_inf = " << settings_.eps_dual_inf << ";\n";
    // g << "settings.alpha = " << settings_.alpha << ";\n";
    // g << "settings.delta = " << settings_.delta << ";\n";
    // g << "settings.polish = " << settings_.polish << ";\n";
    // g << "settings.polish_refine_iter = " << settings_.polish_refine_iter << ";\n";
    // g << "settings.verbose = " << settings_.verbose << ";\n";
    // g << "settings.scaled_termination = " << settings_.scaled_termination << ";\n";
    // g << "settings.check_termination = " << settings_.check_termination << ";\n";
    // //g << "settings.time_limit = " << settings_.time_limit << ";\n";

    // g << codegen_mem(g) + " = piqp_setup(&data, &settings);\n";
    // g << "return 0;\n";
  }

  void PiqpInterface::codegen_body(CodeGenerator& g) const {
    // g.add_include("piqp/piqp.h");
    // g.add_auxiliary(CodeGenerator::AUX_INF);

    // g.local("work", "piqp_workspace", "*");
    // g.init_local("work", codegen_mem(g));

    // g.comment("Set objective");
    // g.copy_default(g.arg(CONIC_G), nx_, "w", "0", false);
    // g << "if (piqp_update_lin_cost(work, w)) return 1;\n";

    // g.comment("Set bounds");
    // g.copy_default(g.arg(CONIC_LBX), nx_, "w", "-casadi_inf", false);
    // g.copy_default(g.arg(CONIC_LBA), na_, "w+"+str(nx_), "-casadi_inf", false);
    // g.copy_default(g.arg(CONIC_UBX), nx_, "w+"+str(nx_+na_), "casadi_inf", false);
    // g.copy_default(g.arg(CONIC_UBA), na_, "w+"+str(2*nx_+na_), "casadi_inf", false);
    // g << "if (piqp_update_bounds(work, w, w+" + str(nx_+na_)+ ")) return 1;\n";

    // g.comment("Project Hessian");
    // g << g.tri_project(g.arg(CONIC_H), H_, "w", false);

    // g.comment("Get constraint matrix");
    // std::string A_colind = g.constant(A_.get_colind());
    // g.local("offset", "casadi_int");
    // g.local("n", "casadi_int");
    // g.local("i", "casadi_int");
    // g << "offset = 0;\n";
    // g << "for (i=0; i< " << nx_ << "; ++i) {\n";
    // g << "w[" + str(nnzHupp_) + "+offset] = 1;\n";
    // g << "offset++;\n";
    // g << "n = " + A_colind + "[i+1]-" + A_colind + "[i];\n";
    // g << "casadi_copy(" << g.arg(CONIC_A) << "+" + A_colind + "[i], n, "
    //      "w+offset+" + str(nnzHupp_) + ");\n";
    // g << "offset+= n;\n";
    // g << "}\n";

    // g.comment("Pass Hessian and constraint matrices");
    // g << "if (piqp_update_P_A(work, w, 0, " + str(nnzHupp_) + ", w+" + str(nnzHupp_) +
    //      ", 0, " + str(nnzA_) + ")) return 1;\n";

    // g << "if (piqp_solve(work)) return 1;\n";

    // g.copy_check("&work->result->obj_val", 1, g.res(CONIC_COST), false, true);
    // g.copy_check("work->result->x", nx_, g.res(CONIC_X), false, true);
    // g.copy_check("work->result->y", nx_, g.res(CONIC_LAM_X), false, true);
    // g.copy_check("work->result->y+" + str(nx_), na_, g.res(CONIC_LAM_A), false, true);

    // g << "if (work->info->status_val != PIQP_SOLVED) return 1;\n";
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
    s.unpack("PiqpInterface::settings::reg_finetune_lower_limit", tmp);
    settings_.reg_finetune_lower_limit = tmp;
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
    s.unpack("PiqpInterface::backend", sparse_backend);
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
    s.pack("PiqpInterface::backend", sparse_backend);
  }

} // namespace casadi
