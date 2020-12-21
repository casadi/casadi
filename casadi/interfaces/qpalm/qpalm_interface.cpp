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


#include "qpalm_interface.hpp"
#include "casadi/core/casadi_misc.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_CONIC_QPALM_EXPORT
  casadi_register_conic_qpalm(Conic::Plugin* plugin) {
    plugin->creator = QpalmInterface::creator;
    plugin->name = "qpalm";
    plugin->doc = QpalmInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &QpalmInterface::options_;
    plugin->deserialize = &QpalmInterface::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_QPALM_EXPORT casadi_load_conic_qpalm() {
    Conic::registerPlugin(casadi_register_conic_qpalm);
  }

  QpalmInterface::QpalmInterface(const std::string& name,
                                   const std::map<std::string, Sparsity>& st)
    : Conic(name, st) {

    has_refcount_ = true;
  }

  QpalmInterface::~QpalmInterface() {
    clear_mem();
  }

  const Options QpalmInterface::options_
  = {{&Conic::options_},
     {{"qpalm",
       {OT_DICT,
        "const Options to be passed to qpalm."}},
      {"warm_start_primal",
       {OT_BOOL,
        "Use x0 input to warmstart [Default: true]."}},
      {"warm_start_dual",
       {OT_BOOL,
        "Use lam_a0 and lam_x0 input to warmstart [Default: truw]."}}
     }
  };

  void QpalmInterface::init(const Dict& opts) {
    // Initialize the base classes
    Conic::init(opts);

    qpalm_set_default_settings(&settings_);
    settings_.warm_start = false;

    warm_start_primal_ = true;
    warm_start_dual_ = true;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="warm_start_primal") {
        warm_start_primal_ = op.second;
      } else if (op.first=="warm_start_dual") {
        warm_start_dual_ = op.second;
      } else if (op.first=="qpalm") {
        const Dict& opts = op.second;
        for (auto&& op : opts) {
          if (op.first=="max_iter") {
            settings_.max_iter = op.second;
          } else if (op.first=="inner_max_iter") {
            settings_.inner_max_iter = op.second;
          } else if (op.first=="eps_abs") {
            settings_.eps_abs = op.second;
          } else if (op.first=="eps_rel") {
            settings_.eps_rel = op.second;
          } else if (op.first=="eps_abs_in") {
            settings_.eps_abs_in = op.second;
          } else if (op.first=="eps_rel_in") {
            settings_.eps_rel_in = op.second;
          } else if (op.first=="rho") {
            settings_.rho = op.second;
          } else if (op.first=="eps_prim_inf") {
            settings_.eps_prim_inf = op.second;
          } else if (op.first=="eps_dual_inf") {
            settings_.eps_dual_inf = op.second;
          } else if (op.first=="theta") {
            settings_.theta = op.second;
          } else if (op.first=="delta") {
            settings_.delta = op.second;
          } else if (op.first=="sigma_max") {
            settings_.sigma_max = op.second;
          } else if (op.first=="sigma_init") {
            settings_.sigma_init = op.second;
          } else if (op.first=="proximal") {
            settings_.proximal = op.second;
          } else if (op.first=="gamma_init") {
            settings_.gamma_init = op.second;
          } else if (op.first=="gamma_upd") {
            settings_.gamma_upd = op.second;
          } else if (op.first=="gamma_max") {
            settings_.gamma_max = op.second;
          } else if (op.first=="scaling") {
            settings_.scaling = op.second;
          } else if (op.first=="nonconvex") {
            settings_.nonconvex = op.second;
          } else if (op.first=="warm_start") {
            casadi_error("OSQP's warm_start option is impure and therefore disabled. "
                         "Use CasADi options 'warm_start_primal' and 'warm_start_dual' instead.");
          } else if (op.first=="verbose") {
            settings_.verbose = op.second;
          } else if (op.first=="print_iter") {
            settings_.print_iter = op.second;
          } else if (op.first=="reset_newton_iter") {
            settings_.reset_newton_iter = op.second;
          } else if (op.first=="enable_dual_termination") {
            settings_.enable_dual_termination = op.second;
          } else if (op.first=="dual_objective_limit") {
            settings_.dual_objective_limit = op.second;
          } else if (op.first=="time_limit") {
           settings_.time_limit = op.second;
          } else if (op.first=="ordering") {
            settings_.ordering = op.second;
          } else if (op.first=="factorization_method") {
            settings_.factorization_method = op.second;
          } else if (op.first=="max_rank_update") {
            settings_.max_rank_update = op.second;
          } else if (op.first=="max_rank_update_fraction") {
            settings_.max_rank_update_fraction = op.second;
          } else {
            casadi_error("Not recognised");
          }
        }
      }
    }

    nnzHupp_ = H_.nnz_upper();
    nnzA_ = A_.nnz()+nx_;

    alloc_w(nnzHupp_+nnzA_, false);
    alloc_w(2*nx_+2*na_, false);
  }

  int QpalmInterface::init_mem(void* mem) const {
    if (Conic::init_mem(mem)) return 1;
    auto m = static_cast<QpalmMemory*>(mem);

    Sparsity Asp = vertcat(Sparsity::diag(nx_), A_);
    std::vector<double> dummy(max(nx_+na_, max(Asp.nnz(), H_.nnz())));

    std::vector<c_int> A_row = vector_static_cast<c_int>(Asp.get_row());
    std::vector<c_int> A_colind = vector_static_cast<c_int>(Asp.get_colind());
    std::vector<c_int> H_row = vector_static_cast<c_int>(H_.get_row());
    std::vector<c_int> H_colind = vector_static_cast<c_int>(H_.get_colind());

    ladel_sparse_matrix A;
    A.nrow = nx_ + na_;
    A.ncol = nx_;
    A.nz = NULL;
    A.nzmax = nnzA_;
    A.x = get_ptr(dummy);
    A.i = get_ptr(A_row);
    A.p = get_ptr(A_colind);
    A.symmetry = UNSYMMETRIC;
    A.values = TRUE;

    ladel_sparse_matrix H;
    H.nrow = nx_;
    H.ncol = nx_;
    H.nz = NULL;
    H.nzmax = H_.nnz();
    H.x = get_ptr(dummy);
    H.i = get_ptr(H_row);
    H.p = get_ptr(H_colind);
    H.symmetry = UPPER;
    H.values = TRUE;

    QPALMData data;
    // Populate data
    data.n = nx_;
    data.m = nx_ + na_;
    // csc_matrix in mem
    data.Q = &H;
    data.q = get_ptr(dummy);
    data.c = 0;
    data.A = &A;
    data.bmin = get_ptr(dummy);
    data.bmax = get_ptr(dummy);

    // Setup workspace
    m->work = qpalm_setup(&data, &settings_);

    m->fstats["preprocessing"]  = FStats();
    m->fstats["solver"]         = FStats();
    m->fstats["postprocessing"] = FStats();
    return 0;
  }

  inline const char* return_status_string(casadi_int status) {
    return "Unknown";
  }

  int QpalmInterface::
  solve(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<QpalmMemory*>(mem);

    // Inputs
    const double *a=arg[CONIC_A],
                 *h=arg[CONIC_H],
                 *x_warm,
                 *y_warm;

    // Outputs
    double *x=res[CONIC_X],
           *cost=res[CONIC_COST],
           *lam_a=res[CONIC_LAM_A],
           *lam_x=res[CONIC_LAM_X];

    int ret;

    // Project Hessian
    casadi_tri_project(arg[CONIC_H], H_, w, false);

    // Get contraint matrix
    const casadi_int* colind = A_.colind();
    double* A = w + nnzHupp_;
    // Get constraint matrix
    casadi_int offset = 0;
    // Loop over columns
    for (casadi_int i=0; i<nx_; ++i) {
      A[offset] = 1;
      offset++;
      casadi_int n = colind[i+1]-colind[i];
      casadi_copy(a+colind[i], n, A+offset);
      offset+= n;
    }

    qpalm_update_Q_A(m->work, w, A);

    // Set bounds
    casadi_copy(arg[CONIC_LBX], nx_, w);
    casadi_copy(arg[CONIC_LBA], na_, w+nx_);
    casadi_copy(arg[CONIC_UBX], nx_, w+nx_+na_);
    casadi_copy(arg[CONIC_UBA], na_, w+2*nx_+na_);
    
    // Set objective
    if (arg[CONIC_G]) {
      qpalm_update_q(m->work, arg[CONIC_G]);
      // casadi_assert(ret==0, "Problem in qpalm_update_q");
    }

    qpalm_update_bounds(m->work, w, w+nx_+na_);
    // casadi_assert(ret==0, "Problem in qpalm_update_bounds");

    if (warm_start_primal_) {
      x_warm=arg[CONIC_X0];
    } 
    else
    {
      x_warm=NULL;
    }

    if (warm_start_dual_) {
      casadi_copy(arg[CONIC_LAM_X0], nx_, w);
      casadi_copy(arg[CONIC_LAM_A0], na_, w+nx_);
      y_warm = w;
    } 
    else
    {
      y_warm = NULL;
    }
    
    qpalm_warm_start(m->work, x_warm, y_warm);

    // Solve Problem
    qpalm_solve(m->work);
    // casadi_assert(ret==0, "Problem in qpalm_solve");

    casadi_copy(m->work->solution->x, nx_, res[CONIC_X]);
    casadi_copy(m->work->solution->y, nx_, res[CONIC_LAM_X]);
    casadi_copy(m->work->solution->y+nx_, na_, res[CONIC_LAM_A]);
    if (res[CONIC_COST]) *res[CONIC_COST] = m->work->info->objective;

    m->success = m->work->info->status_val == QPALM_SOLVED;
    if (m->success) m->unified_return_status = SOLVER_RET_SUCCESS;

    return 0;
  }

  void QpalmInterface::codegen_free_mem(CodeGenerator& g) const {
    g << "qpalm_cleanup(" + codegen_mem(g) + ");\n";
  }

  void QpalmInterface::codegen_init_mem(CodeGenerator& g) const {
    Sparsity Asp = vertcat(Sparsity::diag(nx_), A_);
    casadi_int dummy_size = max(nx_+na_, max(Asp.nnz(), H_.nnz()));

    // g.local("A", "ladel_sparse_matrix");
    g.local("dummy[" + str(dummy_size) + "]", "casadi_real");
    g << g.clear("dummy", dummy_size) << "\n";

    g.constant_copy("A_row", Asp.get_row());
    g.constant_copy("A_colind", Asp.get_colind());
    g.constant_copy("H_row", H_.get_row());
    g.constant_copy("H_colind", H_.get_colind());

    g.local("A", "ladel_sparse_matrix");
    g << "A.nrow = " << nx_ + na_ << ";\n";
    g << "A.ncol = " << nx_ << ";\n";
    g << "A.nz = " << "NULL" << ";\n";
    g << "A.nzmax = " << nnzA_ << ";\n";
    g << "A.x = dummy;\n";
    g << "A.i = A_row;\n";
    g << "A.p = A_colind;\n";
    g << "A.symmetry = UNSYMMETRIC;\n";

    g.local("H", "ladel_sparse_matrix");
    g << "H.nrow = " << nx_ << ";\n";
    g << "H.ncol = " << nx_ << ";\n";
    g << "H.nz = " << "NULL" << ";\n";
    g << "H.nzmax = " << H_.nnz() << ";\n";
    g << "H.x = dummy;\n";
    g << "H.i = H_row;\n";
    g << "H.p = H_colind;\n";
    g << "H.symmetry = UPPER;\n";

    g.local("data", "QPALMData");
    g << "data.n = " << nx_ << ";\n";
    g << "data.m = " << nx_ + na_ << ";\n";
    g << "data.Q = &H;\n";
    g << "data.q = dummy;\n";
    g << "data.c = 0;\n";
    g << "data.A = &A;\n";
    g << "data.bmin = dummy;\n";
    g << "data.bmax = dummy;\n";

    g.local("settings", "QPALMSettings");
    g << "qpalm_set_default_settings(&settings);\n";
    g << "settings.max_iter = " << settings_.max_iter << ";\n";
    g << "settings.inner_max_iter = " << settings_.inner_max_iter << ";\n";
    g << "settings.eps_abs = " << settings_.eps_abs << ";\n";
    g << "settings.eps_rel = " << settings_.eps_rel << ";\n";
    g << "settings.eps_abs_in = " << settings_.eps_abs_in << ";\n";
    g << "settings.eps_rel_in = " << settings_.eps_rel_in << ";\n";
    g << "settings.rho = " << settings_.rho << ";\n";
    g << "settings.eps_prim_inf = " << settings_.eps_prim_inf << ";\n";
    g << "settings.eps_dual_inf = " << settings_.eps_dual_inf << ";\n";
    g << "settings.theta = " << settings_.theta << ";\n";
    g << "settings.delta = " << settings_.delta << ";\n";
    g << "settings.sigma_max = " << settings_.sigma_max << ";\n";
    g << "settings.sigma_init = " << settings_.sigma_init << ";\n";
    g << "settings.proximal = " << settings_.proximal << ";\n";
    g << "settings.gamma_init = " << settings_.gamma_init << ";\n";
    g << "settings.gamma_upd = " << settings_.gamma_upd << ";\n";
    g << "settings.gamma_max = " << settings_.gamma_max << ";\n";
    g << "settings.scaling = " << settings_.scaling << ";\n";
    g << "settings.nonconvex = " << settings_.nonconvex << ";\n";
    g << "settings.warm_start = " << settings_.warm_start << ";\n";
    g << "settings.verbose = " << settings_.verbose << ";\n";
    g << "settings.print_iter = " << settings_.print_iter << ";\n";
    g << "settings.reset_newton_iter = " << settings_.reset_newton_iter << ";\n";
    g << "settings.enable_dual_termination = " << settings_.enable_dual_termination << ";\n";
    g << "settings.dual_objective_limit = " << settings_.dual_objective_limit << ";\n";
    g << "settings.time_limit = " << settings_.time_limit << ";\n";
    g << "settings.ordering = " << settings_.ordering << ";\n";
    g << "settings.factorization_method = " << settings_.factorization_method << ";\n";
    g << "settings.max_rank_update = " << settings_.max_rank_update << ";\n";
    g << "settings.max_rank_update_fraction = " << settings_.max_rank_update_fraction << ";\n";

    g << codegen_mem(g) + " = qpalm_setup(&data, &settings);\n";
    g << "return 0;\n";
  }

  void QpalmInterface::codegen_body(CodeGenerator& g) const {
    g.add_include("qpalm/qpalm.h");
    g.add_auxiliary(CodeGenerator::AUX_INF);

    g.local("work", "QPALMWorkspace", "*");
    g.init_local("work", codegen_mem(g));

    g.comment("Project Hessian");
    g << g.tri_project(g.arg(CONIC_H), H_, "w", false);

    g.comment("Get constraint matrix");
    std::string A_colind = g.constant(A_.get_colind());
    g.local("offset", "casadi_int");
    g.local("n", "casadi_int");
    g.local("i", "casadi_int");
    g << "offset = 0;\n";
    g << "for (i=0; i< " << nx_ << "; ++i) {\n";
    g << "w[" + str(nnzHupp_) + "+offset] = 1;\n";
    g << "offset++;\n";
    g << "n = " + A_colind + "[i+1]-" + A_colind + "[i];\n";
    g << "casadi_copy(" << g.arg(CONIC_A) << "+" + A_colind + "[i], n, "
         "w+offset+" + str(nnzHupp_) + ");\n";
    g << "offset+= n;\n";
    g << "}\n";

    g.comment("Pass Hessian and constraint matrices");
    g << "qpalm_update_Q_A(work, w, A);\n";

    g.comment("Set objective");
    g.copy_default(g.arg(CONIC_G), nx_, "w", "0", false);
    g << "qpalm_update_q(work, w);\n";

    g.comment("Set bounds");
    g.copy_default(g.arg(CONIC_LBX), nx_, "w", "-casadi_inf", false);
    g.copy_default(g.arg(CONIC_LBA), na_, "w+"+str(nx_), "-casadi_inf", false);
    g.copy_default(g.arg(CONIC_UBX), nx_, "w+"+str(nx_+na_), "casadi_inf", false);
    g.copy_default(g.arg(CONIC_UBA), na_, "w+"+str(2*nx_+na_), "casadi_inf", false);
    g << "qpalm_update_bounds(work, w, w+" + str(nx_+na_)+ ");\n";
  
    g.copy_default(g.arg(CONIC_LAM_X0), nx_, "w", "0", false);
    g.copy_default(g.arg(CONIC_LAM_A0), na_, "w+"+str(nx_), "0", false);
    g << "qpalm_warm_start(work," + g.arg(CONIC_X0) + ", w);\n";

    g << "qpalm_solve(work);\n";

    g.copy_check("&work->info->objective", 1, g.res(CONIC_COST), false, true);
    g.copy_check("work->solution->x", nx_, g.res(CONIC_X), false, true);
    g.copy_check("work->solution->y", nx_, g.res(CONIC_LAM_X), false, true);
    g.copy_check("work->solution->y+" + str(nx_), na_, g.res(CONIC_LAM_A), false, true);

    g << "if (work->info->status_val != QPALM_SOLVED) return 1;\n";
  }

  Dict QpalmInterface::get_stats(void* mem) const {
    Dict stats = Conic::get_stats(mem);
    auto m = static_cast<QpalmMemory*>(mem);
    stats["return_status"] = m->work->info->status;
    return stats;
  }

  QpalmMemory::QpalmMemory() {
  }

  QpalmMemory::~QpalmMemory() {
    qpalm_cleanup(work);
  }

  QpalmInterface::QpalmInterface(DeserializingStream& s) : Conic(s) {
    s.version("QpalmInterface", 1);
    s.unpack("QpalmInterface::nnzHupp", nnzHupp_);
    s.unpack("QpalmInterface::nnzA", nnzA_);
    s.unpack("QpalmInterface::warm_start_primal", warm_start_primal_);
    s.unpack("QpalmInterface::warm_start_dual", warm_start_dual_);

    qpalm_set_default_settings(&settings_);
    s.unpack("QpalmInterface::settings::max_iter", settings_.max_iter);
    s.unpack("QpalmInterface::settings::inner_max_iter", settings_.inner_max_iter);
    s.unpack("QpalmInterface::settings::eps_abs", settings_.eps_abs);
    s.unpack("QpalmInterface::settings::eps_rel", settings_.eps_rel);
    s.unpack("QpalmInterface::settings::eps_abs_in", settings_.eps_abs_in);
    s.unpack("QpalmInterface::settings::eps_rel_in", settings_.eps_rel_in);
    s.unpack("QpalmInterface::settings::rho", settings_.rho);
    s.unpack("QpalmInterface::settings::eps_prim_inf", settings_.eps_prim_inf);
    s.unpack("QpalmInterface::settings::eps_dual_inf", settings_.eps_dual_inf);
    s.unpack("QpalmInterface::settings::theta", settings_.theta);
    s.unpack("QpalmInterface::settings::delta", settings_.delta);
    s.unpack("QpalmInterface::settings::sigma_max", settings_.sigma_max);
    s.unpack("QpalmInterface::settings::sigma_init", settings_.sigma_init);
    s.unpack("QpalmInterface::settings::proximal", settings_.proximal);
    s.unpack("QpalmInterface::settings::gamma_init", settings_.gamma_init);
    s.unpack("QpalmInterface::settings::gamma_upd", settings_.gamma_upd);
    s.unpack("QpalmInterface::settings::gamma_max", settings_.gamma_max);
    s.unpack("QpalmInterface::settings::scaling", settings_.scaling);
    s.unpack("QpalmInterface::settings::nonconvex", settings_.nonconvex);
    s.unpack("QpalmInterface::settings::warm_start", settings_.warm_start);
    s.unpack("QpalmInterface::settings::verbose", settings_.verbose);
    s.unpack("QpalmInterface::settings::print_iter", settings_.print_iter);
    s.unpack("QpalmInterface::settings::reset_newton_iter", settings_.reset_newton_iter);
    s.unpack("QpalmInterface::settings::enable_dual_termination", settings_.enable_dual_termination);
    s.unpack("QpalmInterface::settings::dual_objective_limit", settings_.dual_objective_limit);
    s.unpack("QpalmInterface::settings::time_limit", settings_.time_limit);
    s.unpack("QpalmInterface::settings::ordering", settings_.ordering);
    s.unpack("QpalmInterface::settings::factorization_method", settings_.factorization_method);
    s.unpack("QpalmInterface::settings::max_rank_update", settings_.max_rank_update);
    s.unpack("QpalmInterface::settings::max_rank_update_fraction", settings_.max_rank_update_fraction);
  }

  void QpalmInterface::serialize_body(SerializingStream &s) const {
    Conic::serialize_body(s);
    s.version("QpalmInterface", 1);
    s.pack("QpalmInterface::nnzHupp", nnzHupp_);
    s.pack("QpalmInterface::nnzA", nnzA_);
    s.pack("QpalmInterface::warm_start_primal", warm_start_primal_);
    s.pack("QpalmInterface::warm_start_dual", warm_start_dual_);
    s.pack("QpalmInterface::settings::max_iter", settings_.max_iter);
    s.pack("QpalmInterface::settings::inner_max_iter", settings_.inner_max_iter);
    s.pack("QpalmInterface::settings::eps_abs", settings_.eps_abs);
    s.pack("QpalmInterface::settings::eps_rel", settings_.eps_rel);
    s.pack("QpalmInterface::settings::eps_abs_in", settings_.eps_abs_in);
    s.pack("QpalmInterface::settings::eps_rel_in", settings_.eps_rel_in);
    s.pack("QpalmInterface::settings::rho", settings_.rho);
    s.pack("QpalmInterface::settings::eps_prim_inf", settings_.eps_prim_inf);
    s.pack("QpalmInterface::settings::eps_dual_inf", settings_.eps_dual_inf);
    s.pack("QpalmInterface::settings::theta", settings_.theta);
    s.pack("QpalmInterface::settings::delta", settings_.delta);
    s.pack("QpalmInterface::settings::sigma_max", settings_.sigma_max);
    s.pack("QpalmInterface::settings::sigma_init", settings_.sigma_init);
    s.pack("QpalmInterface::settings::proximal", settings_.proximal);
    s.pack("QpalmInterface::settings::gamma_init", settings_.gamma_init);
    s.pack("QpalmInterface::settings::gamma_upd", settings_.gamma_upd);
    s.pack("QpalmInterface::settings::gamma_max", settings_.gamma_max);
    s.pack("QpalmInterface::settings::scaling", settings_.scaling);
    s.pack("QpalmInterface::settings::nonconvex", settings_.nonconvex);
    s.pack("QpalmInterface::settings::warm_start", settings_.warm_start);
    s.pack("QpalmInterface::settings::verbose", settings_.verbose);
    s.pack("QpalmInterface::settings::print_iter", settings_.print_iter);
    s.pack("QpalmInterface::settings::reset_newton_iter", settings_.reset_newton_iter);
    s.pack("QpalmInterface::settings::enable_dual_termination", settings_.enable_dual_termination);
    s.pack("QpalmInterface::settings::dual_objective_limit", settings_.dual_objective_limit);
    s.pack("QpalmInterface::settings::time_limit", settings_.time_limit);
    s.pack("QpalmInterface::settings::ordering", settings_.ordering);
    s.pack("QpalmInterface::settings::factorization_method", settings_.factorization_method);
    s.pack("QpalmInterface::settings::max_rank_update", settings_.max_rank_update);
    s.pack("QpalmInterface::settings::max_rank_update_fraction", settings_.max_rank_update_fraction);
  }

} // namespace casadi
