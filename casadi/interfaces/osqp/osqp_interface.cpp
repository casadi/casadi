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


#include "osqp_interface.hpp"
#include "casadi/core/casadi_misc.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_CONIC_OSQP_EXPORT
  casadi_register_conic_osqp(Conic::Plugin* plugin) {
    plugin->creator = OsqpInterface::creator;
    plugin->name = "osqp";
    plugin->doc = OsqpInterface::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &OsqpInterface::options_;
    plugin->deserialize = &OsqpInterface::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_OSQP_EXPORT casadi_load_conic_osqp() {
    Conic::registerPlugin(casadi_register_conic_osqp);
  }

  OsqpInterface::OsqpInterface(const std::string& name,
                                   const std::map<std::string, Sparsity>& st)
    : Conic(name, st) {

    has_refcount_ = true;
  }

  OsqpInterface::~OsqpInterface() {
    clear_mem();
  }

  const Options OsqpInterface::options_
  = {{&Conic::options_},
     {{"osqp",
       {OT_DICT,
        "const Options to be passed to osqp."}},
      {"warm_start_primal",
       {OT_BOOL,
        "Use x0 input to warmstart [Default: true]."}},
      {"warm_start_dual",
       {OT_BOOL,
        "Use lam_a0 and lam_x0 input to warmstart [Default: truw]."}}
     }
  };

  void OsqpInterface::init(const Dict& opts) {
    // Initialize the base classes
    Conic::init(opts);

    osqp_set_default_settings(&settings_);
    settings_.warm_start = false;

    warm_start_primal_ = true;
    warm_start_dual_ = true;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="warm_start_primal") {
        warm_start_primal_ = op.second;
      } else if (op.first=="warm_start_dual") {
        warm_start_dual_ = op.second;
      } else if (op.first=="osqp") {
        const Dict& opts = op.second;
        for (auto&& op : opts) {
          if (op.first=="rho") {
            settings_.rho = op.second;
          } else if (op.first=="sigma") {
            settings_.sigma = op.second;
          } else if (op.first=="scaling") {
            settings_.scaling = op.second;
          } else if (op.first=="adaptive_rho") {
            settings_.adaptive_rho = op.second;
          } else if (op.first=="adaptive_rho_interval") {
            settings_.adaptive_rho_interval = op.second;
          } else if (op.first=="adaptive_rho_tolerance") {
            settings_.adaptive_rho_tolerance = op.second;
          //} else if (op.first=="adaptive_rho_fraction") {
          //  settings_.adaptive_rho_fraction = op.second;
          } else if (op.first=="max_iter") {
            settings_.max_iter = op.second;
          } else if (op.first=="eps_abs") {
            settings_.eps_abs = op.second;
          } else if (op.first=="eps_rel") {
            settings_.eps_rel = op.second;
          } else if (op.first=="eps_prim_inf") {
            settings_.eps_prim_inf = op.second;
          } else if (op.first=="eps_dual_inf") {
            settings_.eps_dual_inf = op.second;
          } else if (op.first=="alpha") {
            settings_.alpha = op.second;
          } else if (op.first=="delta") {
            settings_.delta = op.second;
          } else if (op.first=="polish") {
            settings_.polish = op.second;
          } else if (op.first=="polish_refine_iter") {
            settings_.polish_refine_iter = op.second;
          } else if (op.first=="verbose") {
            settings_.verbose = op.second;
          } else if (op.first=="scaled_termination") {
            settings_.scaled_termination = op.second;
          } else if (op.first=="check_termination") {
            settings_.check_termination = op.second;
          } else if (op.first=="warm_start") {
            casadi_error("OSQP's warm_start option is impure and therefore disabled. "
                         "Use CasADi options 'warm_start_primal' and 'warm_start_dual' instead.");
          //} else if (op.first=="time_limit") {
          //  settings_.time_limit = op.second;
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

  int OsqpInterface::init_mem(void* mem) const {
    if (Conic::init_mem(mem)) return 1;
    auto m = static_cast<OsqpMemory*>(mem);

    Sparsity Asp = vertcat(Sparsity::diag(nx_), A_);
    std::vector<double> dummy(max(nx_+na_, max(Asp.nnz(), H_.nnz())));

    std::vector<c_int> A_row = vector_static_cast<c_int>(Asp.get_row());
    std::vector<c_int> A_colind = vector_static_cast<c_int>(Asp.get_colind());
    std::vector<c_int> H_row = vector_static_cast<c_int>(H_.get_row());
    std::vector<c_int> H_colind = vector_static_cast<c_int>(H_.get_colind());

    csc A;
    A.m = nx_ + na_;
    A.n = nx_;
    A.nz = nnzA_;
    A.nzmax = A.nz;
    A.x = get_ptr(dummy);
    A.i = get_ptr(A_row);
    A.p = get_ptr(A_colind);

    csc H;
    H.m = nx_;
    H.n = nx_;
    H.nz = H_.nnz();
    H.nzmax = H_.nnz();
    H.x = get_ptr(dummy);
    H.i = get_ptr(H_row);
    H.p = get_ptr(H_colind);

    OSQPData data;
    // Populate data
    data.n = nx_;
    data.m = nx_ + na_;
    // csc_matrix in mem
    data.P = &H;
    data.q = get_ptr(dummy);
    data.A = &A;
    data.l = get_ptr(dummy);
    data.u = get_ptr(dummy);

    // Setup workspace
    m->work = osqp_setup(&data, &settings_);

    m->fstats["preprocessing"]  = FStats();
    m->fstats["solver"]         = FStats();
    m->fstats["postprocessing"] = FStats();
    return 0;
  }

  inline const char* return_status_string(casadi_int status) {
    return "Unknown";
  }

  int OsqpInterface::
  solve(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<OsqpMemory*>(mem);

    // Inputs
    const double *a=arg[CONIC_A],
                 *h=arg[CONIC_H];

    // Outputs
    double *x=res[CONIC_X],
           *cost=res[CONIC_COST],
           *lam_a=res[CONIC_LAM_A],
           *lam_x=res[CONIC_LAM_X];

    int ret;

    // Set objective
    if (arg[CONIC_G]) {
      ret = osqp_update_lin_cost(m->work, arg[CONIC_G]);
      casadi_assert(ret==0, "Problem in osqp_update_lin_cost");
    }

    // Set bounds
    casadi_copy(arg[CONIC_LBX], nx_, w);
    casadi_copy(arg[CONIC_LBA], na_, w+nx_);
    casadi_copy(arg[CONIC_UBX], nx_, w+nx_+na_);
    casadi_copy(arg[CONIC_UBA], na_, w+2*nx_+na_);

    ret = osqp_update_bounds(m->work, w, w+nx_+na_);
    casadi_assert(ret==0, "Problem in osqp_update_bounds");

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

    // Pass Hessian and constraint matrices
    ret = osqp_update_P_A(m->work, w, nullptr, nnzHupp_, A, nullptr, nnzA_);
    casadi_assert(ret==0, "Problem in osqp_update_P_A");


    if (warm_start_primal_) {
      ret = osqp_warm_start_x(m->work, arg[CONIC_X0]);
      casadi_assert(ret==0, "Problem in osqp_warm_start_x");
    }

    if (warm_start_dual_) {
      casadi_copy(arg[CONIC_LAM_X0], nx_, w);
      casadi_copy(arg[CONIC_LAM_A0], na_, w+nx_);
      ret = osqp_warm_start_y(m->work, w);
      casadi_assert(ret==0, "Problem in osqp_warm_start_y");
    }

    // Solve Problem
    ret = osqp_solve(m->work);
    casadi_assert(ret==0, "Problem in osqp_solve");

    casadi_copy(m->work->solution->x, nx_, res[CONIC_X]);
    casadi_copy(m->work->solution->y, nx_, res[CONIC_LAM_X]);
    casadi_copy(m->work->solution->y+nx_, na_, res[CONIC_LAM_A]);
    if (res[CONIC_COST]) *res[CONIC_COST] = m->work->info->obj_val;

    m->success = m->work->info->status_val == OSQP_SOLVED;
    if (m->success) m->unified_return_status = SOLVER_RET_SUCCESS;

    return 0;
  }

  void OsqpInterface::codegen_free_mem(CodeGenerator& g) const {
    g << "osqp_cleanup(" + codegen_mem(g) + ");\n";
  }

  void OsqpInterface::codegen_init_mem(CodeGenerator& g) const {
    Sparsity Asp = vertcat(Sparsity::diag(nx_), A_);
    casadi_int dummy_size = max(nx_+na_, max(Asp.nnz(), H_.nnz()));

    g.local("A", "csc");
    g.local("dummy[" + str(dummy_size) + "]", "casadi_real");
    g << g.clear("dummy", dummy_size) << "\n";

    g.constant_copy("A_row", Asp.get_row());
    g.constant_copy("A_colind", Asp.get_colind());
    g.constant_copy("H_row", H_.get_row());
    g.constant_copy("H_colind", H_.get_colind());

    g.local("A", "csc");
    g << "A.m = " << nx_ + na_ << ";\n";
    g << "A.n = " << nx_ << ";\n";
    g << "A.nz = " << nnzA_ << ";\n";
    g << "A.nzmax = " << nnzA_ << ";\n";
    g << "A.x = dummy;\n";
    g << "A.i = A_row;\n";
    g << "A.p = A_colind;\n";

    g.local("H", "csc");
    g << "H.m = " << nx_ << ";\n";
    g << "H.n = " << nx_ << ";\n";
    g << "H.nz = " << H_.nnz() << ";\n";
    g << "H.nzmax = " << H_.nnz() << ";\n";
    g << "H.x = dummy;\n";
    g << "H.i = H_row;\n";
    g << "H.p = H_colind;\n";

    g.local("data", "OSQPData");
    g << "data.n = " << nx_ << ";\n";
    g << "data.m = " << nx_ + na_ << ";\n";
    g << "data.P = &H;\n";
    g << "data.q = dummy;\n";
    g << "data.A = &A;\n";
    g << "data.l = dummy;\n";
    g << "data.u = dummy;\n";

    g.local("settings", "OSQPSettings");
    g << "osqp_set_default_settings(&settings);\n";
    g << "settings.rho = " << settings_.rho << ";\n";
    g << "settings.sigma = " << settings_.sigma << ";\n";
    g << "settings.scaling = " << settings_.scaling << ";\n";
    g << "settings.adaptive_rho = " << settings_.adaptive_rho << ";\n";
    g << "settings.adaptive_rho_interval = " << settings_.adaptive_rho_interval << ";\n";
    g << "settings.adaptive_rho_tolerance = " << settings_.adaptive_rho_tolerance << ";\n";
    //g << "settings.adaptive_rho_fraction = " << settings_.adaptive_rho_fraction << ";\n";
    g << "settings.max_iter = " << settings_.max_iter << ";\n";
    g << "settings.eps_abs = " << settings_.eps_abs << ";\n";
    g << "settings.eps_rel = " << settings_.eps_rel << ";\n";
    g << "settings.eps_prim_inf = " << settings_.eps_prim_inf << ";\n";
    g << "settings.eps_dual_inf = " << settings_.eps_dual_inf << ";\n";
    g << "settings.alpha = " << settings_.alpha << ";\n";
    g << "settings.delta = " << settings_.delta << ";\n";
    g << "settings.polish = " << settings_.polish << ";\n";
    g << "settings.polish_refine_iter = " << settings_.polish_refine_iter << ";\n";
    g << "settings.verbose = " << settings_.verbose << ";\n";
    g << "settings.scaled_termination = " << settings_.scaled_termination << ";\n";
    g << "settings.check_termination = " << settings_.check_termination << ";\n";
    g << "settings.warm_start = " << settings_.warm_start << ";\n";
    //g << "settings.time_limit = " << settings_.time_limit << ";\n";

    g << codegen_mem(g) + " = osqp_setup(&data, &settings);\n";
    g << "return 0;\n";
  }

  void OsqpInterface::codegen_body(CodeGenerator& g) const {
    g.add_include("osqp/osqp.h");
    g.add_auxiliary(CodeGenerator::AUX_INF);

    g.local("work", "OSQPWorkspace", "*");
    g.init_local("work", codegen_mem(g));

    g.comment("Set objective");
    g.copy_default(g.arg(CONIC_G), nx_, "w", "0", false);
    g << "if (osqp_update_lin_cost(work, w)) return 1;\n";

    g.comment("Set bounds");
    g.copy_default(g.arg(CONIC_LBX), nx_, "w", "-casadi_inf", false);
    g.copy_default(g.arg(CONIC_LBA), na_, "w+"+str(nx_), "-casadi_inf", false);
    g.copy_default(g.arg(CONIC_UBX), nx_, "w+"+str(nx_+na_), "casadi_inf", false);
    g.copy_default(g.arg(CONIC_UBA), na_, "w+"+str(2*nx_+na_), "casadi_inf", false);
    g << "if (osqp_update_bounds(work, w, w+" + str(nx_+na_)+ ")) return 1;\n";

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
    g << "if (osqp_update_P_A(work, w, 0, " + str(nnzHupp_) + ", w+" + str(nnzHupp_) +
         ", 0, " + str(nnzA_) + ")) return 1;\n";

    g << "if (osqp_warm_start_x(work, " + g.arg(CONIC_X0) + ")) return 1;\n";
    g.copy_default(g.arg(CONIC_LAM_X0), nx_, "w", "0", false);
    g.copy_default(g.arg(CONIC_LAM_A0), na_, "w+"+str(nx_), "0", false);
    g << "if (osqp_warm_start_y(work, w)) return 1;\n";

    g << "if (osqp_solve(work)) return 1;\n";

    g.copy_check("&work->info->obj_val", 1, g.res(CONIC_COST), false, true);
    g.copy_check("work->solution->x", nx_, g.res(CONIC_X), false, true);
    g.copy_check("work->solution->y", nx_, g.res(CONIC_LAM_X), false, true);
    g.copy_check("work->solution->y+" + str(nx_), na_, g.res(CONIC_LAM_A), false, true);

    g << "if (work->info->status_val != OSQP_SOLVED) return 1;\n";
  }

  Dict OsqpInterface::get_stats(void* mem) const {
    Dict stats = Conic::get_stats(mem);
    auto m = static_cast<OsqpMemory*>(mem);
    stats["return_status"] = m->work->info->status;
    return stats;
  }

  OsqpMemory::OsqpMemory() {
  }

  OsqpMemory::~OsqpMemory() {
    osqp_cleanup(work);
  }

  OsqpInterface::OsqpInterface(DeserializingStream& s) : Conic(s) {
    s.version("OsqpInterface", 1);
    s.unpack("OsqpInterface::nnzHupp", nnzHupp_);
    s.unpack("OsqpInterface::nnzA", nnzA_);
    s.unpack("OsqpInterface::warm_start_primal", warm_start_primal_);
    s.unpack("OsqpInterface::warm_start_dual", warm_start_dual_);

    osqp_set_default_settings(&settings_);
    s.unpack("OsqpInterface::settings::rho", settings_.rho);
    s.unpack("OsqpInterface::settings::sigma", settings_.sigma);
    s.unpack("OsqpInterface::settings::scaling", settings_.scaling);
    s.unpack("OsqpInterface::settings::adaptive_rho", settings_.adaptive_rho);
    s.unpack("OsqpInterface::settings::adaptive_rho_interval", settings_.adaptive_rho_interval);
    s.unpack("OsqpInterface::settings::adaptive_rho_tolerance", settings_.adaptive_rho_tolerance);
    //s.unpack("OsqpInterface::settings::adaptive_rho_fraction", settings_.adaptive_rho_fraction);
    s.unpack("OsqpInterface::settings::max_iter", settings_.max_iter);
    s.unpack("OsqpInterface::settings::eps_abs", settings_.eps_abs);
    s.unpack("OsqpInterface::settings::eps_rel", settings_.eps_rel);
    s.unpack("OsqpInterface::settings::eps_prim_inf", settings_.eps_prim_inf);
    s.unpack("OsqpInterface::settings::eps_dual_inf", settings_.eps_dual_inf);
    s.unpack("OsqpInterface::settings::alpha", settings_.alpha);
    s.unpack("OsqpInterface::settings::delta", settings_.delta);
    s.unpack("OsqpInterface::settings::polish", settings_.polish);
    s.unpack("OsqpInterface::settings::polish_refine_iter", settings_.polish_refine_iter);
    s.unpack("OsqpInterface::settings::verbose", settings_.verbose);
    s.unpack("OsqpInterface::settings::scaled_termination", settings_.scaled_termination);
    s.unpack("OsqpInterface::settings::check_termination", settings_.check_termination);
    s.unpack("OsqpInterface::settings::warm_start", settings_.warm_start);
    //s.unpack("OsqpInterface::settings::time_limit", settings_.time_limit);
  }

  void OsqpInterface::serialize_body(SerializingStream &s) const {
    Conic::serialize_body(s);
    s.version("OsqpInterface", 1);
    s.pack("OsqpInterface::nnzHupp", nnzHupp_);
    s.pack("OsqpInterface::nnzA", nnzA_);
    s.pack("OsqpInterface::warm_start_primal", warm_start_primal_);
    s.pack("OsqpInterface::warm_start_dual", warm_start_dual_);
    s.pack("OsqpInterface::settings::rho", settings_.rho);
    s.pack("OsqpInterface::settings::sigma", settings_.sigma);
    s.pack("OsqpInterface::settings::scaling", settings_.scaling);
    s.pack("OsqpInterface::settings::adaptive_rho", settings_.adaptive_rho);
    s.pack("OsqpInterface::settings::adaptive_rho_interval", settings_.adaptive_rho_interval);
    s.pack("OsqpInterface::settings::adaptive_rho_tolerance", settings_.adaptive_rho_tolerance);
    //s.pack("OsqpInterface::settings::adaptive_rho_fraction", settings_.adaptive_rho_fraction);
    s.pack("OsqpInterface::settings::max_iter", settings_.max_iter);
    s.pack("OsqpInterface::settings::eps_abs", settings_.eps_abs);
    s.pack("OsqpInterface::settings::eps_rel", settings_.eps_rel);
    s.pack("OsqpInterface::settings::eps_prim_inf", settings_.eps_prim_inf);
    s.pack("OsqpInterface::settings::eps_dual_inf", settings_.eps_dual_inf);
    s.pack("OsqpInterface::settings::alpha", settings_.alpha);
    s.pack("OsqpInterface::settings::delta", settings_.delta);
    s.pack("OsqpInterface::settings::polish", settings_.polish);
    s.pack("OsqpInterface::settings::polish_refine_iter", settings_.polish_refine_iter);
    s.pack("OsqpInterface::settings::verbose", settings_.verbose);
    s.pack("OsqpInterface::settings::scaled_termination", settings_.scaled_termination);
    s.pack("OsqpInterface::settings::check_termination", settings_.check_termination);
    s.pack("OsqpInterface::settings::warm_start", settings_.warm_start);
    //s.pack("OsqpInterface::settings::time_limit", settings_.time_limit);
  }

} // namespace casadi
