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
        "Use lam_a0 and lam_x0 input to warmstart [Default: true]."}}
     }
  };

  void QpalmInterface::init(const Dict& opts) {
    // Initialize the base classes
    Conic::init(opts);

    warm_start_dual_ = true;
    warm_start_primal_ = true;

    Sparsity Asp = vertcat(Sparsity::diag(nx_), A_);

    qpalm_set_default_settings(&settings_);

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
          } else if (op.first=="tau_init") {
            settings_.tau_init = op.second;
          } else if (op.first=="proximal") {
            settings_.proximal = op.second;
          } else if (op.first=="gamma") {
            settings_.gamma = op.second;
          } else if (op.first=="gamma_upd") {
            settings_.gamma_upd = op.second;
          } else if (op.first=="gamma_max") {
            settings_.gamma_max = op.second;
          } else if (op.first=="scaling") {
            settings_.scaling = op.second;
          } else if (op.first=="verbose") {
            settings_.verbose = op.second;
          } else if (op.first=="warm_start") {
            settings_.warm_start = true;
            casadi_error("Use CasADi options 'warm_start_primal' and 'warm_start_dual' instead.");
          } else {
            casadi_error("Not recognised");
          }
        }
      }
    }

    alloc_w(na_+nx_, true); // lb
    alloc_w(na_+nx_, true); // ub
    alloc_w(nx_, true); // q
    alloc_w(na_+nx_, true); // y0
    alloc_w(nx_+A_.nnz(), true); // Adata
    alloc_w(H_.nnz(), true); // Qdata
  }

  int QpalmInterface::init_mem(void* mem) const {
    if (Conic::init_mem(mem)) return 1;
    auto m = static_cast<QpalmMemory*>(mem);

    Sparsity Asp = vertcat(Sparsity::diag(nx_), A_);
    
    m->A_row = to_int(Asp.get_row());
    m->A_colind = to_int(Asp.get_colind());
    if (m->A_row.empty()) m->A_row.push_back(0);
    if (m->A_colind.empty()) m->A_colind.push_back(0);

    m->A.nrow = nx_ + na_;
    m->A.ncol = nx_;
    m->A.nzmax = A_.nnz()+nx_;
    m->A.nz = nullptr;
    m->A.itype = CHOLMOD_INT;
    m->A.stype = 0;
    m->A.xtype = CHOLMOD_REAL;
    m->A.dtype = CHOLMOD_DOUBLE;
    m->A.sorted = 1;
    m->A.packed = 1;
    m->Q.z = nullptr;
    m->A.i = static_cast<void*>(get_ptr(m->A_row));
    m->A.p = static_cast<void*>(get_ptr(m->A_colind));

    m->Q_row = to_int(H_.get_row());
    m->Q_colind = to_int(H_.get_colind());
    if (m->Q_row.empty()) m->Q_row.push_back(0);
    if (m->Q_colind.empty()) m->Q_colind.push_back(0);

    m->Q.nrow = nx_;
    m->Q.ncol = nx_;
    m->Q.nzmax = H_.nnz();
    m->Q.nz = nullptr;
    m->Q.itype = CHOLMOD_INT;
    m->Q.stype = -1; // Lower triangular symmetric
    m->Q.xtype = CHOLMOD_REAL;
    m->Q.dtype = CHOLMOD_DOUBLE;
    m->Q.sorted = 1;
    m->Q.packed = 1;
    m->Q.z = nullptr;
    m->Q.i = static_cast<void*>(get_ptr(m->Q_row));
    m->Q.p = static_cast<void*>(get_ptr(m->Q_colind));

    m->data.n = nx_;
    m->data.m = na_+nx_;
    m->data.A = &m->A;
    m->data.Q = &m->Q;

    // Setup workspace

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

    m->fstats.at("preprocessing").tic();
    // Inputs
    const double *a=arg[CONIC_A];

    double* lb = w; w+= nx_+na_;
    double* ub = w; w+= nx_+na_;
    double* q = w; w+= nx_;
    double* y0 = w; w+= nx_+na_;
    double* Adata = w; w+= nx_+A_.nnz();
    double* Qdata = w; w+= H_.nnz();

    m->data.A->x = Adata;
    m->data.Q->x = Qdata;

    m->data.bmin = lb;
    m->data.bmax = ub;
    m->data.q = q;

    // Copy bounds
    casadi_copy(arg[CONIC_LBX], nx_, lb);
    casadi_copy(arg[CONIC_LBA], na_, lb+nx_);
    casadi_copy(arg[CONIC_UBX], nx_, ub);
    casadi_copy(arg[CONIC_UBA], na_, ub+nx_);

    for (casadi_int i=0;i<nx_+na_;++i) {
      lb[i] = fmax(lb[i], -QPALM_INFTY);
      ub[i] = fmin(ub[i], QPALM_INFTY);
    }

    // Get linear part of objective
    casadi_copy(arg[CONIC_G], nx_, q);

    // Project Hessian
    casadi_copy(arg[CONIC_H], H_.nnz(), Qdata);

    // Get contraint matrix
    const casadi_int* colind = A_.colind();
    casadi_int offset = 0;
    // Loop over columns
    for (casadi_int i=0; i<nx_; ++i) {
      Adata[offset] = 1;
      offset++;
      casadi_int n = colind[i+1]-colind[i];
      casadi_copy(a+colind[i], n, Adata+offset);
      offset+= n;
    }

    m->work = qpalm_setup(&m->data, &settings_, &m->c);

    if (warm_start_dual_) {
      casadi_copy(arg[CONIC_LAM_X0], nx_, y0);
      casadi_copy(arg[CONIC_LAM_A0], na_, y0+nx_);
    }
    
    qpalm_warm_start(m->work, warm_start_primal_ ? arg[CONIC_X0] : nullptr, warm_start_dual_ ? y0 : nullptr);

    m->fstats.at("preprocessing").toc();
    m->fstats.at("solver").tic();
  
    // Solve Problem
    qpalm_solve(m->work);
    m->fstats.at("solver").toc();
    m->fstats.at("postprocessing").tic();
    printf("Solver status: ");
    printf("%s", m->work->info->status);
    printf(" \n");
    printf("Iter: %d\n", (int) m->work->info->iter);
    printf("Iter Out: %d\n", (int) m->work->info->iter_out);

    m->success = m->work->info->status_val == QPALM_SOLVED;
    if (m->success) m->unified_return_status = SOLVER_RET_SUCCESS;

    casadi_copy(m->work->solution->x, nx_, res[CONIC_X]);
    casadi_copy(m->work->solution->y, nx_, res[CONIC_LAM_X]);
    casadi_copy(m->work->solution->y+nx_, na_, res[CONIC_LAM_A]);

    if (res[CONIC_COST])
      *res[CONIC_COST] = casadi_dot(nx_, m->work->solution->x, arg[CONIC_G]);
    if (res[CONIC_COST])
      *res[CONIC_COST] += 0.5*casadi_bilin(arg[CONIC_H], H_, m->work->solution->x, m->work->solution->x);

    m->fstats.at("postprocessing").toc();

    return 0;
  }

  void QpalmInterface::codegen_free_mem(CodeGenerator& g) const {
 
  }

  void QpalmInterface::codegen_declarations(CodeGenerator& g) const {
    g.auxiliaries << "struct QpalmMemoryStruct {\n";
    g.auxiliaries << "QPALMData data;\n";
    g.auxiliaries << "QPALMWorkspace *work;\n";
    g.auxiliaries << "cholmod_common c;\n";
    g.auxiliaries << "cholmod_sparse A, Q;\n";
    g.auxiliaries << "QPALMSettings settings;\n";
    g.auxiliaries << "int A_colind[];\n";
    g.auxiliaries << "int A_row[];\n";
    g.auxiliaries << "int Q_colind[];\n";
    g.auxiliaries << "int Q_row[];\n";
    g.auxiliaries << "};\n";

    g.auxiliaries << "typedef struct QpalmMemoryStruct QpalmMemory;\n";
  }

  void QpalmInterface::codegen_init_mem(CodeGenerator& g) const {

    std::string m = codegen_mem(g);

    g << m << "->A.nrow = " << nx_ + na_ << ";\n";
    g << m << "->A.ncol = " << nx_ << ";\n";
    g << m << "->A.nzmax = " << A_.nnz()+nx_ << ";\n";
    g << m << "->A.nz = 0;\n";
    g << m << "->A.itype = CHOLMOD_INT;\n";
    g << m << "->A.stype = 0;\n";
    g << m << "->A.xtype = CHOLMOD_REAL;\n";
    g << m << "->A.dtype = CHOLMOD_DOUBLE;\n";
    g << m << "->A.sorted = 1;\n";
    g << m << "->A.packed = 1;\n";
    g << m << "->A.z = 0;\n";
    //g << m << "->A.i = //;\n";
    //g << m << "->A.p = //;\n";

    g << m << "->Q.nrow = " << nx_ << ";\n";
    g << m << "->Q.ncol = " << nx_ << ";\n";
    g << m << "->Q.nzmax = " << H_.nnz() << ";\n";
    g << m << "->Q.nz = 0;\n";
    g << m << "->Q.itype = CHOLMOD_INT;\n";
    g << m << "->Q.stype = -1;\n";
    g << m << "->Q.xtype = CHOLMOD_REAL;\n";
    g << m << "->Q.dtype = CHOLMOD_DOUBLE;\n";
    g << m << "->Q.sorted = 1;\n";
    g << m << "->Q.packed = 1;\n";
    g << m << "->Q.z = 0;\n";
    //g << m << "->Q.i = //;\n";
    //g << m << "->Q.p = //;\n";

    g << m << "->data.n = " << nx_ << ";\n";
    g << m << "->data.m = " << na_ + nx_ << ";\n";
    g << m << "->data.A = &" << m << "->A;\n";
    g << m << "->data.Q = &" << m << "->Q;\n";

    g << "qpalm_set_default_settings(&" << m << "->settings" << ");\n";
    g << m << "->settings.max_iter = " << settings_.max_iter << ";\n";
    g << m << "->settings.eps_abs = " << settings_.eps_abs << ";\n";
    g << m << "->settings.eps_rel = " << settings_.eps_rel << ";\n";
    g << m << "->settings.eps_abs_in = " << settings_.eps_abs_in << ";\n";
    g << m << "->settings.eps_rel_in = " << settings_.eps_rel_in << ";\n";
    g << m << "->settings.rho = " << settings_.rho << ";\n";
    g << m << "->settings.eps_prim_inf = " << settings_.eps_prim_inf << ";\n";
    g << m << "->settings.eps_dual_inf = " << settings_.eps_dual_inf << ";\n";
    g << m << "->settings.theta = " << settings_.theta << ";\n";
    g << m << "->settings.delta = " << settings_.delta << ";\n";
    g << m << "->settings.tau_init = " << settings_.tau_init << ";\n";
    g << m << "->settings.proximal = " << settings_.proximal << ";\n";
    g << m << "->settings.gamma = " << settings_.gamma << ";\n";
    g << m << "->settings.gamma_upd = " << settings_.gamma_upd << ";\n";
    g << m << "->settings.gamma_max = " << settings_.gamma_max << ";\n";
    g << m << "->settings.scaling = " << settings_.scaling << ";\n";
    g << m << "->settings.verbose = " << settings_.verbose << ";\n";
    g << m << "->settings.warm_start = " << settings_.warm_start << ";\n";

    g << "return 0;\n";
  }

  void QpalmInterface::codegen_body(CodeGenerator& g) const {
    g.add_include("qpalm/qpalm.h");
    g.add_include("qpalm/constants.h");
  }

  Dict QpalmInterface::get_stats(void* mem) const {
    Dict stats = Conic::get_stats(mem);
    auto m = static_cast<QpalmMemory*>(mem);
    stats["return_status"] = m->work->info->status;
    return stats;
  }

  QpalmMemory::QpalmMemory() : work(nullptr) {
  }

  QpalmMemory::~QpalmMemory() {
    if (work) qpalm_cleanup(work);
  }

  QpalmInterface::QpalmInterface(DeserializingStream& s) : Conic(s) {
    s.version("QpalmInterface", 1);
    s.unpack("QpalmInterface::warm_start_primal", warm_start_primal_);
    s.unpack("QpalmInterface::warm_start_dual", warm_start_dual_);
    s.unpack("QpalmInterface::settings::max_iter", settings_.max_iter);
    s.unpack("QpalmInterface::settings::eps_abs", settings_.eps_abs);
    s.unpack("QpalmInterface::settings::eps_rel", settings_.eps_rel);
    s.unpack("QpalmInterface::settings::eps_abs_in", settings_.eps_abs_in);
    s.unpack("QpalmInterface::settings::eps_rel_in", settings_.eps_rel_in);
    s.unpack("QpalmInterface::settings::rho", settings_.rho);
    s.unpack("QpalmInterface::settings::eps_prim_inf", settings_.eps_prim_inf);
    s.unpack("QpalmInterface::settings::eps_dual_inf", settings_.eps_dual_inf);
    s.unpack("QpalmInterface::settings::theta", settings_.theta);
    s.unpack("QpalmInterface::settings::delta", settings_.delta);
    s.unpack("QpalmInterface::settings::tau_init", settings_.tau_init);
    s.unpack("QpalmInterface::settings::proximal", settings_.proximal);
    s.unpack("QpalmInterface::settings::gamma", settings_.gamma);
    s.unpack("QpalmInterface::settings::gamma_upd", settings_.gamma_upd);
    s.unpack("QpalmInterface::settings::gamma_max", settings_.gamma_max);
    s.unpack("QpalmInterface::settings::scaling", settings_.scaling);
    s.unpack("QpalmInterface::settings::verbose", settings_.verbose);
    s.unpack("QpalmInterface::settings::warm_start", settings_.warm_start);
  }

  void QpalmInterface::serialize_body(SerializingStream &s) const {
    Conic::serialize_body(s);
    s.version("QpalmInterface", 1);
    s.pack("QpalmInterface::warm_start_primal", warm_start_primal_);
    s.pack("QpalmInterface::warm_start_dual", warm_start_dual_);
    s.pack("QpalmInterface::settings::max_iter", settings_.max_iter);
    s.pack("QpalmInterface::settings::eps_abs", settings_.eps_abs);
    s.pack("QpalmInterface::settings::eps_rel", settings_.eps_rel);
    s.pack("QpalmInterface::settings::eps_abs_in", settings_.eps_abs_in);
    s.pack("QpalmInterface::settings::eps_rel_in", settings_.eps_rel_in);
    s.pack("QpalmInterface::settings::rho", settings_.rho);
    s.pack("QpalmInterface::settings::eps_prim_inf", settings_.eps_prim_inf);
    s.pack("QpalmInterface::settings::eps_dual_inf", settings_.eps_dual_inf);
    s.pack("QpalmInterface::settings::theta", settings_.theta);
    s.pack("QpalmInterface::settings::delta", settings_.delta);
    s.pack("QpalmInterface::settings::tau_init", settings_.tau_init);
    s.pack("QpalmInterface::settings::proximal", settings_.proximal);
    s.pack("QpalmInterface::settings::gamma", settings_.gamma);
    s.pack("QpalmInterface::settings::gamma_upd", settings_.gamma_upd);
    s.pack("QpalmInterface::settings::gamma_max", settings_.gamma_max);
    s.pack("QpalmInterface::settings::scaling", settings_.scaling);
    s.pack("QpalmInterface::settings::verbose", settings_.verbose);
    s.pack("QpalmInterface::settings::warm_start", settings_.warm_start);
  }

} // namespace casadi
