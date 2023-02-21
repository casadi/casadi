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


#include "qrqp.hpp"
#include "casadi/core/nlpsol.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_CONIC_QRQP_EXPORT
  casadi_register_conic_qrqp(Conic::Plugin* plugin) {
    plugin->creator = Qrqp::creator;
    plugin->name = "qrqp";
    plugin->doc = Qrqp::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &Qrqp::options_;
    plugin->deserialize = &Qrqp::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_QRQP_EXPORT casadi_load_conic_qrqp() {
    Conic::registerPlugin(casadi_register_conic_qrqp);
  }

  Qrqp::Qrqp(const std::string& name, const std::map<std::string, Sparsity> &st)
    : Conic(name, st) {
  }

  Qrqp::~Qrqp() {
    clear_mem();
  }

  const Options Qrqp::options_
  = {{&Conic::options_},
     {{"max_iter",
       {OT_INT,
        "Maximum number of iterations [1000]."}},
      {"constr_viol_tol",
       {OT_DOUBLE,
        "Constraint violation tolerance [1e-8]."}},
      {"dual_inf_tol",
       {OT_DOUBLE,
        "Dual feasibility violation tolerance [1e-8]"}},
      {"print_header",
       {OT_BOOL,
        "Print header [true]."}},
      {"print_iter",
       {OT_BOOL,
        "Print iterations [true]."}},
      {"print_info",
       {OT_BOOL,
        "Print info [true]."}},
      {"print_lincomb",
       {OT_BOOL,
        "Print dependant linear combinations of constraints [false]. "
        "Printed numbers are 0-based indices into the vector of [simple bounds;linear bounds]"}},
      {"min_lam",
       {OT_DOUBLE,
        "Smallest multiplier treated as inactive for the initial active set [0]."}}
     }
  };

  void Qrqp::init(const Dict& opts) {
    // Initialize the base classes
    Conic::init(opts);

    // Transpose of the Jacobian
    AT_ = A_.T();

    // Assemble KKT system sparsity
    kkt_ = Sparsity::kkt(H_, A_, true, true);

    // Symbolic QR factorization
    kkt_.qr_sparse(sp_v_, sp_r_, prinv_, pc_);

    // Setup memory structure
    set_qp_prob();

    // Default options
    print_iter_ = true;
    print_header_ = true;
    print_info_ = true;
    print_lincomb_ = false;

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="max_iter") {
        p_.max_iter = op.second;
      } else if (op.first=="constr_viol_tol") {
        p_.constr_viol_tol = op.second;
      } else if (op.first=="dual_inf_tol") {
        p_.dual_inf_tol = op.second;
      } else if (op.first=="min_lam") {
        p_.min_lam = op.second;
      } else if (op.first=="print_iter") {
        print_iter_ = op.second;
      } else if (op.first=="print_header") {
        print_header_ = op.second;
      } else if (op.first=="print_info") {
        print_info_ = op.second;
      } else if (op.first=="print_lincomb") {
        print_lincomb_ = op.second;
      }
    }

    // Allocate memory
    casadi_int sz_w, sz_iw;
    casadi_qp_work(&p_, &sz_iw, &sz_w);
    alloc_iw(sz_iw, true);
    alloc_w(sz_w, true);

    if (print_header_) {
      // Print summary
      print("-------------------------------------------\n");
      print("This is casadi::QRQP\n");
      print("Number of variables:                       %9d\n", nx_);
      print("Number of constraints:                     %9d\n", na_);
      print("Number of nonzeros in H:                   %9d\n", H_.nnz());
      print("Number of nonzeros in A:                   %9d\n", A_.nnz());
      print("Number of nonzeros in KKT:                 %9d\n", kkt_.nnz());
      print("Number of nonzeros in QR(V):               %9d\n", sp_v_.nnz());
      print("Number of nonzeros in QR(R):               %9d\n", sp_r_.nnz());
    }
  }

  void Qrqp::set_qp_prob() {
    p_.sp_a = A_;
    p_.sp_h = H_;
    p_.sp_at = AT_;
    p_.sp_kkt = kkt_;
    p_.sp_v = sp_v_;
    p_.sp_r = sp_r_;
    p_.prinv = get_ptr(prinv_);
    p_.pc = get_ptr(pc_);
    casadi_qp_setup(&p_);
  }

  int Qrqp::init_mem(void* mem) const {
    if (Conic::init_mem(mem)) return 1;
    auto m = static_cast<QrqpMemory*>(mem);
    m->return_status = "";
    return 0;
  }



// SYMBOL "qp_singular_step"
// C-REPLACE "static_cast<T1*>(0)" "0"
template<typename T1>
int casadi_qp_singular_step2(casadi_qp_data<T1>* d) {
  // Local variables
  T1 tau_test, tau;
  casadi_int nnz_kkt, nk, k, i, best_k, best_neg, neg;
  const casadi_qp_prob<T1>* p = d->prob;
  // Find the columns that take part in any linear combination
  for (i = 0; i < p->nz; ++i) d->lincomb[i] = 0;
  for (k = 0; k < d->sing; ++k) {
    if (!d->has_search_dir) {
      casadi_qr_colcomb(d->dlam, d->nz_r, p->sp_r, p->pc, 1e-12, k);
    }
    for (i = 0; i < p->nz; ++i) if (fabs(d->dlam[i]) >= 1e-12) d->lincomb[i]++;
  }

  if (d->has_search_dir) {
    casadi_message("search dir");
    // One, given search direction
    nk = 1;
  } else {
    // QR factorization of the transpose
    casadi_trans(d->nz_kkt, p->sp_kkt, d->nz_v, p->sp_kkt, d->iw);
    nnz_kkt = p->sp_kkt[2+p->nz]; // kkt_colind[nz]
    casadi_copy(d->nz_v, nnz_kkt, d->nz_kkt);
    casadi_qr(p->sp_kkt, d->nz_kkt, d->w, p->sp_v, d->nz_v, p->sp_r, d->nz_r,
              d->beta, p->prinv, p->pc);
    // For all nullspace vectors
    nk = casadi_qr_singular(static_cast<T1*>(0), 0, d->nz_r, p->sp_r, p->pc, 1e-12);
    casadi_message("no search dir, nk = " + to_string(nk));
  }
  // Best flip
  best_k = best_neg = -1;
  tau = p->inf;
  for (k=0; k<nk; ++k) {
    if (!d->has_search_dir) {
      // Get a linear combination of the rows in kkt
      casadi_qr_colcomb(d->dz, d->nz_r, p->sp_r, p->pc, 1e-12, k);
    }
    // Which constraints can be flipped in order to increase rank?
    for (i=0; i<p->nz; ++i) {
      d->iw[i] = d->lincomb[i] && fabs(casadi_qp_kkt_dot(d, d->dz, i)) > 1e-12;
    }
    // Calculate step, dz and dlam
    casadi_qp_expand_step(d);


    if (d->ico >= 0) {
      // Try to reduce complementary slackness
      casadi_message("slackness: ico = " + str(d->ico) + ", co = " + str(d->co) + ", lbz = " + str(d->lbz[d->ico]) + 
      ", ubz = " + str(d->ubz[d->ico]) + ", z = " + str(d->z[d->ico])  + ", dz = " + str(d->dz[d->ico]) + 
      + ", lam = " + str(d->lam[d->ico])  + ", dlam = " + str(d->dlam[d->ico])
      )

      // Remove slackness violation by dropping constraint
      if (d->dlam[d->ico] != 0.) {
        tau = d->lam[d->ico] / d->dlam[d->ico];
        d->r_index = d->ico;
        d->r_sign = 0;
        best_k = k;

        casadi_message("dropping " + str(d->ico) + ", tau = " + str(tau));

        
      } else {
        casadi_message("not implemented");
        return 1;
      }

    }

    #if 0

    // Try both positive and negative direction
    for (neg = 0; neg < 2; ++neg) {
      // Negate direction
      if (neg) {
        casadi_scal(p->nz, -1., d->dz);
        casadi_scal(p->nz, -1., d->dlam);
        casadi_scal(p->nx, -1., d->tinfeas);
      }
      // Make sure primal infeasibility doesn't exceed limits
      if (casadi_qp_pr_direction(d)) continue;
      // Make sure dual infeasibility doesn't exceed limits
      if (casadi_qp_du_direction(d)) continue;
      // Loop over potential active set changes
      for (i=0; i<p->nz; ++i) {
        // Skip if no rank increase
        if (!d->iw[i]) continue;
        // Enforced or not?
        if (d->lam[i]==0.) {
          if (d->z[i] <= d->ubz[i] && (d->z[i] >= d->lbz[i] ?
              d->dz[i] < -1e-12 : d->dz[i] > 1e-12)) {
            // Enforce lower bound?
            if (!d->neverlower[i]
                && (tau_test = (d->lbz[i] - d->z[i]) / d->dz[i]) < tau
                && casadi_qp_enforceable(d, i, -1)) {
              tau = tau_test;
              d->r_index = i;
              d->r_sign = -1;
              best_k = k;
              best_neg = neg;
            }
          } else if (d->z[i] >= d->lbz[i] && (d->z[i] <= d->ubz[i] ?
              d->dz[i] > 1e-12 : d->dz[i] < -1e-12)) {
            // Enforce upper bound?
            if (!d->neverupper[i]
                && (tau_test = (d->ubz[i] - d->z[i]) / d->dz[i]) < tau
                && casadi_qp_enforceable(d, i, 1)) {
              tau = tau_test;
              d->r_index = i;
              d->r_sign = 1;
              best_k = k;
              best_neg = neg;
            }
          }
        } else if (!d->neverzero[i]) {
          // Drop a constraint?
          if (d->lam[i] > 0 ? d->dlam[i] < -1e-12 : d->dlam[i] > 1e-12) {
            if ((tau_test = -d->lam[i] / d->dlam[i]) < tau) {
              tau = tau_test;
              d->r_index = i;
              d->r_sign = 0;
              best_k = k;
              best_neg = neg;
            }
          }
        }
      }
    }
  #endif
  }


  // Can we restore feasibility?
  if (d->r_index < 0) return 1;
  // Recalculate direction, if needed
  if (--k != best_k) {
    // Need to recalculate direction
    casadi_qr_colcomb(d->dz, d->nz_r, p->sp_r, p->pc, 1e-12, best_k);
    casadi_qp_expand_step(d);
    //if (best_neg) tau *= -1;
  //} else if (--neg != best_neg) {
    // No need to recalculate, but opposite direction
    //tau *= -1;
  }
  // Scale step so that that tau=1 corresponds to a full step
  casadi_scal(p->nz, tau, d->dz);
  casadi_scal(p->nz, tau, d->dlam);
  casadi_scal(p->nx, tau, d->tinfeas);
  return 0;
}

// SYMBOL "qp_calc_step"
template<typename T1>
int casadi_qp_calc_step2(casadi_qp_data<T1>* d) {
  // Local variables
  const casadi_qp_prob<T1>* p = d->prob;
  // Reset returns
  d->r_index = -1;
  d->r_sign = 0;
  // Handle singularity
  if (d->sing) return casadi_qp_singular_step2(d);
  // Negative KKT residual
  casadi_qp_kkt_residual(d, d->dz);
  // Solve to get step in z[:nx] and lam[nx:]
  casadi_qr_solve(d->dz, 1, 1, p->sp_v, d->nz_v, p->sp_r, d->nz_r, d->beta,
                  p->prinv, p->pc, d->w);
  // Have step in dz[:nx] and dlam[nx:]. Calculate complete dz and dlam
  casadi_qp_expand_step(d);
  // Successful return
  return 0;
}

// SYMBOL "qp_linesearch"
template<typename T1>
void casadi_qp_linesearch2(casadi_qp_data<T1>* d) {
  // Local variables
  casadi_int du_index;
  // Start with a full step and no active set change
  d->sign = 0;
  d->index = -1;
  d->tau = 1.;


  if (d->sing && d->r_index >= 0) {
    //if (d->r_sign != 0 || casadi_qp_du_check(d, d->r_index)) {
      d->index = d->r_index;
      d->sign = d->r_sign;
      if (d->sign > 0) {
        d->msg = "Enforced ubz for regularity";
      } else if (d->sign < 0) {
        d->msg = "Enforced lbz for regularity";
      } else if (d->lam[d->index] > 0) {
        d->msg = "Dropped ubz for regularity";
      } else {
        d->msg = "Dropped lbz for regularity";
      }
      d->msg_ind = d->index;
//    }
  }


  // Find largest possible step without exceeding acceptable |pr|
  casadi_qp_primal_blocking(d);
  // Find largest possible step without exceeding acceptable |du|
  du_index = casadi_qp_dual_blocking(d);
  // Take primal-dual step, avoiding accidental sign changes for lam
  casadi_qp_take_step(d);
  // Handle dual blocking constraints
  if (du_index >= 0) {
    // Sensititivity in decreasing du_index
    casadi_qp_calc_sens(d, du_index);
    // Find corresponding index
    casadi_qp_du_index(d);
  }
}

// SYMBOL "qp_iterate"
template<typename T1>
int casadi_qp_iterate2(casadi_qp_data<T1>* d) {
  // Reset message flag
  d->msg = 0;
  // Start a new iteration
  d->iter++;
  // Calculate search direction
  if (casadi_qp_calc_step2(d)) {
    d->status = QP_NO_SEARCH_DIR;
    return 1;
  }
  // Line search in the calculated direction
  casadi_qp_linesearch2(d);
  // Keep iterating
  return 0;
}



  int Qrqp::
  solve(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<QrqpMemory*>(mem);
    // Message buffer
    char buf[121];
    // Setup data structure
    casadi_qp_data<double> d;
    d.prob = &p_;
    d.nz_h = arg[CONIC_H];
    d.g = arg[CONIC_G];
    d.nz_a = arg[CONIC_A];
    casadi_qp_init(&d, &iw, &w);
    // Pass bounds on z
    casadi_copy(arg[CONIC_LBX], nx_, d.lbz);
    casadi_copy(arg[CONIC_LBA], na_, d.lbz+nx_);
    casadi_copy(arg[CONIC_UBX], nx_, d.ubz);
    casadi_copy(arg[CONIC_UBA], na_, d.ubz+nx_);
    // Pass initial guess
    casadi_copy(arg[CONIC_X0], nx_, d.z);
    casadi_fill(d.z+nx_, na_, nan);
    casadi_copy(arg[CONIC_LAM_X0], nx_, d.lam);
    casadi_copy(arg[CONIC_LAM_A0], na_, d.lam+nx_);



    uout() << std::vector<double>(d.z, d.z + nx_ + na_) << "\n";
    uout() << std::vector<double>(d.lam, d.lam + nx_ + na_) << "\n";

    // Reset solver
    if (casadi_qp_reset(&d)) return 1;
    while (true) {
      // Prepare QP
      int flag = casadi_qp_prepare(&d);
      // Print iteration progress
      if (print_iter_) {
        if (d.iter % 10 == 0) {
          // Print header
          if (casadi_qp_print_header(&d, buf, sizeof(buf))) break;
          uout() << buf << "\n";
        }
        // Print iteration
        if (casadi_qp_print_iteration(&d, buf, sizeof(buf))) break;
        uout() << buf << "\n";
      }
      uout() << "z = " << std::vector<double>(d.z, d.z + nx_ + na_) << "\n";
      uout() << "lam = " << std::vector<double>(d.lam, d.lam + nx_ + na_) << "\n";
      uout() << "dz = " << std::vector<double>(d.dz, d.dz + nx_ + na_) << "\n";
      uout() << "dlam = " << std::vector<double>(d.dlam, d.dlam + nx_ + na_) << "\n";
      uout() << "lbz = " << std::vector<double>(d.lbz, d.lbz + nx_ + na_) << "\n";
      uout() << "ubz = " << std::vector<double>(d.ubz, d.ubz + nx_ + na_) << "\n";
      // Make an iteration
      flag = flag || casadi_qp_iterate2(&d);
      // Print debug info
      if (print_lincomb_) {
        for (casadi_int k=0;k<d.sing;++k) {
          uout() << "lincomb: ";
          casadi_qp_print_colcomb(&d, buf, sizeof(buf), k);
          uout() << buf << "\n";
        }
      }
      if (flag) break;

      // User interrupt
      InterruptHandler::check();
    }
    // Check return flag
    switch (d.status) {
      case QP_SUCCESS:
        m->return_status = "success";
        break;
      case QP_MAX_ITER:
        m->return_status = "Maximum number of iterations reached";
        m->unified_return_status = SOLVER_RET_LIMITED;
        break;
      case QP_NO_SEARCH_DIR:
        m->return_status = "Failed to calculate search direction";
        break;
      case QP_PRINTING_ERROR:
        m->return_status = "Printing error";
        break;
    }
    // Get solution
    casadi_copy(&d.f, 1, res[CONIC_COST]);
    casadi_copy(d.z, nx_, res[CONIC_X]);
    casadi_copy(d.lam, nx_, res[CONIC_LAM_X]);
    casadi_copy(d.lam+nx_, na_, res[CONIC_LAM_A]);
    // Return
    if (verbose_) casadi_warning(m->return_status);
    m->success = d.status == QP_SUCCESS;
    return 0;
  }

  void Qrqp::codegen_body(CodeGenerator& g) const {
    g.add_auxiliary(CodeGenerator::AUX_QP);
    if (print_iter_) g.add_auxiliary(CodeGenerator::AUX_PRINTF);
    g.local("d", "struct casadi_qp_data");
    g.local("p", "struct casadi_qp_prob");
    g.local("flag", "int");
    if (print_iter_ || print_lincomb_) g.local("buf[121]", "char");
    if (print_lincomb_) g.local("k", "casadi_int");

    // Setup memory structure
    g << "p.sp_a = " << g.sparsity(A_) << ";\n";
    g << "p.sp_h = " << g.sparsity(H_) << ";\n";
    g << "p.sp_at = " << g.sparsity(AT_) << ";\n";
    g << "p.sp_kkt = " << g.sparsity(kkt_) << ";\n";
    g << "p.sp_v = " << g.sparsity(sp_v_) << ";\n";
    g << "p.sp_r = " << g.sparsity(sp_r_) << ";\n";
    g << "p.prinv = " << g.constant(prinv_) << ";\n";
    g << "p.pc =  " << g.constant(pc_) << ";\n";
    g << "casadi_qp_setup(&p);\n";

    // Copy options
    g << "p.max_iter = " << p_.max_iter << ";\n";
    g << "p.min_lam = " << p_.min_lam << ";\n";
    g << "p.constr_viol_tol = " << p_.constr_viol_tol << ";\n";
    g << "p.dual_inf_tol = " << p_.dual_inf_tol << ";\n";

    // Setup data structure
    g << "d.prob = &p;\n";
    g << "d.nz_h = arg[" << CONIC_H << "];\n";
    g << "d.g = arg[" << CONIC_G << "];\n";
    g << "d.nz_a = arg[" << CONIC_A << "];\n";
    g << "casadi_qp_init(&d, &iw, &w);\n";

    g.comment("Pass bounds on z");
    g.copy_default(g.arg(CONIC_LBX), nx_, "d.lbz", "-casadi_inf", false);
    g.copy_default(g.arg(CONIC_LBA), na_, "d.lbz+" + str(nx_), "-casadi_inf", false);
    g.copy_default(g.arg(CONIC_UBX), nx_, "d.ubz", "casadi_inf", false);
    g.copy_default(g.arg(CONIC_UBA), na_, "d.ubz+" + str(nx_), "casadi_inf", false);

    g.comment("Pass initial guess");
    g.copy_default(g.arg(CONIC_X0), nx_, "d.z", "0", false);
    g << g.fill("d.z+"+str(nx_), na_, g.constant(nan)) << "\n";
    g.copy_default(g.arg(CONIC_LAM_X0), nx_, "d.lam", "0", false);
    g.copy_default(g.arg(CONIC_LAM_A0), na_, "d.lam+" + str(nx_), "0", false);

    g.comment("Solve QP");
    g << "if (casadi_qp_reset(&d)) return 1;\n";
    g << "while (1) {\n";
    g << "flag = casadi_qp_prepare(&d);\n";
    if (print_iter_) {
      // Print header
      g << "if (d.iter % 10 == 0) {\n";
      g << "if (casadi_qp_print_header(&d, buf, sizeof(buf))) break;\n";
      g << g.printf("%s\\n", "buf") << "\n";
      g << "}\n";
      // Print iteration
      g << "if (casadi_qp_print_iteration(&d, buf, sizeof(buf))) break;\n";
      g << g.printf("%s\\n", "buf") << "\n";
    }
    g << "if (flag || casadi_qp_iterate(&d)) break;\n";
    if (print_lincomb_) {
      g << "for (k=0;k<d.sing;++k) {\n";
      g << "casadi_qp_print_colcomb(&d, buf, sizeof(buf), k);\n";
      g << g.printf("lincomb: %s\\n", "buf") << "\n";
      g << "}\n";
    }
    g << "}\n";

    g.comment("Get solution");
    g.copy_check("&d.f", 1, g.res(CONIC_COST), false, true);
    g.copy_check("d.z", nx_, g.res(CONIC_X), false, true);
    g.copy_check("d.lam", nx_, g.res(CONIC_LAM_X), false, true);
    g.copy_check("d.lam+"+str(nx_), na_, g.res(CONIC_LAM_A), false, true);

    g << "return d.status != QP_SUCCESS;\n";
  }

  Dict Qrqp::get_stats(void* mem) const {
    Dict stats = Conic::get_stats(mem);
    auto m = static_cast<QrqpMemory*>(mem);
    stats["return_status"] = m->return_status;
    return stats;
  }

  Qrqp::Qrqp(DeserializingStream& s) : Conic(s) {
    s.version("Qrqp", 1);
    s.unpack("Qrqp::AT", AT_);
    s.unpack("Qrqp::kkt", kkt_);
    s.unpack("Qrqp::sp_v", sp_v_);
    s.unpack("Qrqp::sp_r", sp_r_);
    s.unpack("Qrqp::prinv", prinv_);
    s.unpack("Qrqp::pc", pc_);
    s.unpack("Qrqp::print_iter", print_iter_);
    s.unpack("Qrqp::print_header", print_header_);
    s.unpack("Qrqp::print_info", print_info_);
    s.unpack("Qrqp::print_lincomb_", print_lincomb_);
    set_qp_prob();
    s.unpack("Qrqp::max_iter", p_.max_iter);
    s.unpack("Qrqp::min_lam", p_.min_lam);
    s.unpack("Qrqp::constr_viol_tol", p_.constr_viol_tol);
    s.unpack("Qrqp::dual_inf_tol", p_.dual_inf_tol);
  }

  void Qrqp::serialize_body(SerializingStream &s) const {
    Conic::serialize_body(s);

    s.version("Qrqp", 1);
    s.pack("Qrqp::AT", AT_);
    s.pack("Qrqp::kkt", kkt_);
    s.pack("Qrqp::sp_v", sp_v_);
    s.pack("Qrqp::sp_r", sp_r_);
    s.pack("Qrqp::prinv", prinv_);
    s.pack("Qrqp::pc", pc_);
    s.pack("Qrqp::print_iter", print_iter_);
    s.pack("Qrqp::print_header", print_header_);
    s.pack("Qrqp::print_info", print_info_);
    s.pack("Qrqp::print_lincomb_", print_lincomb_);
    s.pack("Qrqp::max_iter", p_.max_iter);
    s.pack("Qrqp::min_lam", p_.min_lam);
    s.pack("Qrqp::constr_viol_tol", p_.constr_viol_tol);
    s.pack("Qrqp::dual_inf_tol", p_.dual_inf_tol);
  }

} // namespace casadi
