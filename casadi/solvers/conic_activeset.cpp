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


#include "conic_activeset.hpp"
#include "casadi/core/nlpsol.hpp"

using namespace std;
namespace casadi {

  template<typename T1>
  void casadi_qp_work(casadi_qp_prob<T1>* p, casadi_int* sz_w, casadi_int* sz_iw) {
    // Local variables
    casadi_int nnz_a, nnz_kkt, nnz_v, nnz_r;
    // Get matrix number of nonzeros
    nnz_a = p->sp_a[2+p->sp_a[1]];
    nnz_kkt = p->sp_kkt[2+p->sp_kkt[1]];
    nnz_v = p->sp_v[2+p->sp_v[1]];
    nnz_r = p->sp_r[2+p->sp_r[1]];
    // Reset sz_w, sz_iw
    *sz_w = *sz_iw = 0;
    // Temporary work vectors
    *sz_w = max(*sz_w, p->nz); // casadi_project, tau memory
    *sz_iw = max(*sz_iw, p->nz); // casadi_trans, tau type, allzero
    *sz_w = max(*sz_w, 2*p->nz); // casadi_qr
    // Persistent work vectors
    *sz_w += nnz_kkt; // kkt
    *sz_w += p->nz; // z=[xk,gk]
    *sz_w += p->nz; // lbz
    *sz_w += p->nz; // ubz
    *sz_w += p->nz; // lam
    *sz_w += nnz_a; // trans(a)
    *sz_w += p->nz; // dz
    *sz_w += p->nz; // dlam
    *sz_w += p->nx; // infeas
    *sz_w += p->nx; // tinfeas
    *sz_iw += p->nz; // neverzero
    *sz_iw += p->nz; // neverupper
    *sz_iw += p->nz; // neverlower
    *sz_w += max(nnz_v + nnz_r, nnz_kkt); // either v & r or trans(kkt)
    *sz_w += p->nz; // beta
  }

  extern "C"
  int CASADI_CONIC_ACTIVESET_EXPORT
  casadi_register_conic_activeset(Conic::Plugin* plugin) {
    plugin->creator = ConicActiveSet::creator;
    plugin->name = "activeset";
    plugin->doc = ConicActiveSet::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &ConicActiveSet::options_;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_ACTIVESET_EXPORT casadi_load_conic_activeset() {
    Conic::registerPlugin(casadi_register_conic_activeset);
  }

  ConicActiveSet::ConicActiveSet(const std::string& name, const std::map<std::string, Sparsity> &st)
    : Conic(name, st) {
  }

  ConicActiveSet::~ConicActiveSet() {
    clear_mem();
  }

  Options ConicActiveSet::options_
  = {{&Conic::options_},
     {{"max_iter",
       {OT_INT,
        "Maximum number of iterations [1000]."}},
      {"tol",
       {OT_DOUBLE,
        "Tolerance [1e-8]."}},
      {"du_to_pr",
       {OT_DOUBLE,
        "How much larger dual than primal error is acceptable [1000]"}},
      {"print_header",
       {OT_BOOL,
        "Print header [true]."}},
      {"print_iter",
       {OT_BOOL,
        "Print iterations [true]."}}
     }
  };

  void ConicActiveSet::init(const Dict& opts) {
    // Initialize the base classes
    Conic::init(opts);

    // Default options
    max_iter_ = 1000;
    tol_ = 1e-8;
    print_iter_ = true;
    print_header_ = true;
    du_to_pr_ = 1000.;

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="max_iter") {
        max_iter_ = op.second;
      } else if (op.first=="tol") {
        tol_ = op.second;
      } else if (op.first=="print_iter") {
        print_iter_ = op.second;
      } else if (op.first=="print_header") {
        print_header_ = op.second;
      } else if (op.first=="du_to_pr") {
        du_to_pr_ = op.second;
      }
    }

    // Transpose of the Jacobian
    AT_ = A_.T();

    // Assemble KKT system sparsity
    kkt_ = Sparsity::kkt(H_, A_, true, true);

    // Symbolic QR factorization
    kkt_.qr_sparse(sp_v_, sp_r_, prinv_, pc_);

    // Setup memory structure
    p_.du_to_pr = du_to_pr_;
    p_.print_iter = print_iter_;
    p_.sp_a = A_;
    p_.sp_h = H_;
    p_.sp_at = AT_;
    p_.sp_kkt = kkt_;
    p_.sp_v = sp_v_;
    p_.sp_r = sp_r_;
    p_.prinv = get_ptr(prinv_);
    p_.pc = get_ptr(pc_);
    p_.dmin = std::numeric_limits<double>::min();
    p_.inf = inf;
    p_.nx = nx_;
    p_.na = na_;
    p_.nz = nx_+na_;

    // Allocate memory
    casadi_int sz_w, sz_iw;
    casadi_qp_work(&p_, &sz_w, &sz_iw);
    alloc_w(sz_w, true);
    alloc_iw(sz_iw, true);

    if (print_header_) {
      // Print summary
      print("-------------------------------------------\n");
      print("This is casadi::ConicActiveSet.\n");
      print("Number of variables:                       %9d\n", nx_);
      print("Number of constraints:                     %9d\n", na_);
      print("Work in progress!\n");
    }
  }

  template<typename T1>
  void casadi_qp_pr(casadi_qp_data<T1>* d) {
    // Calculate largest constraint violation
    casadi_int i;
    const casadi_qp_prob<T1>* p = d->prob;
    d->pr = 0;
    d->ipr = -1;
    for (i=0; i<p->nz; ++i) {
      if (d->z[i] > d->ubz[i]+d->pr) {
        d->pr = d->z[i]-d->ubz[i];
        d->ipr = i;
      } else if (d->z[i] < d->lbz[i]-d->pr) {
        d->pr = d->lbz[i]-d->z[i];
        d->ipr = i;
      }
    }
  }

  template<typename T1>
  void casadi_qp_du(casadi_qp_data<T1>* d) {
    // Calculate largest constraint violation
    casadi_int i;
    const casadi_qp_prob<T1>* p = d->prob;
    d->du = 0;
    d->idu = -1;
    for (i=0; i<p->nx; ++i) {
      if (d->infeas[i] > d->du) {
        d->du = d->infeas[i];
        d->idu = i;
      } else if (d->infeas[i] < -d->du) {
        d->du = -d->infeas[i];
        d->idu = i;
      }
    }
  }

  template<typename T1>
  void casadi_qp_log(casadi_qp_data<T1>* d, const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    vsnprintf(d->msg, sizeof(d->msg), fmt, args);
    va_end(args);
  }

  template<typename T1>
  casadi_int casadi_qp_pr_index(casadi_qp_data<T1>* d, casadi_int* sign) {
    // Try to improve primal feasibility by adding a constraint
    if (d->lam[d->ipr]==0.) {
      // Add the most violating constraint
      *sign = d->z[d->ipr]<d->lbz[d->ipr] ? -1 : 1;
      casadi_qp_log(d, "Added %lld to reduce |pr|", d->ipr);
      return d->ipr;
    }
    return -1;
  }

  template<typename T1>
  T1 casadi_qp_du_check(casadi_qp_data<T1>* d, casadi_int i) {
    // Local variables
    casadi_int k;
    T1 new_du;
    const casadi_int *at_colind, *at_row;
    const casadi_qp_prob<T1>* p = d->prob;
    // AT sparsity
    at_colind = p->sp_at + 2;
    at_row = at_colind + p->na + 1;
    // Maximum infeasibility from setting from setting lam[i]=0
    if (i<p->nx) {
      new_du = fabs(d->infeas[i]-d->lam[i]);
    } else {
      new_du = 0.;
      for (k=at_colind[i-p->nx]; k<at_colind[i-p->nx+1]; ++k) {
        new_du = fmax(new_du, fabs(d->infeas[at_row[k]]-d->nz_at[k]*d->lam[i]));
      }
    }
    return new_du;
  }

  template<typename T1>
  casadi_int casadi_qp_du_index(casadi_qp_data<T1>* d, casadi_int* sign) {
    // Try to improve dual feasibility by removing a constraint
    // Local variables
    casadi_int best_ind, i;
    T1 best_w;
    const casadi_qp_prob<T1>* p = d->prob;
    // We need to increase or decrease infeas[idu]. Sensitivity:
    casadi_fill(d->w, p->nz, 0.);
    d->w[d->idu] = d->infeas[d->idu]>0 ? -1. : 1.;
    casadi_mv(d->nz_a, p->sp_a, d->w, d->w+p->nx, 0);
    // Find the best lam[i] to make zero
    best_ind = -1;
    best_w = 0.;
    for (i=0; i<p->nz; ++i) {
      // Make sure variable influences du
      if (d->w[i]==0.) continue;
      // Make sure removing the constraint decreases dual infeasibility
      if (d->w[i]>0. ? d->lam[i]>=0. : d->lam[i]<=0.) continue;
      // Skip if maximum infeasibility increases
      if (casadi_qp_du_check(d, i)>d->du) continue;
      // Check if best so far
      if (fabs(d->w[i])>best_w) {
        best_w = fabs(d->w[i]);
        best_ind = i;
      }
    }
    // Accept, if any
    if (best_ind>=0) {
      *sign = 0;
      casadi_qp_log(d, "Removed %lld to reduce |du|", best_ind);
      return best_ind;
    } else {
      return -1;
    }
  }

  template<typename T1>
  void casadi_qp_kkt(casadi_qp_data<T1>* d) {
    // Local variables
    casadi_int i, k;
    const casadi_int *h_colind, *h_row, *a_colind, *a_row, *at_colind, *at_row,
                     *kkt_colind, *kkt_row;
    const casadi_qp_prob<T1>* p = d->prob;
    // Extract sparsities
    a_row = (a_colind = p->sp_a+2) + p->nx + 1;
    at_row = (at_colind = p->sp_at+2) + p->na + 1;
    h_row = (h_colind = p->sp_h+2) + p->nx + 1;
    kkt_row = (kkt_colind = p->sp_kkt+2) + p->nz + 1;
    // Reset w to zero
    casadi_fill(d->w, p->nz, 0.);
    // Loop over rows of the (transposed) KKT
    for (i=0; i<p->nz; ++i) {
      // Copy row of KKT to w
      if (i<p->nx) {
        if (d->lam[i]==0) {
          for (k=h_colind[i]; k<h_colind[i+1]; ++k) d->w[h_row[k]] = d->nz_h[k];
          for (k=a_colind[i]; k<a_colind[i+1]; ++k) d->w[p->nx+a_row[k]] = d->nz_a[k];
        } else {
          d->w[i] = 1.;
        }
      } else {
        if (d->lam[i]==0) {
          d->w[i] = -1.;
        } else {
          for (k=at_colind[i-p->nx]; k<at_colind[i-p->nx+1]; ++k) {
            d->w[at_row[k]] = d->nz_at[k];
          }
        }
      }
      // Copy row to KKT, zero out w
      for (k=kkt_colind[i]; k<kkt_colind[i+1]; ++k) {
        d->nz_kkt[k] = d->w[kkt_row[k]];
        d->w[kkt_row[k]] = 0;
      }
    }
  }

  template<typename T1>
  void casadi_qp_kkt_column(casadi_qp_data<T1>* d, T1* kkt_i, casadi_int i,
                            casadi_int sign) {
    // Local variables
    casadi_int k;
    const casadi_int *h_colind, *h_row, *a_colind, *a_row, *at_colind, *at_row;
    const casadi_qp_prob<T1>* p = d->prob;
    // Extract sparsities
    a_row = (a_colind = p->sp_a+2) + p->nx + 1;
    at_row = (at_colind = p->sp_at+2) + p->na + 1;
    h_row = (h_colind = p->sp_h+2) + p->nx + 1;
    // Reset kkt_i to zero
    casadi_fill(kkt_i, p->nz, 0.);
    // Copy column of KKT to kkt_i
    if (i<p->nx) {
      if (sign==0) {
        for (k=h_colind[i]; k<h_colind[i+1]; ++k) kkt_i[h_row[k]] = d->nz_h[k];
        for (k=a_colind[i]; k<a_colind[i+1]; ++k) kkt_i[p->nx+a_row[k]] = d->nz_a[k];
      } else {
        kkt_i[i] = 1.;
      }
    } else {
      if (sign==0) {
        kkt_i[i] = -1.;
      } else {
        for (k=at_colind[i-p->nx]; k<at_colind[i-p->nx+1]; ++k) {
          kkt_i[at_row[k]] = d->nz_at[k];
        }
      }
    }
  }

  template<typename T1>
  T1 casadi_qp_kkt_dot(casadi_qp_data<T1>* d, const T1* v, casadi_int i, casadi_int sign) {
    // Local variables
    casadi_int k;
    const casadi_int *h_colind, *h_row, *a_colind, *a_row, *at_colind, *at_row;
    T1 r;
    const casadi_qp_prob<T1>* p = d->prob;
    // Extract sparsities
    a_row = (a_colind = p->sp_a+2) + p->nx + 1;
    at_row = (at_colind = p->sp_at+2) + p->na + 1;
    h_row = (h_colind = p->sp_h+2) + p->nx + 1;
    // Scalar product with the desired column
    if (i<p->nx) {
      if (sign==0) {
        r = 0.;
        for (k=h_colind[i]; k<h_colind[i+1]; ++k) r += v[h_row[k]] * d->nz_h[k];
        for (k=a_colind[i]; k<a_colind[i+1]; ++k) r += v[p->nx+a_row[k]] * d->nz_a[k];
      } else {
        r = v[i];
      }
    } else {
      if (sign==0) {
        r = -v[i];
      } else {
        r = 0.;
        for (k=at_colind[i-p->nx]; k<at_colind[i-p->nx+1]; ++k) {
          r += v[at_row[k]] * d->nz_at[k];
        }
      }
    }
    return r;
  }

  template<typename T1>
  void casadi_qp_kkt_residual(casadi_qp_data<T1>* d, T1* r) {
    casadi_int i;
    const casadi_qp_prob<T1>* p = d->prob;
    for (i=0; i<p->nz; ++i) {
      if (d->lam[i]>0.) {
        r[i] = d->ubz[i]-d->z[i];
      } else if (d->lam[i]<0.) {
        r[i] = d->lbz[i]-d->z[i];
      } else if (i<p->nx) {
        r[i] = d->lam[i]-d->infeas[i];
      } else {
        r[i] = d->lam[i];
      }
    }
  }

  template<typename T1>
  int casadi_qp_zero_blocking(casadi_qp_data<T1>* d, T1 e,
                              casadi_int* index, casadi_int* sign) {
    // Local variables
    T1 dz_max;
    casadi_int i;
    int ret;
    const casadi_qp_prob<T1>* p = d->prob;
    ret = 0;
    dz_max = 0.;
    for (i=0; i<p->nz; ++i) {
      if (-d->dz[i]>dz_max && d->z[i]<=d->lbz[i]-e) {
        ret = 1;
        if (index) *index = i;
        if (sign) *sign = -1;
        casadi_qp_log(d, "lbz[%lld] violated at 0", i);
      } else if (d->dz[i]>dz_max && d->z[i]>=d->ubz[i]+e) {
        ret = 1;
        if (index) *index = i;
        if (sign) *sign = 1;
        casadi_qp_log(d, "ubz[%lld] violated at 0", i);
      }
    }
    return ret;
  }

  template<typename T1>
  void casadi_qp_primal_blocking(casadi_qp_data<T1>* d, T1 e,
                                 casadi_int* index, casadi_int* sign) {
    // Local variables
    casadi_int i;
    T1 trial_z;
    const casadi_qp_prob<T1>* p = d->prob;
    // Check if violation with tau=0 and not improving
    if (casadi_qp_zero_blocking(d, e, index, sign)) {
      d->tau = 0.;
      return;
    }
    // Loop over all primal variables
    for (i=0; i<p->nz; ++i) {
      if (d->dz[i]==0.) continue; // Skip zero steps
      // Trial primal step
      trial_z=d->z[i] + d->tau*d->dz[i];
      if (d->dz[i]<0 && trial_z<d->lbz[i]-e) {
        // Trial would increase maximum infeasibility
        d->tau = (d->lbz[i]-e-d->z[i])/d->dz[i];
        if (index) *index = d->lam[i]<0. ? -1 : i;
        if (sign) *sign = -1;
        casadi_qp_log(d, "Enforcing lbz[%lld]", i);
      } else if (d->dz[i]>0 && trial_z>d->ubz[i]+e) {
        // Trial would increase maximum infeasibility
        d->tau = (d->ubz[i]+e-d->z[i])/d->dz[i];
        if (index) *index = d->lam[i]>0. ? -1 : i;
        if (sign) *sign = 1;
        casadi_qp_log(d, "Enforcing ubz[%lld]", i);
      }
      if (d->tau<=0) return;
    }
  }

  template<typename T1>
  casadi_int casadi_qp_dual_breakpoints(casadi_qp_data<T1>* d, T1* tau_list,
                                        casadi_int* ind_list, T1 e, T1 tau) {
    // Local variables
    casadi_int i, n_tau, loc, next_ind, tmp_ind, j;
    T1 trial_lam, new_tau, next_tau, tmp_tau;
    const casadi_qp_prob<T1>* p = d->prob;
    // Dual feasibility is piecewise linear. Start with one interval [0,tau]:
    tau_list[0] = tau;
    ind_list[0] = -1; // no associated index
    n_tau = 1;
    // Find the taus corresponding to lam crossing zero and insert into list
    for (i=0; i<p->nz; ++i) {
      if (d->dlam[i]==0.) continue; // Skip zero steps
      if (d->lam[i]==0.) continue; // Skip inactive constraints
      // Trial dual step
      trial_lam = d->lam[i] + tau*d->dlam[i];
      // Skip if no sign change
      if (d->lam[i]>0 ? trial_lam>=0 : trial_lam<=0) continue;
      // Location of the sign change
      new_tau = -d->lam[i]/d->dlam[i];
      // Where to insert the w[i]
      for (loc=0; loc<n_tau-1; ++loc) {
        if (new_tau<tau_list[loc]) break;
      }
      // Insert element
      n_tau++;
      next_tau=new_tau;
      next_ind=i;
      for (j=loc; j<n_tau; ++j) {
        tmp_tau = tau_list[j];
        tau_list[j] = next_tau;
        next_tau = tmp_tau;
        tmp_ind = ind_list[j];
        ind_list[j] = next_ind;
        next_ind = tmp_ind;
      }
    }
    return n_tau;
  }

  template<typename T1>
  casadi_int casadi_qp_dual_blocking(casadi_qp_data<T1>* d, T1 e) {
    // Local variables
    casadi_int i, n_tau, j, k, du_index;
    T1 tau_k, dtau, new_infeas, tau1;
    const casadi_int *at_colind, *at_row;
    const casadi_qp_prob<T1>* p = d->prob;
    // Extract sparsities
    at_row = (at_colind = p->sp_at+2) + p->na + 1;
    // Dual feasibility is piecewise linear in tau. Get the intervals:
    n_tau = casadi_qp_dual_breakpoints(d, d->w, d->iw, e, d->tau);
    // No dual blocking yet
    du_index = -1;
    // How long step can we take without exceeding e?
    tau_k = 0.;
    for (j=0; j<n_tau; ++j) {
      // Distance to the next tau (may be zero)
      dtau = d->w[j] - tau_k;
      // Check if maximum dual infeasibilty gets exceeded
      for (k=0; k<p->nx; ++k) {
        new_infeas = d->infeas[k]+dtau*d->tinfeas[k];
        if (fabs(new_infeas)>e) {
          tau1 = fmax(0., tau_k + ((new_infeas>0 ? e : -e)-d->infeas[k])/d->tinfeas[k]);
          if (tau1 < d->tau) {
            // Smallest tau found so far
            d->tau = tau1;
            du_index = k;
          }
        }
      }
      // Update infeasibility
      casadi_axpy(p->nx, fmin(d->tau - tau_k, dtau), d->tinfeas, d->infeas);
      // Stop here if dual blocking constraint
      if (du_index>=0) return du_index;
      // Continue to the next tau
      tau_k = d->w[j];
      // Get component, break if last
      i = d->iw[j];
      if (i<0) break;
      // Update sign or tinfeas
      if (!d->neverzero[i]) {
        // lam becomes zero, update the infeasibility tangent
        if (i<p->nx) {
          // Set a lam_x to zero
          d->tinfeas[i] -= d->dlam[i];
        } else {
          // Set a lam_a to zero
          for (k=at_colind[i-p->nx]; k<at_colind[i-p->nx+1]; ++k) {
            d->tinfeas[at_row[k]] -= d->nz_at[k]*d->dlam[i];
          }
        }
      }
    }
    return du_index;
  }

  template<typename T1>
  void casadi_qp_take_step(casadi_qp_data<T1>* d) {
    // Local variables
    casadi_int i;
    const casadi_qp_prob<T1>* p = d->prob;
    // Get current sign
    for (i=0; i<p->nz; ++i) d->iw[i] = d->lam[i]>0. ? 1 : d->lam[i]<0 ? -1 : 0;
    // Take primal-dual step
    casadi_axpy(p->nz, d->tau, d->dz, d->z);
    casadi_axpy(p->nz, d->tau, d->dlam, d->lam);
    // Update sign
    for (i=0; i<p->nz; ++i) {
      // Allow sign changes for certain components
      if (d->neverzero[i] && (d->iw[i]<0 ? d->lam[i]>0 : d->lam[i]<0)) {
        d->iw[i]=-d->iw[i];
      }
      // Ensure correct sign
      switch (d->iw[i]) {
        case -1: d->lam[i] = fmin(d->lam[i], -p->dmin); break;
        case  1: d->lam[i] = fmax(d->lam[i],  p->dmin); break;
        case  0: d->lam[i] = 0.; break;
      }
    }
  }

  template<typename T1>
  int casadi_qp_flip_check(casadi_qp_data<T1>* d, casadi_int index, casadi_int sign,
                           casadi_int* r_index, casadi_int* r_sign, T1 e) {
    // Local variables
    casadi_int new_sign, i;
    T1 best_slack, new_slack;
    const casadi_qp_prob<T1>* p = d->prob;
    // New column that we're trying to add
    casadi_qp_kkt_column(d, d->dz, index, sign);
    // Express it using the other columns
    casadi_qr_solve(d->dz, 1, 0, p->sp_v, d->nz_v, p->sp_r, d->nz_r, d->beta,
                    p->prinv, p->pc, d->w);
    // Quick return if columns are linearly independent
    if (fabs(d->dz[index])>=1e-12) return 0;
    // Column that we're removing
    casadi_qp_kkt_column(d, d->w, index, !sign);
    // Find best constraint we can flip, if any
    *r_index=-1;
    *r_sign=0;
    best_slack = -inf;
    for (i=0; i<p->nz; ++i) {
      // Can't be the same
      if (i==index) continue;
      // Make sure constraint is flippable
      if (d->lam[i]==0 ? d->neverlower[i] && d->neverupper[i] : d->neverzero[i]) continue;
      // If dz[i]!=0, column i is redundant
      if (fabs(d->dz[i])<1e-12) continue;
      // We want to make sure that the new, flipped column i is linearly
      // independent with other columns. We have:
      // flipped_column[index] = dz[0]*column[0] + ... + dz[N-1]*column[N-1]
      // We also require that flipped_column[i] isn't othogonal to
      // (old) column[index], as this will surely lead to singularity
      // This will not cover all cases of singularity, but many important
      // ones. General singularity handling is done below.
      if (fabs(casadi_qp_kkt_dot(d, d->w, i, d->lam[i]==0.)) < 1e-12) continue;
      // Dual infeasibility
      // Check if the best so far
      if (d->lam[i]==0.) {
        // Which bound is closer?
        new_sign = d->lbz[i]-d->z[i] >= d->z[i]-d->ubz[i] ? -1 : 1;
        // Better than negative slack, worse than positive slack
        new_slack = 0;
      } else {
        // Skip if flipping would result in too large |du|
        if (casadi_qp_du_check(d, i)>e) continue;
        // Slack to the bound
        new_slack = d->lam[i]>0 ? d->ubz[i]-d->z[i] : d->z[i]-d->lbz[i];
        new_sign = 0;
      }
      // Best so far?
      if (new_slack > best_slack) {
        best_slack = new_slack;
        *r_index = i;
        *r_sign = new_sign;
      }
    }
    // Accept, if any
    return *r_index>=0 ? 0 : 1;
  }

  template<typename T1>
  void casadi_qp_factorize(casadi_qp_data<T1>* d) {
    const casadi_qp_prob<T1>* p = d->prob;
    // Construct the KKT matrix
    casadi_qp_kkt(d);
    // QR factorization
    casadi_qr(p->sp_kkt, d->nz_kkt, d->w, p->sp_v, d->nz_v, p->sp_r,
              d->nz_r, d->beta, p->prinv, p->pc);
    // Check singularity
    d->sing = casadi_qr_singular(&d->mina, &d->imina, d->nz_r, p->sp_r, p->pc, 1e-12);
  }

  template<typename T1>
  int casadi_qp_scale_step(casadi_qp_data<T1>* d, casadi_int* r_index, casadi_int* r_sign) {
    // Local variables
    T1 tpr, tdu, terr, tau_test, minat_tr, tau;
    int pos_ok, neg_ok;
    casadi_int nnz_kkt, nullity_tr, nulli, imina_tr, i;
    const casadi_qp_prob<T1>* p = d->prob;
    // Reset r_index, r_sign, quick return if non-singular
    *r_index = -1;
    *r_sign = 0;
    if (!d->sing) {
      return 0;
    }
    // Change in pr, du in the search direction
    tpr = d->ipr<0 ? 0.
                   : d->z[d->ipr]>d->ubz[d->ipr] ? d->dz[d->ipr]/d->pr
                                                 : -d->dz[d->ipr]/d->pr;
    tdu = d->idu<0 ? 0. : d->tinfeas[d->idu]/d->infeas[d->idu];
    // Change in max(pr, du) in the search direction
    pos_ok=1, neg_ok=1;
    if (d->pr>d->du) {
      // |pr|>|du|
      if (tpr<0) {
        neg_ok = 0;
      } else if (tpr>0) {
        pos_ok = 0;
      }
      terr = tpr;
    } else if (d->pr<d->du) {
      // |pr|<|du|
      if (tdu<0) {
        neg_ok = 0;
      } else if (tdu>0) {
        pos_ok = 0;
      }
      terr = tdu;
    } else {
      // |pr|==|du|
      if ((tpr>0 && tdu<0) || (tpr<0 && tdu>0)) {
        // |pr|==|du| cannot be decreased along the search direction
        pos_ok = neg_ok = 0;
        terr = 0;
      } else if (fmin(tpr, tdu)<0) {
        // |pr|==|du| decreases for positive tau
        neg_ok = 0;
        terr = fmax(tpr, tdu);
      } else if (fmax(tpr, tdu)>0) {
        // |pr|==|du| decreases for negative tau
        pos_ok = 0;
        terr = fmin(tpr, tdu);
      } else {
        terr = 0;
      }
    }
    // If primal error is dominating and constraint is active,
    // then only allow the multiplier to become larger
    if (p->du_to_pr*d->pr>=d->du && d->lam[d->ipr]!=0 && fabs(d->dlam[d->ipr])>1e-12) {
      if ((d->lam[d->ipr]>0)==(d->dlam[d->ipr]>0)) {
        neg_ok = 0;
      } else {
        pos_ok = 0;
      }
    }
    // QR factorization of the transpose
    casadi_trans(d->nz_kkt, p->sp_kkt, d->nz_v, p->sp_kkt, d->iw);
    nnz_kkt = p->sp_kkt[2+p->nz]; // kkt_colind[nz]
    casadi_copy(d->nz_v, nnz_kkt, d->nz_kkt);
    casadi_qr(p->sp_kkt, d->nz_kkt, d->w, p->sp_v, d->nz_v, p->sp_r, d->nz_r,
              d->beta, p->prinv, p->pc);
    // Best flip
    tau = p->inf;
    // For all nullspace vectors
    nullity_tr = casadi_qr_singular(&minat_tr, &imina_tr, d->nz_r, p->sp_r,
                                    p->pc, 1e-12);
    for (nulli=0; nulli<nullity_tr; ++nulli) {
      // Get a linear combination of the rows in kkt
      casadi_qr_colcomb(d->w, d->nz_r, p->sp_r, p->pc, imina_tr, nulli);
      // Look for the best constraint for increasing rank
      for (i=0; i<p->nz; ++i) {
        // Check if old column can be removed without decreasing rank
        if (fabs(i<p->nx ? d->dz[i] : d->dlam[i])<1e-12) continue;
        // If dot(w, kkt(i)-kkt_flipped(i))==0, rank won't increase
        if (fabs(casadi_qp_kkt_dot(d, d->w, i, 0)
                 - casadi_qp_kkt_dot(d, d->w, i, 1))<1e-12) continue;
        // Is constraint active?
        if (d->lam[i]==0.) {
          // Make sure that step is nonzero
          if (fabs(d->dz[i])<1e-12) continue;
          // Step needed to bring z to lower bound
          if (!d->neverlower[i]) {
            tau_test = (d->lbz[i]-d->z[i])/d->dz[i];
            // Ensure nonincrease in max(pr, du)
            if (!((terr>0. && tau_test>0.) || (terr<0. && tau_test<0.))) {
              // Only allow removing constraints if tau_test==0
              if (fabs(tau_test)>=1e-16) {
                // Check if best so far
                if (fabs(tau_test)<fabs(tau)) {
                  tau = tau_test;
                  *r_index = i;
                  *r_sign = -1;
                  casadi_qp_log(d, "Enforced lbz[%lld] for regularity", i);
                }
              }
            }
          }
          // Step needed to bring z to upper bound
          if (!d->neverupper[i]) {
            tau_test = (d->ubz[i]-d->z[i])/d->dz[i];
            // Ensure nonincrease in max(pr, du)
            if (!((terr>0. && tau_test>0.) || (terr<0. && tau_test<0.))) {
              // Only allow removing constraints if tau_test==0
              if (fabs(tau_test)>=1e-16) {
                // Check if best so far
                if (fabs(tau_test)<fabs(tau)) {
                  tau = tau_test;
                  *r_index = i;
                  *r_sign = 1;
                  casadi_qp_log(d, "Enforced ubz[%lld] for regularity", i);
                }
              }
            }
          }
        } else {
          // Make sure that step is nonzero
          if (fabs(d->dlam[i])<1e-12) continue;
          // Step needed to bring lam to zero
          if (!d->neverzero[i]) {
            tau_test = -d->lam[i]/d->dlam[i];
            // Ensure nonincrease in max(pr, du)
            if ((terr>0. && tau_test>0.) || (terr<0. && tau_test<0.)) continue;
            // Make sure direction is permitted
            if ((tau_test>0 && !pos_ok) || (tau_test<0 && !neg_ok)) continue;
            // Check if best so far
            if (fabs(tau_test)<fabs(tau)) {
              tau = tau_test;
              *r_index = i;
              *r_sign = 0;
              casadi_qp_log(d, "Dropped %s[%lld] for regularity",
                     d->lam[i]>0 ? "lbz" : "ubz", i);
            }
          }
        }
      }
    }
    // Can we restore feasibility?
    if (*r_index<0) return 1;
    // Scale step so that that tau=1 corresponds to a full step
    casadi_scal(p->nz, tau, d->dz);
    casadi_scal(p->nz, tau, d->dlam);
    casadi_scal(p->nx, tau, d->tinfeas);
    return 0;
  }

  template<typename T1>
  int casadi_qp_calc_step(casadi_qp_data<T1>* d, casadi_int* r_index, casadi_int* r_sign) {
    // Local variables
    casadi_int i;
    const casadi_qp_prob<T1>* p = d->prob;
    // Calculate step in z[:nx] and lam[nx:]
    if (!d->sing) {
      // Negative KKT residual
      casadi_qp_kkt_residual(d, d->dz);
      // Solve to get primal-dual step
      casadi_qr_solve(d->dz, 1, 1, p->sp_v, d->nz_v, p->sp_r, d->nz_r, d->beta,
                      p->prinv, p->pc, d->w);
    } else {
      // Get a linear combination of the columns in KKT
      casadi_qr_colcomb(d->dz, d->nz_r, p->sp_r, p->pc, d->imina, 0);
    }
    // Calculate change in Lagrangian gradient
    casadi_fill(d->dlam, p->nx, 0.);
    casadi_mv(d->nz_h, p->sp_h, d->dz, d->dlam, 0); // gradient of the objective
    casadi_mv(d->nz_a, p->sp_a, d->dz+p->nx, d->dlam, 1); // gradient of the Lagrangian
    // Step in lam[:nx]
    casadi_scal(p->nx, -1., d->dlam);
    // For inactive constraints, lam(x) step is zero
    for (i=0; i<p->nx; ++i) if (d->lam[i]==0.) d->dlam[i] = 0.;
    // Step in lam[nx:]
    casadi_copy(d->dz+p->nx, p->na, d->dlam+p->nx);
    // Step in z[nx:]
    casadi_fill(d->dz+p->nx, p->na, 0.);
    casadi_mv(d->nz_a, p->sp_a, d->dz, d->dz+p->nx, 0);
    // Avoid steps that are nonzero due to numerics
    for (i=0; i<p->nz; ++i) if (fabs(d->dz[i])<1e-14) d->dz[i] = 0.;
    // Tangent of the dual infeasibility at tau=0
    casadi_fill(d->tinfeas, p->nx, 0.);
    casadi_mv(d->nz_h, p->sp_h, d->dz, d->tinfeas, 0);
    casadi_mv(d->nz_a, p->sp_a, d->dlam+p->nx, d->tinfeas, 1);
    casadi_axpy(p->nx, 1., d->dlam, d->tinfeas);
    // Calculate step length
    return casadi_qp_scale_step(d, r_index, r_sign);
  }

  template<typename T1>
  void casadi_qp_calc_dependent(casadi_qp_data<T1>* d) {
    // Local variables
    casadi_int i;
    const casadi_qp_prob<T1>* p = d->prob;
    // Calculate f
    d->f = casadi_bilin(d->nz_h, p->sp_h, d->z, d->z)/2.
         + casadi_dot(p->nx, d->z, d->g);
    // Calculate z[nx:]
    casadi_fill(d->z+p->nx, p->na, 0.);
    casadi_mv(d->nz_a, p->sp_a, d->z, d->z+p->nx, 0);
    // Calculate gradient of the Lagrangian
    casadi_copy(d->g, p->nx, d->infeas);
    casadi_mv(d->nz_h, p->sp_h, d->z, d->infeas, 0);
    casadi_mv(d->nz_a, p->sp_a, d->lam+p->nx, d->infeas, 1);
    // Calculate lam[:nx] without changing the sign, dual infeasibility
    for (i=0; i<p->nx; ++i) {
      if (d->lam[i]>0) {
        d->lam[i] = fmax(-d->infeas[i], p->dmin);
      } else if (d->lam[i]<0) {
        d->lam[i] = fmin(-d->infeas[i], -p->dmin);
      }
      d->infeas[i] += d->lam[i];
    }
    // Calculate primal and dual error
    casadi_qp_pr(d);
    casadi_qp_du(d);
  }

  template<typename T1>
  void casadi_qp_linesearch(casadi_qp_data<T1>* d, casadi_int* index, casadi_int* sign) {
    const casadi_qp_prob<T1>* p = d->prob;
    // Start with a full step and no active set change
    *sign=0;
    *index=-1;
    d->tau = 1.;
    // Find largest possible step without exceeding acceptable |pr|
    casadi_qp_primal_blocking(d, fmax(d->pr, d->du/p->du_to_pr), index, sign);
    // Find largest possible step without exceeding acceptable |du|
    if (casadi_qp_dual_blocking(d, fmax(d->pr*p->du_to_pr, d->du))>=0) {
      *index = -1;
      *sign=0;
    }
    // Take primal-dual step, avoiding accidental sign changes for lam
    casadi_qp_take_step(d);
  }

  template<typename T1>
  void casadi_qp_flip(casadi_qp_data<T1>* d, casadi_int *index, casadi_int *sign,
                                            casadi_int r_index, casadi_int r_sign) {
    // Local variables
    T1 e;
    const casadi_qp_prob<T1>* p = d->prob;
    // acceptable dual error
    e = fmax(p->du_to_pr*d->pr, d->du);
    // Try to restore regularity if possible
    if (r_index>=0 && (r_sign!=0 || casadi_qp_du_check(d, r_index)<=e)) {
      *index = r_index;
      *sign = r_sign;
      casadi_qp_log(d, "%lld->%lld for regularity", *index, *sign);
    }
    // Improve primal or dual feasibility
    if (*index==-1 && d->tau>1e-16 && (d->ipr>=0 || d->idu>=0)) {
      if (p->du_to_pr*d->pr >= d->du) {
        // Try to improve primal feasibility
        *index = casadi_qp_pr_index(d, sign);
      } else {
        // Try to improve dual feasibility
        *index = casadi_qp_du_index(d, sign);
      }
    }
    // If a constraint was added
    if (*index>=0) {
      // Try to maintain non-singularity if possible
      if (!d->sing && casadi_qp_flip_check(d, *index, *sign, &r_index, &r_sign, e)) {
        if (r_index>=0) {
          // Also flip r_index to avoid singularity
          d->lam[r_index] = r_sign==0 ? 0 : r_sign>0 ? p->dmin : -p->dmin;
          casadi_qp_log(d, "%lld->%lld, %lld->%lld", *index, *sign, r_index, r_sign);
        }
      }
      d->lam[*index] = *sign==0 ? 0 : *sign>0 ? p->dmin : -p->dmin;
      // Recalculate primal and dual infeasibility
      casadi_qp_calc_dependent(d);
      // Reset index
      *index=-2;
    }
  }


  int ConicActiveSet::init_mem(void* mem) const {
    //auto m = static_cast<ConicActiveSetMemory*>(mem);
    return 0;
  }

  int ConicActiveSet::
  eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<ConicActiveSetMemory*>(mem);

    // Statistics
    for (auto&& s : m->fstats) s.second.reset();

    if (inputs_check_) {
      check_inputs(arg[CONIC_LBX], arg[CONIC_UBX], arg[CONIC_LBA], arg[CONIC_UBA]);
    }

    // Setup memory structure
    casadi_qp_data<double> qp_m;
    qp_m.prob = &p_;

    // Local variables
    int flag;
    casadi_int i;
    // Get input pointers
    const double *h, *g, *a, *lba, *uba, *lbx, *ubx, *x0, *lam_x0, *lam_a0;
    h = arg[CONIC_H];
    g = arg[CONIC_G];
    a = arg[CONIC_A];
    lba = arg[CONIC_LBA];
    uba = arg[CONIC_UBA];
    lbx = arg[CONIC_LBX];
    ubx = arg[CONIC_UBX];
    x0 = arg[CONIC_X0];
    lam_x0 = arg[CONIC_LAM_X0];
    lam_a0 = arg[CONIC_LAM_A0];

    // Get output pointers
    double *x, *f, *lam_a, *lam_x;
    x = res[CONIC_X];
    f = res[CONIC_COST];
    lam_a = res[CONIC_LAM_A];
    lam_x = res[CONIC_LAM_X];

    // Work vectors
    double *kkt, *z, *lam, *v, *r, *beta, *dz, *dlam, *lbz, *ubz,
           *trans_a, *infeas, *tinfeas, *vr;
    kkt = w; w += kkt_.nnz();
    z = w; w += nx_+na_;
    lbz = w; w += nx_+na_;
    ubz = w; w += nx_+na_;
    lam = w; w += nx_+na_;
    dz = w; w += nx_+na_;
    dlam = w; w += nx_+na_;
    vr = w; w += max(sp_v_.nnz()+sp_r_.nnz(), kkt_.nnz());
    v = vr;
    r = vr + sp_v_.nnz();
    beta = w; w += nx_+na_;
    trans_a = w; w += AT_.nnz();
    infeas = w; w += nx_;
    tinfeas = w; w += nx_;
    casadi_int *neverzero, *neverupper, *neverlower;
    neverzero = iw; iw += nx_+na_;
    neverupper = iw; iw += nx_+na_;
    neverlower = iw; iw += nx_+na_;

    // Bounds on z
    casadi_copy(lbx, nx_, lbz);
    casadi_copy(lba, na_, lbz+nx_);
    casadi_copy(ubx, nx_, ubz);
    casadi_copy(uba, na_, ubz+nx_);

    // Pass initial guess
    casadi_copy(x0, nx_, z);
    casadi_copy(lam_x0, nx_, lam);
    casadi_copy(lam_a0, na_, lam+nx_);

    // Correct lam if needed, determine permitted signs
    for (i=0; i<nx_+na_; ++i) {
      // Permitted signs for lam
      neverzero[i] = lbz[i]==ubz[i];
      neverupper[i] = isinf(ubz[i]);
      neverlower[i] = isinf(lbz[i]);
      casadi_assert(!neverzero[i] || !neverupper[i] || !neverlower[i],
                    "No sign possible for " + str(i));
      // Correct initial active set if required
      if (neverzero[i] && lam[i]==0.) {
        lam[i] = neverupper[i] || z[i]-lbz[i] <= ubz[i]-z[i] ? -p_.dmin : p_.dmin;
      } else if (neverupper[i] && lam[i]>0.) {
        lam[i] = neverzero[i] ? -p_.dmin : 0.;
      } else if (neverlower[i] && lam[i]<0.) {
        lam[i] = neverzero[i] ? p_.dmin : 0.;
      }
    }
    // Transpose A
    casadi_trans(a, A_, trans_a, AT_, iw);
    qp_m.msg[0] = '\0';
    qp_m.z = z;
    qp_m.lam = lam;
    qp_m.lbz = lbz;
    qp_m.ubz = ubz;
    qp_m.infeas = infeas;
    qp_m.tinfeas = tinfeas;
    qp_m.dz = dz;
    qp_m.dlam = dlam;
    qp_m.w = w;
    qp_m.iw = iw;
    qp_m.neverzero = neverzero;
    qp_m.neverlower = neverlower;
    qp_m.neverupper = neverupper;
    // Matrix nonzeros
    qp_m.g = g;
    qp_m.nz_a = a;
    qp_m.nz_at = trans_a;
    qp_m.nz_h = h;
    qp_m.nz_kkt = kkt;
    // QR factorization
    qp_m.nz_v = v;
    qp_m.nz_r = r;
    qp_m.beta = beta;
    // Misc
    qp_m.tau = 0.;
    qp_m.sing = 0; // set to avoid false positive warning
    // Constraint to be flipped, if any
    casadi_int index=-2, sign=0, r_index=-2, r_sign=0;
    // QP iterations
    casadi_int iter = 0;
    while (true) {
      // Calculate dependent quantities
      casadi_qp_calc_dependent(&qp_m);
      // Make an active set change
      casadi_qp_flip(&qp_m, &index, &sign, r_index, r_sign);
      // Form and factorize the KKT system
      casadi_qp_factorize(&qp_m);
      // Print iteration progress:
      if (print_iter_) {
        if (iter % 10 == 0) {
          print("%5s %5s %9s %9s %5s %9s %5s %9s %5s %9s %40s\n",
                "Iter", "Sing", "fk", "|pr|", "con", "|du|", "var",
                "min_R", "con", "last_tau", "Note");
        }
        print("%5d %5d %9.2g %9.2g %5d %9.2g %5d %9.2g %5d %9.2g %40s\n",
              iter, qp_m.sing, qp_m.f, qp_m.pr, qp_m.ipr, qp_m.du, qp_m.idu,
              qp_m.mina, qp_m.imina, qp_m.tau, qp_m.msg);
        qp_m.msg[0] = '\0';
      }
      // Successful return if still no change
      if (index==-1) {
        flag = 0;
        break;
      }
      // Too many iterations?
      if (++iter>max_iter_) {
        casadi_warning("Maximum number of iterations reached");
        flag = 1;
        break;
      }
      // Start new iteration
      // Calculate search direction
      if (casadi_qp_calc_step(&qp_m, &r_index, &r_sign)) {
        casadi_warning("Failed to calculate search direction");
        flag = 1;
        break;
      }
      // Line search in the calculated direction
      casadi_qp_linesearch(&qp_m, &index, &sign);
    }

    // Calculate optimal cost
    if (f) *f = qp_m.f;

    // Get solution
    casadi_copy(z, nx_, x);
    casadi_copy(lam, nx_, lam_x);
    casadi_copy(lam+nx_, na_, lam_a);

    return flag;
  }

} // namespace casadi
