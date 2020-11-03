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


#include "qpchasm.hpp"
#include "casadi/core/nlpsol.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_CONIC_QPCHASM_EXPORT
  casadi_register_conic_qpchasm(Conic::Plugin* plugin) {
    plugin->creator = Qpchasm::creator;
    plugin->name = "qpchasm";
    plugin->doc = Qpchasm::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &Qpchasm::options_;
    plugin->deserialize = &Qpchasm::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_QPCHASM_EXPORT casadi_load_conic_qpchasm() {
    Conic::registerPlugin(casadi_register_conic_qpchasm);
  }

  void print_vec(const std::string& str, const double* v, casadi_int n) {
    uout() << str << " =";
    for (casadi_int k = 0; k < n; ++k) uout() << " " << v[k];
    uout() << "\n";
  }

  template<typename T1>
  void init_ip(casadi_qpip_data<T1>* d) {
    // Local variables
    casadi_int k;
    T1 margin, mid;
    const casadi_qp_prob<T1>* p = d->prob;
    // Required margin to constraints
    margin = .1;
    // Reset constraint count
    d->n_con = 0;
    // Initialize constraints to zero
    for (k = p->nx; k < p->nz; ++k) d->z[k] = 0;
    // Find interior point
    for (k = 0; k < p->nz; ++k) {
      if (d->lbz[k] > -p->inf) {
        if (d->ubz[k] < p->inf) {
          // Both upper and lower bounds
          mid = .5 * (d->lbz[k] + d->ubz[k]);
          // Ensure margin to boundary, without crossing midpoint
          if (d->z[k] < mid) {
            d->z[k] = fmin(fmax(d->z[k], d->lbz[k] + margin), mid);
          } else if (d->z[k] > mid) {
            d->z[k] = fmax(fmin(d->z[k], d->ubz[k] - margin), mid);
          }
          if (d->ubz[k] > d->lbz[k] + p->dmin) {
            d->lam_lbz[k] = 1;
            d->lam_ubz[k] = 1;
            d->n_con += 2;
          }
        } else {
          // Only lower bound
          d->z[k] = fmax(d->z[k], d->lbz[k] + margin);
          d->lam_lbz[k] = 1;
          d->n_con++;
        }
      } else {
        if (d->ubz[k] < p->inf) {
          // Only upper bound
          d->z[k] = fmin(d->z[k], d->ubz[k] - margin);
          d->lam_ubz[k] = 1;
          d->n_con++;
        }
      }
    }
    // Reset iteration counter
    d->iter = 0;
    // Next task
    d->task = QP_MV;
    d->next = QP_RESIDUAL;
  }

  template<typename T1>
  void calc_diag(casadi_qpip_data<T1>* d) {
    // Local variables
    casadi_int k;
    const casadi_qp_prob<T1>* p = d->prob;
    // Diagonal entries corresponding to variables
    for (k = 0; k < p->nx; ++k) {
      if (d->ubz[k] <= d->lbz[k] + p->dmin) {
        // Fixed variable (eliminate)
        d->D[k] = -1;
      } else {
        d->D[k] = d->lam_lbz[k] * d->dinv_lbz[k]
          + d->lam_ubz[k] * d->dinv_ubz[k];
      }
    }
    // Diagonal entries corresponding to constraints
    for (; k < p->nz; ++k) {
      if (d->lbz[k] <= -p->inf && d->ubz[k] >= p->inf) {
        // Unconstrained (eliminate)
        d->D[k] = -1;
      } else if (d->ubz[k] <= d->lbz[k] + p->dmin) {
        // Equality constrained
        d->D[k] = 0;
      } else {
        d->D[k] = 1. / (d->lam_lbz[k] * d->dinv_lbz[k]
          + d->lam_ubz[k] * d->dinv_ubz[k]);
      }
    }
    // Scale diagonal entries
    for (k = 0; k < p->nz; ++k) {
      if (d->D[k] < 0) {
        // Eliminate
        d->S[k] = 0;
        d->D[k] = 1;
      } else {
        // Scale
        d->S[k] = fmin(1., std::sqrt(1. / d->D[k]));
        d->D[k] = fmin(1., d->D[k]);
      }
    }
  }

  template<typename T1>
  int qp_ip_iter(casadi_qpip_data<T1>* d) {
    // Local variables
    const casadi_qp_prob<T1>* p = d->prob;
    // Stop, if converged or max iter
    if (d->iter >= p->max_iter) return 1;
    // Start new iteration
    d->iter++;
    // Calculate diagonal entries and scaling factors
    calc_diag(d);
    // Success
    return 0;
  }

  template<typename T1>
  void calc_res(casadi_qpip_data<T1>* d) {
    // Local variables
    casadi_int k;
    T1 bdiff;
    const casadi_qp_prob<T1>* p = d->prob;
    // Gradient of the Lagrangian
    casadi_axpy(p->nx, 1., d->g, d->rz);
    for (k = 0; k < p->nx; ++k) {
      if (d->ubz[k] <= d->lbz[k] + p->dmin) {
        // Fixed variable: Solve to get multiplier explicitly
        d->lam[k] = -d->rz[k];
        d->rz[k] = 0;
      } else {
        // Residual
        d->rz[k] += d->lam[k];
      }
    }
    // Constraint violation (only possible for linear constraints)
    d->ipr = -1;
    d->pr = 0;
    for (k = p->na; k < p->nz; ++k) {
      if (d->rz[k] + d->pr < d->lbz[k]) {
        d->pr = d->lbz[k] - d->rz[k];
        d->ipr = k;
      } else if (d->rz[k] - d->pr > d->ubz[k]) {
        d->pr = d->rz[k] - d->ubz[k];
        d->ipr = k;
      }
    }
    // Dual infeasibility
    d->idu = -1;
    d->du = 0;
    // Linear constraint
    casadi_axpy(p->na, -1., d->z + p->nx, d->rz + p->nx);
    // Multiplier consistency
    for (k = 0; k < p->nz; ++k) {
      if (d->ubz[k] <= d->lbz[k] + p->dmin) {
        // Fixed variable: Solve to get lam_lbz, lam_ubz
        d->lam_ubz[k] = fmax(d->lam[k], 0.);
        d->lam_lbz[k] = fmax(-d->lam[k], 0.);
        d->rlam[k] = 0;
      } else {
        // Residual
        d->rlam[k] = d->lam_ubz[k] - d->lam_lbz[k] - d->lam[k];
        // Largest dual infeasibility
        if (fabs(d->rlam[k]) > d->du) {
          d->du = fabs(d->rlam[k]);
          d->idu = k;
        }
      }
    }
    // Complementarity conditions, mu
    d->mu = 0;
    for (k = 0; k < p->nz; ++k) {
      // Lower bound
      if (d->lbz[k] > -p->inf && d->ubz[k] > d->lbz[k] + p->dmin) {
        bdiff = d->z[k] - d->lbz[k];
        d->mu += d->rlam_lbz[k] = d->lam_lbz[k] * bdiff;
        d->dinv_lbz[k] = 1. / bdiff;
      } else {
        d->rlam_lbz[k] = 0;
        d->dinv_lbz[k] = 0;
      }
      // Upper bound
      if (d->ubz[k] < p->inf && d->ubz[k] > d->lbz[k] + p->dmin) {
        bdiff = d->ubz[k] - d->z[k];
        d->mu += d->rlam_ubz[k] = d->lam_ubz[k] * bdiff;
        d->dinv_ubz[k] = 1. / bdiff;
      } else {
        d->rlam_ubz[k] = 0;
        d->dinv_ubz[k] = 0;
      }
    }
    // Divide mu by total number of finite constraints
    if (d->n_con > 0) d->mu /= d->n_con;
  }

  template<typename T1>
  void qp_factorize(casadi_qpip_data<T1>* d) {
    // Local variables
    casadi_int i, k, j;
    const casadi_int *h_colind, *h_row, *a_colind, *a_row, *at_colind, *at_row,
                     *kkt_colind, *kkt_row;
    const casadi_qp_prob<T1>* p = d->prob;
    // Extract sparsities
    a_row = (a_colind = p->sp_a+2) + p->nx + 1;
    at_row = (at_colind = p->sp_at+2) + p->na + 1;
    h_row = (h_colind = p->sp_h+2) + p->nx + 1;
    kkt_row = (kkt_colind = p->sp_kkt+2) + p->nz + 1;
    // Reset w to zero
    casadi_clear(d->w, p->nz);
    // Loop over rows of the (transposed) KKT
    for (i=0; i<p->nz; ++i) {
      // Copy row of KKT to w
      if (i<p->nx) {
        for (k=h_colind[i]; k<h_colind[i+1]; ++k) d->w[h_row[k]] = d->nz_h[k];
        for (k=a_colind[i]; k<a_colind[i+1]; ++k) d->w[p->nx+a_row[k]] = d->nz_a[k];
      } else {
        for (k=at_colind[i-p->nx]; k<at_colind[i-p->nx+1]; ++k) {
          d->w[at_row[k]] = d->nz_at[k];
        }
      }
      // Copy row to KKT, scale, zero out w
      for (k=kkt_colind[i]; k<kkt_colind[i+1]; ++k) {
        j = kkt_row[k];
        d->nz_kkt[k] = d->S[j] * d->w[j] * d->S[i];
        d->w[j] = 0;
        if (i == j) {
          d->nz_kkt[k] += i<p->nx ? d->D[i] : -d->D[i];
        }
      }
    }
    // QR factorization
    casadi_qr(p->sp_kkt, d->nz_kkt, d->w, p->sp_v, d->nz_v, p->sp_r,
              d->nz_r, d->beta, p->prinv, p->pc);
    // Check singularity
    d->sing = casadi_qr_singular(&d->mina, &d->imina, d->nz_r, p->sp_r, p->pc, 1e-12);
  }

  template<typename T1>
  void qp_predictor_prepare(casadi_qpip_data<T1>* d) {
    // Local variables
    casadi_int k;
    const casadi_qp_prob<T1>* p = d->prob;
    // Store r_lam - dinv_lbz * rlam_lbz + dinv_ubz * rlam_ubz in dz
    casadi_copy(d->rlam, p->nz, d->dz);
    for (k=0; k<p->nz; ++k) d->dz[k] += d->dinv_lbz[k] * d->rlam_lbz[k];
    for (k=0; k<p->nz; ++k) d->dz[k] -= d->dinv_ubz[k] * d->rlam_ubz[k];
    // Finish calculating x-component of right-hand-side and store in dz[:nx]
    for (k=0; k<p->nx; ++k) d->dz[k] += d->rz[k];
    // Copy tilde{r}_lam to dlam[nx:] (needed to calculate step in g later)
    for (k=p->nx; k<p->nz; ++k) d->dlam[k] = d->dz[k];
    // Finish calculating g-component of right-hand-side and store in dz[nx:]
    for (k=p->nx; k<p->nz; ++k) {
      d->dz[k] *= d->D[k] / (d->S[k] * d->S[k]);
      d->dz[k] += d->rz[k];
    }
    // Scale and negate right-hand-side
    for (k=0; k<p->nz; ++k) d->dz[k] *= -d->S[k];
    // dlam_lbz := -rlam_lbz, dlam_ubz := -rlam_ubz
    for (k=0; k<p->nz; ++k) d->dlam_lbz[k] = -d->rlam_lbz[k];
    for (k=0; k<p->nz; ++k) d->dlam_ubz[k] = -d->rlam_ubz[k];
    // dlam_x := rlam_x
    for (k=0; k<p->nx; ++k) d->dlam[k] = d->rlam[k];
    // Solve to get step
    d->linsys = d->dz;
  }

  template<typename T1>
  void qp_predictor(casadi_qpip_data<T1>* d) {
    // Local variables
    casadi_int k;
    T1 t, alpha, sigma;
    const casadi_qp_prob<T1>* p = d->prob;
    // Scale results
    for (k=0; k<p->nz; ++k) d->dz[k] *= d->S[k];
    // Calculate step in z(g), lam(g)
    for (k=p->nx; k<p->nz; ++k) {
      t = d->D[k] / (d->S[k] * d->S[k]) * (d->dz[k] - d->dlam[k]);
      d->dlam[k] = d->dz[k];
      d->dz[k] = t;
    }
    // Finish calculation in dlam_lbz, dlam_ubz
    for (k=0; k<p->nz; ++k) {
      d->dlam_lbz[k] -= d->lam_lbz[k] * d->dz[k];
      d->dlam_lbz[k] *= d->dinv_lbz[k];
    }
    for (k=0; k<p->nz; ++k) {
      d->dlam_ubz[k] += d->lam_ubz[k] * d->dz[k];
      d->dlam_ubz[k] *= d->dinv_ubz[k];
    }
    // Finish calculation of dlam(x)
    for (k=0; k<p->nx; ++k) d->dlam[k] += d->dlam_ubz[k] - d->dlam_lbz[k];
    // Maximum primal and dual step
    qp_stepsize(d, &alpha, &alpha, 1.);
    // Calculate sigma
    sigma = qp_calc_sigma(d, alpha);
    // Prepare corrector step
    qp_corrector_prepare(d, sigma * d->mu);
    // Solve to get step
    d->linsys = d->rz;
  }

  template<typename T1>
  void qp_ipstep(casadi_qpip_data<T1>* d, T1 alpha_pr, T1 alpha_du) {
    // Local variables
    casadi_int k;
    const casadi_qp_prob<T1>* p = d->prob;
    // Primal step
    for (k=0; k<p->nz; ++k) d->z[k] += alpha_pr * d->dz[k];
    // Dual step
    for (k=0; k<p->nz; ++k) d->lam[k] += alpha_du * d->dlam[k];
    for (k=0; k<p->nz; ++k) d->lam_lbz[k] += alpha_du * d->dlam_lbz[k];
    for (k=0; k<p->nz; ++k) d->lam_ubz[k] += alpha_du * d->dlam_ubz[k];
  }

  template<typename T1>
  void qp_stepsize(casadi_qpip_data<T1>* d, T1* alpha_pr, T1* alpha_du,
      T1 tau) {
    // Local variables
    casadi_int k;
    const casadi_qp_prob<T1>* p = d->prob;
    // Primal step
    *alpha_pr = 1. / tau;
    for (k=0; k<p->nz; ++k) {
      if (d->dz[k] > 0 && d->ubz[k] < p->inf) {
        *alpha_pr = fmin(*alpha_pr, (d->ubz[k] - d->z[k]) / d->dz[k]);
      } else if (d->dz[k] < 0 && d->lbz[k] > -p->inf) {
        *alpha_pr = fmin(*alpha_pr, (d->lbz[k] - d->z[k]) / d->dz[k]);
      }
    }
    *alpha_pr *= tau;
    // Dual step
    *alpha_du = 1. / tau;
    for (k=0; k<p->nz; ++k) {
      if (d->dlam_lbz[k] < 0. && d->lam_lbz[k] > 0.) {
        *alpha_du = fmin(*alpha_du, -d->lam_lbz[k] / d->dlam_lbz[k]);
      } else if (d->dlam_ubz[k] < 0. && d->lam_ubz[k] > 0.) {
        *alpha_du = fmin(*alpha_du, -d->lam_ubz[k] / d->dlam_ubz[k]);
      }
    }
    *alpha_du *= tau;
  }

  template<typename T1>
  T1 qp_calc_sigma(casadi_qpip_data<T1>* d, T1 alpha) {
    // Local variables
    T1 sigma;
    casadi_int k;
    const casadi_qp_prob<T1>* p = d->prob;
    // Calculate projected mu (and save to sigma variable)
    sigma = 0;
    for (k = 0; k < p->nz; ++k) {
      // Lower bound
      if (d->lbz[k] > -p->inf && d->ubz[k] > d->lbz[k] + p->dmin) {
        sigma += (d->lam_lbz[k] + alpha * d->dlam_lbz[k])
          * (d->z[k] + alpha * d->dz[k] - d->lbz[k]);
      }
      // Upper bound
      if (d->ubz[k] < p->inf && d->ubz[k] > d->lbz[k] + p->dmin) {
        sigma += (d->lam_ubz[k] + alpha * d->dlam_ubz[k])
          * (d->ubz[k] + alpha * d->dz[k] - d->z[k]);
      }
    }
    // Divide mu by total number of finite constraints
    if (d->n_con > 0) sigma /= d->n_con;
    // Finish calculation of sigma := (mu_aff / mu)^3
    sigma /= d->mu;
    sigma *= sigma * sigma;
    return sigma;
  }

  template<typename T1>
  void qp_corrector_prepare(casadi_qpip_data<T1>* d, T1 shift) {
    // Local variables
    casadi_int k;
    const casadi_qp_prob<T1>* p = d->prob;
    // Modified residual in lam_lbz, lam_ubz
    for (k=0; k<p->nz; ++k) d->rlam_lbz[k] = d->dlam_lbz[k] * d->dz[k] - shift;
    for (k=0; k<p->nz; ++k) d->rlam_ubz[k] = d->dlam_ubz[k] * d->dz[k] - shift;
    // Difference in tilde(r)_x, tilde(r)_lamg
    for (k=0; k<p->nz; ++k)
      d->rz[k] = d->dinv_lbz[k] * d->rlam_lbz[k]
        - d->dinv_ubz[k] * d->rlam_ubz[k];
    // Difference in tilde(r)_g
    for (k=p->nx; k<p->nz; ++k) {
      d->rlam[k] = d->rz[k];
      d->rz[k] *= d->D[k] / (d->S[k] * d->S[k]);
    }
    // Scale and negate right-hand-side
    for (k=0; k<p->nz; ++k) d->rz[k] *= -d->S[k];
  }

  template<typename T1>
  void qp_corrector(casadi_qpip_data<T1>* d) {
    // Local variables
    T1 t, alpha;
    casadi_int k;
    const casadi_qp_prob<T1>* p = d->prob;
    // Scale results
    for (k=0; k<p->nz; ++k) d->rz[k] *= d->S[k];
    // Calculate step in z(g), lam(g)
    for (k=p->nx; k<p->nz; ++k) {
      t = d->D[k] / (d->S[k] * d->S[k]) * (d->rz[k] - d->rlam[k]);
      d->rlam[k] = d->rz[k];
      d->rz[k] = t;
    }
    // Update step in dz
    for (k=0; k<p->nz; ++k) d->dz[k] += d->rz[k];
    // Update step in lam_lbz
    for (k=0; k<p->nz; ++k) {
      t = d->dinv_lbz[k] * (d->rlam_lbz[k] + d->lam_lbz[k] * d->rz[k]);
      d->dlam_lbz[k] -= t;
      if (k<p->nx) d->dlam[k] -= t;
    }
    // Update step in lam_ubz
    for (k=0; k<p->nz; ++k) {
      t = d->dinv_ubz[k] * (d->rlam_ubz[k] - d->lam_ubz[k] * d->rz[k]);
      d->dlam_ubz[k] -= t;
      if (k<p->nx) d->dlam[k] += t;
    }
    // Maximum primal and dual step
    qp_stepsize(d, &alpha, &alpha, .9);
    // Take step
    qp_ipstep(d, alpha, alpha);
  }

  template<typename T1>
  int qp_ip(casadi_qpip_data<T1>* d) {
    // Local variables
    const casadi_qp_prob<T1>* p = d->prob;
    switch (d->next) {
      case QP_RESIDUAL:
        // Calculate residual
        calc_res(d);
        d->task = QP_PROGRESS;
        d->next = QP_NEWITER;
        return 1;
      case QP_NEWITER:
        // New iteration
        if (qp_ip_iter(d)) return 0;
        d->task = QP_FACTOR;
        d->next = QP_PREPARE;
        return 1;
      case QP_PREPARE:
        // Prepare predictor step
        qp_predictor_prepare(d);
        d->task = QP_SOLVE;
        d->next = QP_PREDICTOR;
        return 1;
      case QP_PREDICTOR:
        // Complete predictor step
        qp_predictor(d);
        d->task = QP_SOLVE;
        d->next = QP_CORRECTOR;
        return 1;
      case QP_CORRECTOR:
        // Complete predictor step
        qp_corrector(d);
        d->task = QP_MV;
        d->next = QP_RESIDUAL;
        return 1;
      default:
        return 0;
    }
    // Error
    return 0;
  }

  Qpchasm::Qpchasm(const std::string& name, const std::map<std::string, Sparsity> &st)
    : Conic(name, st) {
  }

  Qpchasm::~Qpchasm() {
    clear_mem();
  }

  const Options Qpchasm::options_
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

  void Qpchasm::init(const Dict& opts) {
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

    // D_x, D_g
    sz_w += p_.nz;
    // d.S
    sz_w += p_.nz;
    // lam_lbz, lam_ubz
    sz_w += p_.nz;
    sz_w += p_.nz;
    // dlam_lbz, dlam_ubz
    sz_w += p_.nz;
    sz_w += p_.nz;
    // Residual
    sz_w += p_.nz;
    sz_w += p_.nz;
    sz_w += p_.nz;
    sz_w += p_.nz;
    // Inverse of distance to bounds
    sz_w += p_.nz;
    sz_w += p_.nz;

    alloc_iw(sz_iw, true);
    alloc_w(sz_w, true);

    if (print_header_) {
      // Print summary
      print("-------------------------------------------\n");
      print("This is casadi::QPCHASM\n");
      print("Number of variables:                       %9d\n", nx_);
      print("Number of constraints:                     %9d\n", na_);
      print("Number of nonzeros in H:                   %9d\n", H_.nnz());
      print("Number of nonzeros in A:                   %9d\n", A_.nnz());
      print("Number of nonzeros in KKT:                 %9d\n", kkt_.nnz());
      print("Number of nonzeros in QR(V):               %9d\n", sp_v_.nnz());
      print("Number of nonzeros in QR(R):               %9d\n", sp_r_.nnz());
    }
  }

  void Qpchasm::set_qp_prob() {
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

  int Qpchasm::init_mem(void* mem) const {
    if (Conic::init_mem(mem)) return 1;
    auto m = static_cast<QpchasmMemory*>(mem);
    m->return_status = "";
    return 0;
  }

  int Qpchasm::
  solve(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<QpchasmMemory*>(mem);
    // Message buffer
    char buf[121];
    // Setup data structure
    casadi_qpip_data<double> d;
    d.prob = &p_;
    d.nz_h = arg[CONIC_H];
    d.g = arg[CONIC_G];
    d.nz_a = arg[CONIC_A];
    // D_x, D_a
    d.D = w; w += p_.nz;
    // d.S
    d.S = w; w += p_.nz;
    // lam_lbx, lam_ubz
    d.lam_lbz = w; w += p_.nz;
    d.lam_ubz = w; w += p_.nz;
    // dlam_lbx, dlam_ubz
    d.dlam_lbz = w; w += p_.nz;
    d.dlam_ubz = w; w += p_.nz;
    // Residual
    d.rz = w; w += p_.nz;
    d.rlam = w; w += p_.nz;
    d.rlam_lbz = w; w += p_.nz;
    d.rlam_ubz = w; w += p_.nz;
    // Inverse of distance to bounds
    d.dinv_lbz = w; w += p_.nz;
    d.dinv_ubz = w; w += p_.nz;

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
    casadi_fill(d.lam_lbz, p_.nz, 0.);
    casadi_fill(d.lam_ubz, p_.nz, 0.);

    // Find interior point
    init_ip(&d);
    uout() << d.n_con << " finite constraints\n";

    // Transpose A
    casadi_trans(d.nz_a, p_.sp_a, d.nz_at, p_.sp_at, d.iw);

    print_vec("init z", d.z, p_.nz);
    print_vec("init lam", d.lam, p_.nz);
    print_vec("init lam_lbz", d.lam_lbz, p_.nz);
    print_vec("init lam_ubz", d.lam_ubz, p_.nz);

    // Reverse communication loop
    while (qp_ip(&d)) {
      switch (d.task) {
      case QP_MV:
        // Matrix-vector multiplication
        casadi_clear(d.rz, p_.nz);
        casadi_mv(d.nz_h, p_.sp_h, d.z, d.rz, 0);
        casadi_mv(d.nz_a, p_.sp_a, d.lam + p_.nx, d.rz, 1);
        casadi_mv(d.nz_a, p_.sp_a, d.z, d.rz + p_.nx, 0);
        break;
      case QP_PROGRESS:
        // Print progress
        if (print_iter_) {
          if (d.iter % 10 == 0) {
            // Print header
            if (casadi_qp_print_header(&d, buf, sizeof(buf))) break;
            uout() << buf << "\n";
          }
          // Print iteration
          d.sing = 0;
          d.f = nan;
          d.mina = nan;
          d.imina = 0;
          d.tau = nan;
          d.msg = 0;
          if (casadi_qp_print_iteration(&d, buf, sizeof(buf))) break;
          uout() << buf << "\n";
        }
        break;
      case QP_FACTOR:
        // Factorize KKT
        qp_factorize(&d);
        break;
      case QP_SOLVE:
        // Calculate step
        casadi_qr_solve(d.linsys, 1, 1, p_.sp_v, d.nz_v, p_.sp_r, d.nz_r,
          d.beta, p_.prinv, p_.pc, d.w);
        break;
      }
    }

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
      // Make an iteration
      flag = flag || casadi_qp_iterate(&d);
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

  void Qpchasm::codegen_body(CodeGenerator& g) const {
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

  Dict Qpchasm::get_stats(void* mem) const {
    Dict stats = Conic::get_stats(mem);
    auto m = static_cast<QpchasmMemory*>(mem);
    stats["return_status"] = m->return_status;
    return stats;
  }

  Qpchasm::Qpchasm(DeserializingStream& s) : Conic(s) {
    s.version("Qpchasm", 1);
    s.unpack("Qpchasm::AT", AT_);
    s.unpack("Qpchasm::kkt", kkt_);
    s.unpack("Qpchasm::sp_v", sp_v_);
    s.unpack("Qpchasm::sp_r", sp_r_);
    s.unpack("Qpchasm::prinv", prinv_);
    s.unpack("Qpchasm::pc", pc_);
    s.unpack("Qpchasm::print_iter", print_iter_);
    s.unpack("Qpchasm::print_header", print_header_);
    s.unpack("Qpchasm::print_info", print_info_);
    s.unpack("Qpchasm::print_lincomb_", print_lincomb_);
    set_qp_prob();
    s.unpack("Qpchasm::max_iter", p_.max_iter);
    s.unpack("Qpchasm::min_lam", p_.min_lam);
    s.unpack("Qpchasm::constr_viol_tol", p_.constr_viol_tol);
    s.unpack("Qpchasm::dual_inf_tol", p_.dual_inf_tol);
  }

  void Qpchasm::serialize_body(SerializingStream &s) const {
    Conic::serialize_body(s);

    s.version("Qpchasm", 1);
    s.pack("Qpchasm::AT", AT_);
    s.pack("Qpchasm::kkt", kkt_);
    s.pack("Qpchasm::sp_v", sp_v_);
    s.pack("Qpchasm::sp_r", sp_r_);
    s.pack("Qpchasm::prinv", prinv_);
    s.pack("Qpchasm::pc", pc_);
    s.pack("Qpchasm::print_iter", print_iter_);
    s.pack("Qpchasm::print_header", print_header_);
    s.pack("Qpchasm::print_info", print_info_);
    s.pack("Qpchasm::print_lincomb_", print_lincomb_);
    s.pack("Qpchasm::max_iter", p_.max_iter);
    s.pack("Qpchasm::min_lam", p_.min_lam);
    s.pack("Qpchasm::constr_viol_tol", p_.constr_viol_tol);
    s.pack("Qpchasm::dual_inf_tol", p_.dual_inf_tol);
  }

} // namespace casadi
