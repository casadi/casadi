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


#include "conic_as.hpp"
#include "casadi/core/nlpsol.hpp"

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_CONIC_AS_EXPORT
  casadi_register_conic_as(Conic::Plugin* plugin) {
    plugin->creator = ConicAs::creator;
    plugin->name = "as";
    plugin->doc = ConicAs::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &ConicAs::options_;
    return 0;
  }

  extern "C"
  void CASADI_CONIC_AS_EXPORT casadi_load_conic_as() {
    Conic::registerPlugin(casadi_register_conic_as);
  }

  ConicAs::ConicAs(const std::string& name, const std::map<std::string, Sparsity> &st)
    : Conic(name, st) {
  }

  ConicAs::~ConicAs() {
    clear_mem();
  }

  Options ConicAs::options_
  = {{&Conic::options_},
     {{"nlpsol",
       {OT_STRING,
        "Name of solver."}},
      {"max_iter",
       {OT_INT,
        "Maximum number of iterations [1000]."}}
     }
  };

  void ConicAs::init(const Dict& opts) {
    // Initialize the base classes
    Conic::init(opts);

    // Default options
    string nlpsol_plugin = "ipopt";
    max_iter_ = 1000;

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="nlpsol") {
        nlpsol_plugin = op.second.to_string();
      } else if (op.first=="max_iter") {
        max_iter_ = op.second;
      }
    }

    // Assemble KKT system sparsity
    kkt_ = Sparsity::kkt(H_, A_, false);

    // Transpose of the Jacobian
    AT_ = A_.T();

    // KKT with diagonal
    kktd_ = kkt_ + Sparsity::diag(nx_ + na_);

    // Symbolic QR factorization
    kktd_.qr_sparse(sp_v_, sp_r_, prinv_, pc_);

    // Allocate memory
    alloc_w(kkt_.nnz(), true); // kkt
    alloc_w(kktd_.nnz(), true); // kktd
    alloc_w(nx_, true); // xk
    alloc_w(na_, true); // gk
    alloc_w(nx_, true); // lam_xk
    alloc_w(na_, true); // lam_ak
    alloc_w(A_.nnz()); // trans(A)
    alloc_iw(nx_, na_); // casadi_trans, tau type
    alloc_w(nx_ + na_); // casadi_project, [alpha_x, -lambda_g], [lambda_x, alpha_g], tau memory
    alloc_w(nx_, true); // alpha_x
    alloc_w(na_, true); // alpha_a
    alloc_w(nx_+na_, true); // step
    alloc_w(nx_, true); // dlam_x
    alloc_w(na_, true); // dg

    // Memory for numerical solution
    alloc_w(sp_v_.nnz(), true); // v
    alloc_w(sp_r_.nnz(), true); // r
    alloc_w(nx_+na_, true); // beta
    alloc_w(2*na_+2*nx_); // casadi_qr

    // Print summary
    print("-------------------------------------------\n");
    print("This is casadi::ConicAs.\n");
    print("Number of variables:                       %9d\n", nx_);
    print("Number of constraints:                     %9d\n", na_);
    print("Work in progress!\n");
  }

  template<typename T1>
  void casadi_set_sub(const T1* y, T1* x, const casadi_int* sp_x,
                      casadi_int rbeg, casadi_int rend,
                      casadi_int cbeg, casadi_int cend) {
    // Local variables
    casadi_int r, c, k;
    // Get sparsities
    casadi_int nrow=sp_x[0], ncol=sp_x[1];
    const casadi_int *colind=sp_x+2, *row=sp_x+2+ncol+1;
    // Set elements in subblock
    for (c=cbeg; c<cend; ++c) {
      for (k=colind[c]; k<colind[c+1] && (r=row[k])<rend; ++k) {
        if (r>=rbeg) x[k] = *y++;
      }
    }
  }

  template<typename T1>
  void casadi_fill_sub(T1 y, T1* x, const casadi_int* sp_x,
                      casadi_int rbeg, casadi_int rend,
                      casadi_int cbeg, casadi_int cend) {
    // Local variables
    casadi_int r, c, k;
    // Get sparsities
    casadi_int nrow=sp_x[0], ncol=sp_x[1];
    const casadi_int *colind=sp_x+2, *row=sp_x+2+ncol+1;
    // Set elements in subblock
    for (c=cbeg; c<cend; ++c) {
      for (k=colind[c]; k<colind[c+1] && (r=row[k])<rend; ++k) {
        if (r>=rbeg) x[k] = y;
      }
    }
  }

  template<typename T1>
  void casadi_row_scal(T1* x, const casadi_int* sp_x, const T1* d) {
    // Local variables
    casadi_int c, k;
    // Get sparsities
    casadi_int ncol=sp_x[1];
    const casadi_int *colind=sp_x+2, *row=sp_x+2+ncol+1;
    // Scale entries
    for (c=0; c<ncol; ++c) {
      for (k=colind[c]; k<colind[c+1]; ++k) {
        x[k] *= d[row[k]];
      }
    }
  }

  void print_vector(const double* x, casadi_int n) {
    cout << vector<double>(x, x+n) << endl;
  }

  void print_matrix(const double* x, const casadi_int* sp_x) {
    Sparsity sp = Sparsity::compressed(sp_x);
    vector<double> nz(x, x+sp.nnz());
    DM(sp, nz).print_dense(cout, false);
    cout << endl;
  }

  template<typename T1>
  void casadi_col_scal(T1* x, const casadi_int* sp_x, const T1* d) {
    // Local variables
    casadi_int c, k;
    // Get sparsities
    casadi_int ncol=sp_x[1];
    const casadi_int *colind=sp_x+2;
    // Scale entries
    for (c=0; c<ncol; ++c) {
      for (k=colind[c]; k<colind[c+1]; ++k) {
        x[k] *= d[c];
      }
    }
  }

  template<typename T1>
  void casadi_add_diag(T1* x, const casadi_int* sp_x, const T1* d) {
    // Local variables
    casadi_int c, k;
    // Get sparsities
    casadi_int ncol=sp_x[1];
    const casadi_int *colind=sp_x+2, *row=sp_x+2+ncol+1;
    // Add to diagonal entry
    for (c=0; c<ncol; ++c) {
      for (k=colind[c]; k<colind[c+1]; ++k) {
        if (row[k]==c) {
          x[k] += d[c];
          break;
        }
      }
    }
  }

  int ConicAs::init_mem(void* mem) const {
    auto m = static_cast<ConicAsMemory*>(mem);
    return 0;
  }

  int ConicAs::
  eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    auto m = static_cast<ConicAsMemory*>(mem);

    // Statistics
    for (auto&& s : m->fstats) s.second.reset();

    if (inputs_check_) {
      check_inputs(arg[CONIC_LBX], arg[CONIC_UBX], arg[CONIC_LBA], arg[CONIC_UBA]);
    }

    // Local variables
    casadi_int i;
    double lb, ub, trial, fk;
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
    double *kkt, *kktd, *xk, *lam_xk, *lam_ak, *v, *r, *beta,
           *alpha_x, *alpha_a, *gk, *step, *dlam_x, *dg;
    kkt = w; w += kkt_.nnz();
    kktd = w; w += kktd_.nnz();
    xk = w; w += nx_;
    lam_xk = w; w += nx_;
    lam_ak = w; w += na_;
    v = w; w += sp_v_.nnz();
    r = w; w += sp_r_.nnz();
    beta = w; w += nx_+na_;
    alpha_x = w; w += nx_;
    alpha_a = w; w += na_;
    gk = w; w += nx_;
    step = w; w += nx_+na_;
    dlam_x = w; w += nx_;
    dg = w; w += na_;

    // Pass initial guess
    casadi_copy(x0, nx_, xk);
    casadi_copy(lam_x0, nx_, lam_xk);
    casadi_copy(lam_a0, na_, lam_ak);

    // Copy A' to w
    casadi_trans(a, A_, w, AT_, iw);

    // Assemble the KKT matrix
    casadi_set_sub(h, kkt, kkt_, 0, nx_, 0, nx_); // h
    casadi_set_sub(a, kkt, kkt_, nx_, nx_+na_, 0, nx_); // a
    casadi_set_sub(w, kkt, kkt_, 0, nx_, nx_, nx_+na_); // a'

    // Calculate g
    casadi_fill(gk, na_, 0.);
    casadi_mv(a, A_, xk, gk, 0);

    // Smallest strictly positive number
    const double DMIN = std::numeric_limits<double>::min();

    // Determine initial active set for simple bounds
    for (i=0; i<nx_; ++i) {
      lb = lbx ? lbx[i] : 0.;
      ub = ubx ? ubx[i] : 0.;
      if (lb!=ub) {
        // All inequality constraints are inactive
        lam_xk[i] = 0;
      } else if (xk[i]<=lb) {
        // Lower bound active (including satisfied bounds)
        lam_xk[i] = fmin(lam_xk[i], -DMIN);
      } else {
        // Upper bound active (excluding satisfied bounds)
        lam_xk[i] = fmax(lam_xk[i],  DMIN);
      }
    }

    // Determine initial active set for simple bounds
    for (i=0; i<na_; ++i) {
      lb = lba ? lba[i] : 0.;
      ub = uba ? uba[i] : 0.;
      if (lb!=ub) {
        // All inequality constraints are inactive
        lam_ak[i] = 0;
      } else if (gk[i]<=ub) {
        // Lower bound active (including satisfied bounds)
        lam_ak[i] = fmin(lam_ak[i], -DMIN);
      } else {
        // Upper bound active (excluding satisfied bounds)
        lam_ak[i] = fmax(lam_ak[i],  DMIN);
      }
    }

    // kktd sparsity
    const casadi_int* kkt_colind = kktd_.colind();
    const casadi_int* kkt_row = kktd_.row();

    // QP iterations
    for (casadi_int iter=0; iter<max_iter_; ++iter) {

      // Debugging
      if (verbose_) {
        print("Current xk = \n");
        print_vector(xk, nx_);
        print("Current gk = \n");
        print_vector(gk, na_);
        print("Current lam_xk = \n");
        print_vector(lam_xk, nx_);
        print("Current lam_ak = \n");
        print_vector(lam_ak, na_);
      }

      // Copy kkt to kktd
      casadi_project(kkt, kkt_, kktd, kktd_, w);

      // Loop over kktd entries (left two blocks of the transposed KKT)
      for (casadi_int c=0; c<nx_; ++c) {
        if (lam_xk[c]!=0) {
          // Zero out column, set diagonal entry to 1
          for (casadi_int k=kkt_colind[c]; k<kkt_colind[c+1]; ++k) {
            kktd[k] = kkt_row[k]==c ? 1. : 0.;
          }
        }
      }

      // Loop over kktd entries (right two blocks of the transposed KKT)
      for (casadi_int c=0; c<na_; ++c) {
        if (lam_ak[c]==0) {
          // Zero out column, set diagonal entry to -1
          for (casadi_int k=kkt_colind[nx_+c]; k<kkt_colind[nx_+c+1]; ++k) {
            kktd[k] = kkt_row[k]==nx_+c ? -1. : 0.;
          }
        }
      }

      // QR factorization
      casadi_qr(kktd_, kktd, w, sp_v_, v, sp_r_, r, beta, get_ptr(prinv_), get_ptr(pc_));

      // Evaluate gradient of the Lagrangian and constraint functions
      casadi_copy(g, nx_, step);
      casadi_mv(h, H_, xk, step, 0); // gradient of the objective
      casadi_mv(a, A_, lam_ak, step, 1); // gradient of the Lagrangian
      casadi_copy(gk, na_, step + nx_); // constraint evaluation

      // Correct for active simple bounds
      for (i=0; i<nx_; ++i) {
        if (lam_xk[i]!=0.) {
          step[i] = xk[i];
          if (lbx && lam_xk[i]<0) step[i] -= lbx[i];
          if (ubx && lam_xk[i]>0) step[i] -= ubx[i];
        }
      }

      // Correct for inactive constraints
      for (i=0; i<na_; ++i) {
        if (lam_ak[i]==0) {
          step[nx_+i] = 0.; // -lam_ak[i]
        } else if (lba && lam_ak[i]<0) {
          step[nx_+i] -= lba[i];
        } else if (uba && lam_ak[i]>0) {
          step[nx_+i] -= uba[i];
        }
      }

      if (verbose_) {
        print("KKT residual = \n");
        print_vector(step, nx_ + na_);
      }

      // Negative residual
      casadi_scal(nx_+na_, -1., step);

      // Solve to get primal-dual step
      casadi_qr_solve(step, 1, 1, sp_v_, v, sp_r_, r, beta,
                      get_ptr(prinv_), get_ptr(pc_), w);

      // Calculate change in Lagrangian gradient
      casadi_fill(dlam_x, nx_, 0.);
      casadi_mv(h, H_, step, dlam_x, 0); // gradient of the objective
      casadi_mv(a, A_, step+nx_, dlam_x, 1); // gradient of the Lagrangian

      // Step in lambda_x
      casadi_scal(nx_, -1., dlam_x);

      // Step in g
      casadi_fill(dg, na_, 0.);
      casadi_mv(a, A_, step, dg, 0);

      if (verbose_) {
        print("dx = \n");
        print_vector(step, nx_);
        print("dg = \n");
        print_vector(dg, na_);
        print("dlam_x = \n");
        print_vector(dlam_x, nx_);
        print("dlam_g = \n");
        print_vector(step+nx_, na_);
      }

      // Get maximum step size
      double tau = 1.;

      // Remember best tau for each constraint
      casadi_fill(w, nx_+na_, -1.);

      // iw will be used to mark the new sign:
      // -1: Lower bound became active
      //  0: Bound became inactive
      //  1: Upper bound became active

      // Loop over primal variables
      for (i=0; i<nx_; ++i) {
        lb = lbx ? lbx[i] : 0.;
        ub = ubx ? ubx[i] : 0.;
        if (lam_xk[i]==0.) {
          // Trial step
          trial=xk[i] + tau*step[i];

          // Constraint is inactive, check for primal blocking constraints
          if (trial<=lb && xk[i]>lb) {
            // Lower bound hit
            tau = (lb-xk[i])/step[i];
            w[i] = tau;
            iw[i] = -1;
          } else if (trial>=ub && xk[i]<ub) {
            // Upper bound hit
            tau = (ub-xk[i])/step[i];
            w[i] = tau;
            iw[i] = 1;
          }
        } else {
          trial = lam_xk[i] + tau*dlam_x[i];
          // Constraint is active, check for dual blocking constraints
          if ((lam_xk[i]<0. && trial>=0) || (lam_xk[i]>0. && trial<=0)) {
            // Sign changes
            tau = -lam_xk[i]/dlam_x[i];
            w[i] = tau;
            iw[i] = 0;
          }
        }
      }

      // Loop over constraints
      for (i=0; i<na_; ++i) {
        lb = lba ? lba[i] : 0.;
        ub = uba ? uba[i] : 0.;
        if (lam_ak[i]==0.) {
          // Trial step
          trial=gk[i] + tau*dg[i];
          // Constraint is inactive, check for primal blocking constraints
          if (trial<lb && gk[i]>=lb) {
            // Lower bound hit
            tau = (lb-gk[i])/dg[i];
            w[nx_+i] = tau;
            iw[nx_+i] = -1;
          } else if (trial>ub && gk[i]<=ub) {
            // Upper bound hit
            tau = (ub-gk[i])/dg[i];
            w[nx_+i] = tau;
            iw[nx_+i] = 1;
          }
        } else {
          trial = lam_ak[i] + tau*step[nx_+i];
          // Constraint is active, check for dual blocking constraints
          if ((lam_ak[i]<0. && trial>=0) || (lam_ak[i]>0. && trial<=0)) {
            // Sign changes
            tau = -lam_ak[i]/step[nx_+i];
            w[nx_+i] = tau;
            iw[nx_+i] = 0;
          }
        }
      }

      // Take primal step
      casadi_axpy(nx_, tau, step, xk);

      // Update lam_xk carefully
      for (i=0; i<nx_; ++i) {
        // Get the current sign
        casadi_int s = lam_xk[i]>0. ? 1 : lam_xk[i]<0. ? -1 : 0;
        // Account for sign changes
        if (w[i]==tau) s = iw[i];
        // Take step
        lam_xk[i] += tau*dlam_x[i];
        // Ensure correct sign
        switch (s) {
          case -1: lam_xk[i] = fmin(lam_xk[i], -DMIN); break;
          case  1: lam_xk[i] = fmax(lam_xk[i],  DMIN); break;
          case  0: lam_xk[i] = 0.; break;
        }
      }

      // Update lam_ak carefully
      for (i=0; i<na_; ++i) {
        // Get the current sign
        casadi_int s = lam_ak[i]>0. ? 1 : lam_ak[i]<0. ? -1 : 0;
        // Account for sign changes
        if (w[i]==tau) s = iw[nx_+i];
        // Take step
        lam_ak[i] += tau*step[nx_+i];
        // Ensure correct sign
        switch (s) {
          case -1: lam_ak[i] = fmin(lam_ak[i], -DMIN); break;
          case  1: lam_ak[i] = fmax(lam_ak[i],  DMIN); break;
          case  0: lam_ak[i] = 0.; break;
        }
      }

      // Recalculate g
      casadi_fill(gk, na_, 0.);
      casadi_mv(a, A_, xk, gk, 0);

      // Calculate cost
      fk = casadi_bilin(h, H_, xk, xk)/2. + casadi_dot(nx_, xk, g);

      // Check if there is any active set change
      bool has_change = false;
      for (i=0; i<nx_+na_ && !has_change; ++i) has_change = w[i]==tau;

      // Look for largest x bound violation
      double maxviol = 0.;
      casadi_int imaxviol;
      for (i=0; i<nx_; ++i) {
        lb = lbx ? lbx[i] : 0.;
        ub = ubx ? ubx[i] : 0.;
        if (xk[i] > ub+maxviol) {
          maxviol = xk[i]-ub;
          imaxviol = i;
        } else if (xk[i] < lb-maxviol) {
          maxviol = lb-xk[i];
          imaxviol = i;
        }
      }

      // Look for largest a bound violation
      for (i=0; i<na_; ++i) {
        lb = lba ? lba[i] : 0.;
        ub = uba ? uba[i] : 0.;
        if (gk[i] > ub+maxviol) {
          maxviol = gk[i]-ub;
          imaxviol = nx_+i;
        } else if (gk[i] < lb-maxviol) {
          maxviol = lb-gk[i];
          imaxviol = nx_+i;
        }
      }

      // Print iteration progress:
      print("Iteration %d: fk=%g, tau=%g, |pr|=%g\n", iter, fk, tau, maxviol);

      // Keep iterating?
      if (has_change) continue;

      // Terminate successfully?
      if (maxviol<1e-10) break;

      // Constraint on x or g?
      if (imaxviol<nx_) {
        // No offset
        i = imaxviol;
        // If already active constraint, terminate
        if (lam_xk[i]!=0.) break;
        // Add constraint to active set
        if (xk[i] < lb) {
          lam_xk[i] = -DMIN;
        } else {
          lam_xk[i] = DMIN;
        }
      } else {
        // Remove offset
        i = imaxviol-nx_;

        // If already active constraint, terminate
        if (lam_ak[i]!=0.) break;
        // Add constraint to active set
        if (gk[i] < lb) {
          lam_ak[i] = -DMIN;
        } else {
          lam_ak[i] = DMIN;
        }
      }
    }

    // Calculate optimal cost
    if (f) *f = fk;

    // Get solution
    casadi_copy(xk, nx_, x);
    casadi_copy(lam_xk, nx_, lam_x);
    casadi_copy(lam_ak, na_, lam_a);

    return 0;
  }

} // namespace casadi
