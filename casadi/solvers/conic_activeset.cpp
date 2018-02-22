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
      {"pr_tol",
       {OT_DOUBLE,
        "Primal tolerance [1e-8]."}},
      {"du_tol",
       {OT_DOUBLE,
        "Dual tolerance [1e-8]."}}
     }
  };

  void ConicActiveSet::init(const Dict& opts) {
    // Initialize the base classes
    Conic::init(opts);

    // Default options
    max_iter_ = 1000;
    pr_tol_ = 1e-8;
    du_tol_ = 1e-8;

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="max_iter") {
        max_iter_ = op.second;
      } else if (op.first=="pr_tol") {
        pr_tol_ = op.second;
      } else if (op.first=="du_tol") {
        du_tol_ = op.second;
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
    alloc_w(nx_+na_, true); // z=[xk,gk]
    alloc_w(nx_+na_, true); // lbz
    alloc_w(nx_+na_, true); // ubz
    alloc_w(nx_+na_, true); // lam
    alloc_w(A_.nnz()); // trans(A)
    alloc_iw(nx_+na_); // casadi_trans, tau type
    alloc_w(nx_+na_); // casadi_project, tau memory
    alloc_w(nx_+na_, true); // dz
    alloc_w(nx_+na_, true); // dlam
    alloc_iw(nx_+na_, true); // ctype

    // Memory for numerical solution
    alloc_w(sp_v_.nnz(), true); // v
    alloc_w(sp_r_.nnz(), true); // r
    alloc_w(nx_+na_, true); // beta
    alloc_w(2*na_+2*nx_); // casadi_qr

    // Print summary
    print("-------------------------------------------\n");
    print("This is casadi::ConicActiveSet.\n");
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
    casadi_int ncol=sp_x[1];
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

  void ConicActiveSet::
  print_vector(const char* id, const double* x, casadi_int n) const {
    print("%s: [", id);
    for (casadi_int i=0; i<n; ++i) {
      if (i!=0) print(", ");
      print("%g", x[i]);
    }
    print("]\n");
  }

  void ConicActiveSet::
  print_ivector(const char* id, const casadi_int* x, casadi_int n) const {
    print("%s: [", id);
    for (casadi_int i=0; i<n; ++i) {
      if (i!=0) print(", ");
      print("%d", x[i]);
    }
    print("]\n");
  }

  void ConicActiveSet::
  print_signs(const char* id, const double* x, casadi_int n) const {
    print("%s: [", id);
    for (casadi_int i=0; i<n; ++i) {
      print(x[i]==0 ? "0" : x[i]>0 ? "+" : "-");
    }
    print("]\n");
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

    // Local variables
    casadi_int i;
    double trial, fk;
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
    double *kkt, *kktd, *z, *lam, *v, *r, *beta,
           *dz, *dlam, *lbz, *ubz;
    casadi_int* ctype;
    kkt = w; w += kkt_.nnz();
    kktd = w; w += kktd_.nnz();
    z = w; w += nx_+na_;
    lbz = w; w += nx_+na_;
    ubz = w; w += nx_+na_;
    lam = w; w += nx_+na_;
    dz = w; w += nx_+na_;
    dlam = w; w += nx_+na_;
    v = w; w += sp_v_.nnz();
    r = w; w += sp_r_.nnz();
    beta = w; w += nx_+na_;
    ctype = iw; iw += nx_+na_;

    // Smallest strictly positive number
    const double DMIN = std::numeric_limits<double>::min();

    // Bounds on z
    casadi_copy(lbx, nx_, lbz);
    casadi_copy(lba, na_, lbz+nx_);
    casadi_copy(ubx, nx_, ubz);
    casadi_copy(uba, na_, ubz+nx_);

    // Get type of constraints
    enum CType {FREE, LOWER, UPPER, FIXED, RANGE};
    for (i=0; i<nx_+na_; ++i) {
      if (isinf(lbz[i]) && isinf(ubz[i])) {
        ctype[i]=FREE; // unconstrained
      } else if (isinf(ubz[i])) {
        ctype[i]=LOWER; // only lower bound
      } else if (isinf(lbz[i])) {
        ctype[i]=UPPER; // only upper bound
      } else if (lbz[i]==ubz[i]) {
        ctype[i]=FIXED; // equality constraints
      } else {
        ctype[i]=RANGE; // range
      }
    }

    if (verbose_) {
      print_ivector("ctype (x)", ctype, nx_);
      print_ivector("ctype (a)", ctype+nx_, na_);
    }

    // Pass initial guess
    casadi_copy(x0, nx_, z);
    casadi_copy(lam_x0, nx_, lam);
    casadi_copy(lam_a0, na_, lam+nx_);

    // Copy A' to w
    casadi_trans(a, A_, w, AT_, iw);

    // Assemble the KKT matrix
    casadi_set_sub(h, kkt, kkt_, 0, nx_, 0, nx_); // h
    casadi_set_sub(a, kkt, kkt_, nx_, nx_+na_, 0, nx_); // a
    casadi_set_sub(w, kkt, kkt_, 0, nx_, nx_, nx_+na_); // a'

    // Calculate g
    casadi_fill(z+nx_, na_, 0.);
    casadi_mv(a, A_, z, z+nx_, 0);

    // Determine initial active set
    for (i=0; i<nx_+na_; ++i) {
      if (lbz[i]!=ubz[i]) {
        // All inequality constraints are inactive
        lam[i] = 0;
      } else if (z[i]<=lbz[i]) {
        // Lower bound active (including satisfied bounds)
        lam[i] = fmin(lam[i], -DMIN);
      } else {
        // Upper bound active (excluding satisfied bounds)
        lam[i] = fmax(lam[i],  DMIN);
      }
    }

    // kktd sparsity
    const casadi_int* kkt_colind = kktd_.colind();
    const casadi_int* kkt_row = kktd_.row();

    // A sparsity
    const casadi_int* a_colind = A_.colind();
    const casadi_int* a_row = A_.row();

    // No change so far
    bool changed_active_set = true;

    // Feasibility phase or optimality phase?
    bool feasibility_phase = true;

    // Print iteration progress:
    print("Entering feasibility phase\n");

    // QP iterations
    casadi_int iter = 0;
    while (true) {
      // Debugging
      if (verbose_) {
        print_vector("z(x)", z, nx_);
        print_vector("z(g)", z+nx_, na_);
        print_vector("lam(x)", lam, nx_);
        print_vector("lam(a)", lam+nx_, na_);
        print_signs("sign(lam_x)", lam, nx_);
        print_signs("sign(lam_a)", lam+nx_, na_);
      }

      // Recalculate g
      casadi_fill(z+nx_, na_, 0.);
      casadi_mv(a, A_, z, z+nx_, 0);

      // Evaluate gradient of the Lagrangian and constraint functions
      casadi_copy(g, nx_, dz);
      casadi_mv(h, H_, z, dz, 0); // gradient of the objective
      casadi_mv(a, A_, lam+nx_, dz, 1); // gradient of the Lagrangian

      // Recalculate lam_x, without changing the sign
      for (i=0; i<nx_; ++i) {
        if (lam[i]>0) {
          lam[i] = fmax(-dz[i], DMIN);
        } else if (lam[i]<0) {
          lam[i] = fmin(-dz[i], -DMIN);
        }
      }

      // Calculate cost
      fk = casadi_bilin(h, H_, z, z)/2. + casadi_dot(nx_, z, g);

      // Look for largest bound violation
      double maxpr = 0.;
      casadi_int imaxpr;
      for (i=0; i<nx_+na_; ++i) {
        if (z[i] > ubz[i]+maxpr) {
          maxpr = z[i]-ubz[i];
          imaxpr = i;
        } else if (z[i] < lbz[i]-maxpr) {
          maxpr = lbz[i]-z[i];
          imaxpr = i;
        }
      }

      // If calculated residual is positive, we need a negative lhs
      bool negative_lhs = dz[i]+lam[i]>0.;

      // Calculate dual infeasibility
      double maxdu = 0.;
      casadi_int imaxdu;
      for (i=0; i<nx_; ++i) {
        double maxdu_trial = fabs(dz[i]+lam[i]);
        if (maxdu_trial>maxdu) {
          maxdu = maxdu_trial;
          imaxdu = i;
        }
      }

      // Print iteration progress:
      print("Iteration %d: fk=%g, |pr|=%g, |du|=%g\n",
            iter, fk, maxpr, maxdu);

      // Feasibility and optimality?
      bool pr_feasible = maxpr<pr_tol_;
      bool du_feasible = maxdu<du_tol_;

      if (pr_feasible && du_feasible) {
        // Successful return
        break;
      } else if (!changed_active_set) {
        // If active set is unchanged, check feasibility and optimality
        if (feasibility_phase) {
          // Feasibility phase
          if (pr_feasible) {
            // Feasibility achieved
            print("Entering optimality phase\n");
            feasibility_phase = false;
          } else {
            // Restore primal feasibility
            if (imaxpr<nx_) {
              i = imaxpr;
              // Add x constraint
              if (lam[i]!=0.) {
                // If already active constraint, terminate
                casadi_warning("Failed to restore primal feasibility");
                break;
              } else {
                // Add constraint to active set
                if (z[i] < lbz[i]) {
                  lam[i] = fmin(-w[i], -DMIN);
                  changed_active_set = true;
                  continue;
                } else if (z[i] > ubz[i]) {
                  lam[i] = fmax(-w[i],  DMIN);
                  changed_active_set = true;
                  continue;
                } else {
                  casadi_warning("Failed to restore primal feasibility");
                  break; // can it happen?
                }
              }
            } else {
              i = imaxpr-nx_;
              // If already active constraint, terminate
              if (lam[nx_+i]!=0.) {
                casadi_warning("Failed to restore primal feasibility");
                break;
              } else {
                // Add constraint to active set
                if (z[nx_+i] < lbz[nx_+i]) {
                  lam[nx_+i] = -DMIN;
                  changed_active_set = true;
                  continue;
                } else if (z[nx_+i] > ubz[nx_+i]) {
                  lam[nx_+i] = DMIN;
                  changed_active_set = true;
                  continue;
                } else {
                  casadi_warning("Failed to restore primal feasibility");
                  break; // can it happen?
                }
              }
            }
          }
        } else {
          // Optimality phase
          // We're feasible but not optimal, let's remove a redundant constraint
          double best_a = 0.;
          casadi_int ibest_a;

          // If calculated residual is positive, we need a negative lhs
          bool negative_lhs = dz[i]+lam[i]>0.;

          // Check redundancy in x bounds with the right sign
          bool negative_lambda = negative_lhs; // coefficient is 1.
          i=imaxdu;
          if (ctype[i]!=FIXED && lam[i]!=0. && negative_lambda==(lam[i]>0.)) {
            best_a = 1.;
            ibest_a = i;
          }

          // Check redundancy in g bounds matching imaxdu with the right sign
          for (casadi_int k=a_colind[imaxdu]; k<a_colind[imaxdu+1]; ++k) {
            i = a_row[k];
            if (ctype[i]!=FIXED && lam[nx_+i]!=0. && fabs(a[k])>best_a) {
              negative_lambda = negative_lhs==a[k]>0.;
              if (negative_lambda==(lam[nx_+i]>0.)) {
                best_a = fabs(a[k]);
                ibest_a = nx_+i;
              }
            }
          }

          // Remove redundant constraint, if any
          if (best_a>0.) {
            lam[ibest_a] = 0.;
            changed_active_set = true;
            continue;
          }
        }

        casadi_warning("Failed to restore dual feasibility");
        break;
      }

      // Too many iterations?
      if (iter>=max_iter_) {
        casadi_warning("Maximum number of iterations reached");
        break;
      }

      // Start new iteration
      iter++;

      // No change so far
      changed_active_set = false;

      // Calculate KKT residual: Correct for active simple bounds
      for (i=0; i<nx_; ++i) {
        if (lam[i]!=0.) {
          dz[i] = z[i];
          if (ctype[i]==FIXED || lam[i]<0) {
            dz[i] -= lbz[i];
          } else if (lam[i]>0) {
            dz[i] -= ubz[i];
          }
        }
      }

      // Calculate KKT residual: Correct for inactive constraints
      casadi_copy(z+nx_, na_, dz + nx_); // constraint evaluation
      for (i=0; i<na_; ++i) {
        if (lam[nx_+i]==0) {
          dz[nx_+i] = 0.; // -lam[nx_+i]
        } else if (lba && lam[nx_+i]<0) {
          dz[nx_+i] -= lba[i];
        } else if (uba && lam[nx_+i]>0) {
          dz[nx_+i] -= uba[i];
        }
      }

      if (verbose_) {
        print_vector("KKT residual", dz, nx_ + na_);
      }

      // Copy kkt to kktd
      casadi_project(kkt, kkt_, kktd, kktd_, w);

      // Loop over kktd entries (left two blocks of the transposed KKT)
      for (casadi_int c=0; c<nx_; ++c) {
        if (lam[c]!=0) {
          // Zero out column, set diagonal entry to 1
          for (casadi_int k=kkt_colind[c]; k<kkt_colind[c+1]; ++k) {
            kktd[k] = kkt_row[k]==c ? 1. : 0.;
          }
        }
      }

      // Loop over kktd entries (right two blocks of the transposed KKT)
      for (casadi_int c=0; c<na_; ++c) {
        if (lam[nx_+c]==0) {
          // Zero out column, set diagonal entry to -1
          for (casadi_int k=kkt_colind[nx_+c]; k<kkt_colind[nx_+c+1]; ++k) {
            kktd[k] = kkt_row[k]==nx_+c ? -1. : 0.;
          }
        }
      }

      // QR factorization
      casadi_qr(kktd_, kktd, w, sp_v_, v, sp_r_, r, beta, get_ptr(prinv_), get_ptr(pc_));

      // Solve to get primal-dual step
      casadi_scal(nx_+na_, -1., dz);
      casadi_qr_solve(dz, 1, 1, sp_v_, v, sp_r_, r, beta,
                      get_ptr(prinv_), get_ptr(pc_), w);

      // Remove NaNs
      for (i=0; i<nx_+na_; ++i) if (isnan(dz[i])) dz[i]=0.;

      // Calculate change in Lagrangian gradient
      casadi_fill(dlam, nx_, 0.);
      casadi_mv(h, H_, dz, dlam, 0); // gradient of the objective
      casadi_mv(a, A_, dz+nx_, dlam, 1); // gradient of the Lagrangian

      // Step in lam(x)
      casadi_scal(nx_, -1., dlam);

      // Step in lam(g)
      casadi_copy(dz+nx_, na_, dlam+nx_);

      // Step in z(g)
      casadi_fill(dz+nx_, na_, 0.);
      casadi_mv(a, A_, dz, dz+nx_, 0);

      if (verbose_) {
        print_vector("dz(x)", dz, nx_);
        print_vector("dz(g)", dz+nx_, na_);
        print_vector("dlam(x)", dlam, nx_);
        print_vector("dlam(g)", dlam+nx_, na_);
      }

      // If we're in the feasibility phase, make sure that the feasibility
      // improves in the imaxpr direction
      if (!pr_feasible) {
        i = imaxpr;
        if (i<nx_) {
          if (dz[i]==0. || (z[i]<lbz[i])==(dz[i]<0)) {
            casadi_message("Direction does not improve feasibility");
            continue;
          }
        } else {
          i -= nx_;
          if (dz[nx_+i]==0. || (z[nx_+i]<lbz[nx_+i])==(dz[nx_+i]<0)) {
            casadi_message("Direction does not improve feasibility");
            continue;
          }
        }
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
        if (lam[i]==0.) {
          // Trial step
          trial=z[i] + tau*dz[i];
          // Constraint is inactive, check for primal blocking constraints
          if (trial<=lbz[i] && z[i]>lbz[i]) {
            // Lower bound hit
            tau = (lbz[i]-z[i])/dz[i];
            w[i] = tau;
            iw[i] = -1;
          } else if (trial>=ubz[i] && z[i]<ubz[i]) {
            // Upper bound hit
            tau = (ubz[i]-z[i])/dz[i];
            w[i] = tau;
            iw[i] = 1;
          }
        } else if (ctype[i]!=FIXED) {
          trial = lam[i] + tau*dlam[i];
          // Constraint is active, check for dual blocking constraints
          if ((lam[i]<0. && trial>=0) || (lam[i]>0. && trial<=0)) {
            // Sign changes
            tau = -lam[i]/dlam[i];
            w[i] = tau;
            iw[i] = 0;
          }
        }
      }

      // Loop over constraints
      for (i=0; i<na_; ++i) {
        if (lam[nx_+i]==0.) {
          // Trial step
          trial=z[nx_+i] + tau*dz[nx_+i];
          // Constraint is inactive, check for primal blocking constraints
          if (trial<lbz[nx_+i] && z[nx_+i]>=lbz[nx_+i]) {
            // Lower bound hit
            tau = (lbz[nx_+i]-z[nx_+i])/dz[nx_+i];
            w[nx_+i] = tau;
            iw[nx_+i] = -1;
          } else if (trial>ubz[nx_+i] && z[nx_+i]<=ubz[nx_+i]) {
            // Upper bound hit
            tau = (ubz[nx_+i]-z[nx_+i])/dz[nx_+i];
            w[nx_+i] = tau;
            iw[nx_+i] = 1;
          }
        } else if (ctype[nx_+i]!=FIXED) {
          trial = lam[nx_+i] + tau*dlam[nx_+i];
          // Constraint is active, check for sign changes
          if (lam[nx_+i]!=0 && ((lam[nx_+i]>0)!=(trial>0))) {
            // Sign changes
            tau = -lam[nx_+i]/dlam[nx_+i];
            w[nx_+i] = tau;
            iw[nx_+i] = 0;
          }
        }
      }

      if (verbose_) {
        print("tau = %g\n", tau);
      }

      // If tau==0, no step to take
      if (tau==0.) continue;

      // Take primal step
      casadi_axpy(nx_, tau, dz, z);

      // Update lam_x carefully
      for (i=0; i<nx_; ++i) {
        // Get the current sign
        casadi_int s = lam[i]>0. ? 1 : lam[i]<0. ? -1 : 0;
        // Account for sign changes
        if (w[i]==tau) {
          changed_active_set = true;
          s = iw[i];
        }
        // Take step
        lam[i] += tau*dlam[i];
        // Ensure correct sign, unless fixed and nonzero
        if (ctype[i]!=FIXED || lam[i]==0.) {
          switch (s) {
            case -1: lam[i] = fmin(lam[i], -DMIN); break;
            case  1: lam[i] = fmax(lam[i],  DMIN); break;
            case  0: lam[i] = 0.; break;
          }
        }
      }

      // Update lam carefully
      for (i=0; i<na_; ++i) {
        // Get the current sign
        casadi_int s = lam[nx_+i]>0. ? 1 : lam[nx_+i]<0. ? -1 : 0;
        // Account for sign changes
        if (w[nx_+i]==tau) {
          changed_active_set = true;
          s = iw[nx_+i];
        }
        // Take step
        lam[nx_+i] += tau*dlam[nx_+i];
        // Ensure correct sign, unless fixed
        if (ctype[nx_+i]!=FIXED || lam[nx_+i]==0.) {
          switch (s) {
            case -1: lam[nx_+i] = fmin(lam[nx_+i], -DMIN); break;
            case  1: lam[nx_+i] = fmax(lam[nx_+i],  DMIN); break;
            case  0: lam[nx_+i] = 0.; break;
          }
        }
      }
    }

    // Calculate optimal cost
    if (f) *f = fk;

    // Get solution
    casadi_copy(z, nx_, x);
    casadi_copy(lam, nx_, lam_x);
    casadi_copy(lam+nx_, na_, lam_a);

    return 0;
  }

} // namespace casadi
