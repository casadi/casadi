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
      {"tol",
       {OT_DOUBLE,
        "Tolerance [1e-8]."}}
     }
  };

  void ConicActiveSet::init(const Dict& opts) {
    // Initialize the base classes
    Conic::init(opts);

    // Default options
    max_iter_ = 1000;
    tol_ = 1e-8;

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="max_iter") {
        max_iter_ = op.second;
      } else if (op.first=="tol") {
        tol_ = op.second;
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
    alloc_w(nx_+na_, true); // kktres
    alloc_w(nx_, true); // glag
    alloc_w(nx_+na_, true); // sens

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
    casadi_int i, flag;
    double fk;
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
           *dz, *dlam, *lbz, *ubz, *kktres, *glag, *sens;
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
    kktres = w; w += nx_+na_;
    glag = w; w += nx_;
    sens = w; w += nx_+na_;

    // Smallest strictly positive number
    const double DMIN = std::numeric_limits<double>::min();

    // Bounds on z
    casadi_copy(lbx, nx_, lbz);
    casadi_copy(lba, na_, lbz+nx_);
    casadi_copy(ubx, nx_, ubz);
    casadi_copy(uba, na_, ubz+nx_);

    if (verbose_) {
      print_vector("lbz(x)", lbz, nx_);
      print_vector("ubz(x)", ubz, nx_);
      print_vector("lbz(g)", lbz+nx_, na_);
      print_vector("ubz(g)", ubz+nx_, na_);
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
    //const casadi_int* a_colind = A_.colind();
    //const casadi_int* a_row = A_.row();

    // No change so far
    bool changed_active_set = true;

    // Stepsize
    double tau = -1.;

    // Primal and dual error (must be non-increasing)
    double err=inf;

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
      casadi_copy(g, nx_, glag);
      casadi_mv(h, H_, z, glag, 0); // gradient of the objective
      casadi_mv(a, A_, lam+nx_, glag, 1); // gradient of the Lagrangian

      // Recalculate lam(x), without changing the sign
      for (i=0; i<nx_; ++i) {
        if (lam[i]>0) {
          lam[i] = fmax(-glag[i], DMIN);
        } else if (lam[i]<0) {
          lam[i] = fmin(-glag[i], -DMIN);
        }
      }

      // Calculate cost
      fk = casadi_bilin(h, H_, z, z)/2. + casadi_dot(nx_, z, g);

      // Look for largest bound violation
      double prerr = 0.;
      casadi_int iprerr = -1;
      for (i=0; i<nx_+na_; ++i) {
        if (z[i] > ubz[i]+prerr) {
          prerr = z[i]-ubz[i];
          iprerr = i;
        } else if (z[i] < lbz[i]-prerr) {
          prerr = lbz[i]-z[i];
          iprerr = i;
        }
      }

      // Calculate dual infeasibility
      double duerr = 0.;
      casadi_int iduerr = -1;
      bool duerr_pos = false; // TODO(jaeandersson): have a look
      for (i=0; i<nx_; ++i) {
        double duerr_trial = fabs(glag[i]+lam[i]);
        if (duerr_trial>duerr) {
          duerr = duerr_trial;
          iduerr = i;
          duerr_pos = glag[i]+lam[i]>0;
        }
      }

      // Should we try to reduce primal or dual infeasibility?
      bool primal_step = prerr>=duerr;

      // Smallest nonzero lambda
      double lam_min = 0.;
      for (i=0; i<nx_+na_; ++i) {
        if (lam[i]!=0. && fabs(lam[i])<lam_min) lam_min = fabs(lam[i]);
      }

      if (iter % 10 == 0) {
        // Print header
        print("%10s %15s %15s %6s %15s %6s %10s %10s\n",
              "Iteration", "fk", "|pr|", "con", "|du|", "var", "tau", "lam_min");
      }

      // Print iteration progress:
      print("%6d (%1s) %15g %15g %6d %15g %6d %10g %10g\n",
            iter, primal_step ? "P" : "D", fk, prerr, iprerr, duerr, iduerr, tau, lam_min);

      // Overall error (must be non-increasing)
      double old_err = err;
      err = fmax(prerr, duerr);
      casadi_assert(err<=old_err, "Consistency check failed");

      // Check if we need to change the active set
      if (!changed_active_set) {
        // Try to reduce either primal or dual infeasibility, whichever is larger
        if (prerr>=duerr && iprerr>=0 && lam[iprerr]==0.) {
          // Reduce primal infeasibility by adding a constraint
          lam[iprerr] = z[iprerr]>ubz[iprerr] ? DMIN : -DMIN;
          changed_active_set = true;
        } else {
          // Reduce infeasibility by removing a constraint
          if (iduerr>=0) {
            // Recalculate sens for the new iduerr as above
            casadi_fill(sens, nx_+na_, 0.);
            sens[iduerr] = 1.;
            casadi_mv(a, A_, sens, sens+nx_, 0);
            // Look for the best constraint to remove
            double best = 0;
            casadi_int index = -1;
            for (i=0; i<nx_+na_; ++i) {
              // Only enforced constraints are candidates
              if (lam[i]==0.) continue;
              // Projected change from *removing* the constraint
              double trial = -sens[i]*lam[i];
              // if duerr_pos is true, we need a decrease. Pick the largest
              if (duerr_pos ? trial<-best : trial>best) {
                best = fabs(trial);
                index = i;
              }
            }
            // Accept if it decreases infeasibility
            if (index>=0 && fabs(lam[index])<duerr) {
              lam[index] = 0.;
              changed_active_set = true;
            }
          }
        }
      }

      // Successful return if still no change
      if (!changed_active_set) {
        flag = 0;
        break;
      }

      // Too many iterations?
      if (iter>=max_iter_) {
        casadi_warning("Maximum number of iterations reached");
        flag = 1;
        break;
      }

      // Start new iteration
      iter++;

      // No change so far
      changed_active_set = false;

      // KKT residual
      for (i=0; i<nx_+na_; ++i) {
        if (lam[i]>0.) {
          kktres[i] = z[i]-ubz[i];
        } else if (lam[i]<0.) {
          kktres[i] = z[i]-lbz[i];
        } else if (i<nx_) {
          kktres[i] = glag[i];
        } else {
          kktres[i] = -lam[i];
        }
      }

      if (verbose_) {
        print_vector("KKT residual", kktres, nx_ + na_);
      }

      // Calculate the sensitivity of dual infeasibility for iduerr
      casadi_fill(sens, nx_+na_, 0.);
      if (iduerr>=0) {
        casadi_fill(sens, nx_+na_, 0.);
        sens[iduerr] = 1.;
        casadi_mv(a, A_, sens, sens+nx_, 0);
        casadi_scal(nx_+na_, 1./sens[iduerr], sens);
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
      casadi_copy(kktres, nx_+na_, dz);
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

      // For inactive constraints, step is zero
      for (i=0; i<nx_; ++i) if (lam[i]==0.) dlam[i] = 0.;

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

      // Get maximum step size and corresponding index and new sign
      tau = 1.;
      casadi_int sign=0, index=-1;

      // Check if the step is nonzero
      bool zero_step = true;
      for (i=0; i<nx_+na_ && zero_step; ++i) zero_step = dz[i]==0.;
      for (i=0; i<nx_+na_ && zero_step; ++i) zero_step = dlam[i]==0.;
      if (zero_step) tau = 0.;

      // Loop over variables and constraints
      for (i=0; i<nx_+na_ && tau>0.; ++i) {
        double tau1 = tau;
        if (lam[i]==0.) {
          // Numerics checks
          casadi_assert(lam[i] + dlam[i] >= -err, "Numerics failure");
          casadi_assert(lam[i] + dlam[i] <= err, "Numerics failure");
          // Acceptable error (to avoid increasing max error)
          double e = prerr;
          if (dz[i]==0.) continue; // Skip zero steps
          // Check if violation with tau=0 and not improving
          if (dz[i]<0 ? z[i]<=lbz[i]-e : z[i]>=ubz[i]+e) {
            tau = 0.;
            index = i;
            sign = dz[i]<0 ? -1 : 1;
            break;
          }
          // Trial primal step
          double trial_z=z[i] + tau*dz[i];
          if (dz[i]<0 && trial_z<lbz[i]-e) {
            // Trial would increase maximum infeasibility
            tau = (lbz[i]-e-z[i])/dz[i];
            index = i;
            sign = -1;
          } else if (dz[i]>0 && trial_z>ubz[i]+e) {
            // Trial would increase maximum infeasibility
            tau = (ubz[i]+e-z[i])/dz[i];
            index = i;
            sign = 1;
          }
        } else {
          // Check numerics
          casadi_assert(z[i] + dz[i] >= lbz[i] - err, "Numerics failure");
          casadi_assert(z[i] + dz[i] <= ubz[i] + err, "Numerics failure");
          if (dlam[i]==0.) continue; // Skip zero steps
          /*
          if lam[i]<0, we need lam[i] + tau*dlam[i] < e as well as -e < tau*dlam[i] < e
          i.e.
               tau*dlam[i] < e           if dlam[i]>0
          -e < tau*dlam[i] < e - lam[i]  if dlam[i]<0
          if lam[i]>0, we need -e < lam[i] + tau*dlam[i] as well as -e < tau*dlam[i] < e
          i.e.
          -e - lam[i] < tau*dlam[i] < e    if dlam[i]>0
                   -e < tau*dlam[i]        if dlam[i]<0
          */
          // Acceptable error (to avoid increasing max dual error)
          //double e = duerr;
          double e = 0.;
          /*
          if (i<nx_) {
            e = duerr;
          } else {
            double max_a = 0.;
            for (casadi_int c=0; c<nx_; ++c) {
              for (casadi_int k=a_colind[c]; k<a_colind[c+1]; ++k) {
                if (i==nx_+a_row[k]) max_a = fmax(max_a, fabs(a[k]));
              }
            }
            e = duerr/max_a;
          }
          // Trial dual stepsize
          double trial_dlam = fabs(tau*dlam[i]);
          if (trial_dlam>e) {
            tau = e/trial_dlam;
            index = i;
            sign = 0;
          }
          */
          // Trial dual step
          double trial_lam = lam[i] + tau*dlam[i];;
          if (lam[i]>0 && trial_lam < -e) {
            tau = -(lam[i]+e)/dlam[i];
            index = i;
            sign = lbz[i]==ubz[i] ? -1 : 0;
          } else if (lam[i]<0 && trial_lam > e) {
            tau = -(lam[i]-e)/dlam[i];
            index = i;
            sign = lbz[i]==ubz[i] ? 1 : 0;
          }
        }
        // Consistency check
        casadi_assert(tau<=tau1, "Inconsistent step size calculation");
      }

      // Ignore sign changes if they happen for a full step
      if (tau==1.) index = -1;

      if (verbose_) {
        print("tau = %g\n", tau);
      }

      // Take primal step
      casadi_axpy(nx_, tau, dz, z);

      // Update lam carefully
      for (i=0; i<nx_+na_; ++i) {
        // Get the current sign
        casadi_int s = lam[i]>0. ? 1 : lam[i]<0. ? -1 : 0;
        // Account for sign changes
        if (i==index) {
          changed_active_set = true;
          s = sign;
        }
        // Take step
        lam[i] += tau*dlam[i];
        // Ensure correct sign
        switch (s) {
          case -1: lam[i] = fmin(lam[i], -DMIN); break;
          case  1: lam[i] = fmax(lam[i],  DMIN); break;
          case  0: lam[i] = 0.; break;
        }
      }
    }

    // Calculate optimal cost
    if (f) *f = fk;

    // Get solution
    casadi_copy(z, nx_, x);
    casadi_copy(lam, nx_, lam_x);
    casadi_copy(lam+nx_, na_, lam_a);

    return flag;
  }

} // namespace casadi
