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
    alloc_w(AT_.nnz(), true); // trans_a
    alloc_iw(nx_+na_); // casadi_trans, tau type
    alloc_w(nx_+na_); // casadi_project, tau memory
    alloc_w(nx_+na_, true); // dz
    alloc_w(nx_+na_, true); // dlam
    alloc_w(nx_+na_, true); // kktres
    alloc_w(nx_, true); // glag
    alloc_w(nx_+na_, true); // sens
    alloc_w(nx_, true); // infeas
    alloc_w(nx_, true); // tinfeas

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
    casadi_int ncol=sp_x[1];
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
      print("%lld", x[i]);
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

  void print_matrix(const char* id, const double* x, const casadi_int* sp_x) {
    cout << id << ": ";
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
    double *kkt, *kktd, *z, *lam, *v, *r, *beta, *dz, *dlam, *lbz, *ubz,
           *kktres, *glag, *sens, *trans_a, *infeas, *tinfeas;
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
    trans_a = w; w += AT_.nnz();
    infeas = w; w += nx_;
    tinfeas = w; w += nx_;

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

    // Transpose A
    casadi_trans(a, A_, trans_a, AT_, iw);

    // Assemble the KKT matrix
    casadi_set_sub(h, kkt, kkt_, 0, nx_, 0, nx_); // h
    casadi_set_sub(a, kkt, kkt_, nx_, nx_+na_, 0, nx_); // a
    casadi_set_sub(trans_a, kkt, kkt_, 0, nx_, nx_, nx_+na_); // a'

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

    // R sparsity
    const casadi_int* r_colind = sp_r_.colind();
    //const casadi_int* r_row = sp_r_.row();

    // A sparsity
    //const casadi_int* a_colind = A_.colind();
    //const casadi_int* a_row = A_.row();

    // AT sparsity
    const casadi_int* at_colind = AT_.colind();
    const casadi_int* at_row = AT_.row();

    // Message buffer
    char msg[40] = "";

    // No change so far
    bool changed_active_set = true;

    // Stepsize
    double tau = -1.;

    // Logarithm of the determinant
    double log_det = -1;

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
      //bool prerr_pos;
      for (i=0; i<nx_+na_; ++i) {
        if (z[i] > ubz[i]+prerr) {
          prerr = z[i]-ubz[i];
          iprerr = i;
          //prerr_pos = true;
        } else if (z[i] < lbz[i]-prerr) {
          prerr = lbz[i]-z[i];
          iprerr = i;
          //prerr_pos = false;
        }
      }

      // Calculate dual infeasibility
      double duerr = 0.;
      casadi_int iduerr = -1;
      bool duerr_pos = false; // TODO(jaeandersson): have a look
      for (i=0; i<nx_; ++i) {
        infeas[i] = glag[i]+lam[i];
        double duerr_trial = fabs(infeas[i]);
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

      // Can any constraint be removed without increasing dual infeasibility?
      if (!changed_active_set) {
        double best = duerr;
        casadi_int index = -1;
        for (i=0; i<nx_+na_; ++i) {
          // Skip non-candidates
          if (i==iprerr || lam[i]==0.) continue;
          // Skip equality constraints, no improvement is possible
          if (lbz[i]==ubz[i]) continue;
          // Check largest dual infeasibility resulting from setting lam[i]=0
          double new_duerr;
          if (i<nx_) {
            // Set a lam_x to zero
            new_duerr = fabs(glag[i]);
          } else {
            // Set a lam_a to zero
            new_duerr = 0.;
            for (casadi_int k=at_colind[i-nx_]; k<at_colind[i-nx_+1]; ++k) {
              casadi_int j = at_row[k];
              new_duerr = fmax(new_duerr, fabs((glag[j]-trans_a[k]*lam[i]) + lam[j]));
            }
          }
          // Is this the best one so far?
          if (new_duerr<best) {
            best = new_duerr;
            index = i;
          }
        }
        // Accept, if any
        if (index>=0) {
          lam[index] = 0.;
          changed_active_set = true;
          sprint(msg, sizeof(msg), "Removed redundant %lld", index);
        }
      }

      // Can we improve primal feasibility?
      if (!changed_active_set && iprerr>=0 && lam[iprerr]==0.) {
        // Constraint is free, enforce
        lam[iprerr] = z[iprerr]<lbz[iprerr] ? -DMIN : DMIN;
        changed_active_set = true;
        sprint(msg, sizeof(msg), "Added %lld to reduce |pr|", iprerr);
      }

      // Can we improve dual feasibility by adding a constraint?
      if (!changed_active_set && iduerr>=0) {
        // Calculate the sensitivity of dual infeasibility for iduerr
        casadi_fill(sens, nx_+na_, 0.);
        sens[iduerr] = duerr_pos ? 1. : -1.;
        casadi_mv(a, A_, sens, sens+nx_, 0);

        // Is it possible to improve dual feasibility by adding a constraint
        double best = 0.;
        casadi_int index = -1;
        for (i=0; i<nx_+na_; ++i) {
          // Skip non-candidates
          if (lam[i]!=0. || sens[i]==0.) continue;
          // We cannot enforce an infinite bound
          if (sens[i]>0. ? isinf(ubz[i]) : isinf(lbz[i])) continue;
          // How far are we from the bounds
          double slack = sens[i]>0 ? ubz[i]-z[i] : z[i]-lbz[i];
          // Check if it's the best so far
          if (slack>best) {
            best = slack;
            index = i;
          }
        }

        // Accept, if any
        if (index>=0 && duerr>1e-8) {
          // Enforce
          lam[index] = sens[index]>0 ? DMIN : -DMIN;
          changed_active_set = true;
          sprint(msg, sizeof(msg), "Added %lld to reduce |du|", iprerr);
        } else if (index>=0) {
#if 0
          i = index;
          print("Improvement still possible by enforcing %s bound %lld. "
                "z=%g, lbz=%g, ubz=%g, slack=%g, sensitivity=%g.\n",
                sens[i]>0 ? "upper": "lower", i, z[i], lbz[i], ubz[i], best, sens[i]);
#endif
        }
      }

      if (iter % 10 == 0) {
        // Print header
        print("%10s %15s %15s %6s %15s %6s %10s %10s %15s %40s\n",
              "Iteration", "fk", "|pr|", "con", "|du|", "var", "tau", "lam_min",
              "log|det(KKT)|", "Note");
      }

      // Print iteration progress:
      print("%6d (%1s) %15g %15g %6d %15g %6d %10g %10g %15g %40s\n",
            iter, primal_step ? "P" : "D", fk, prerr, iprerr, duerr, iduerr, tau,
            lam_min, log_det, msg);

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
      msg[0] = '\0';

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

      if (verbose_) {
        print_matrix("KKT", kktd, kktd_);
      }

      // QR factorization
      casadi_qr(kktd_, kktd, w, sp_v_, v, sp_r_, r, beta, get_ptr(prinv_), get_ptr(pc_));

      // Calculate the logarithm of the determinant
      log_det = 0;
      for (casadi_int c=0; c<nx_+na_; ++c) {
        log_det += log(fabs(r[r_colind[c+1]-1]));
      }

      if (verbose_) {
        print_matrix("QR(V)", v, sp_v_);
        print_matrix("QR(R)", r, sp_r_);
        print_vector("beta", beta, nx_+na_);
      }

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

      // Warning if step becomes zero
      if (zero_step) casadi_warning("No search direction");

      // Check primal feasibility in the search direction
      for (i=0; i<nx_+na_ && tau>0.; ++i) {
        double tau1 = tau;
        // Acceptable primal error (must be non-increasing)
        double e = fmax(prerr, 1e-10);
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
        // Consistency check
        casadi_assert(tau<=tau1, "Inconsistent step size calculation");
      }

      // Calculate and order all tau for which there is a sign change
      casadi_fill(w, nx_+na_, 1.);
      casadi_int n_tau = 0;
      for (i=0; i<nx_+na_; ++i) {
        if (dlam[i]==0.) continue; // Skip zero steps
        if (lam[i]==0.) continue; // Skip inactive constraints
        // Skip full steps
        if (lam[i]>0 ? lam[i]>=-dlam[i] : lam[i]<=-dlam[i]) continue;
        // Trial dual step
        double trial_lam = lam[i] + tau*dlam[i];
        if (lam[i]>0 ? trial_lam < -0. : trial_lam > 0.) {
          w[i] = -lam[i]/dlam[i];
        }
        // Where to insert the w[i]
        casadi_int loc;
        for (loc=0; loc<n_tau; ++loc) {
          if (w[i]<w[iw[loc]]) break;
        }
        // Insert element
        n_tau++;
        casadi_int next=i, tmp;
        for (casadi_int j=loc; j<n_tau; ++j) {
          tmp = iw[j];
          iw[j] = next;
          next = tmp;
        }
      }
      // Acceptable dual error (must be non-increasing)
      double e = fmax(duerr, 1e-10);
      /* With the search direction (dz, dlam) and the restriction that when
         lam=0, it stays at zero, we have the following expressions for the
         updated step in the presence of a zero-crossing
             z(tau)   = z(0) + tau*dz
             lam(tau) = lam(0) + tau*dlam     if tau<=tau1
                        0                     if tau>tau1
          where tau*dlam = -lam(0), z(tau) = [x(tau);g(tau)]
          and lam(tau) = [lam_x(tau);lam_g(tau)]
          This gives the following expression for the dual infeasibility
          as a function of tau<=tau1:
            infeas(tau) = g + H*x(tau) + A'*lam_g(tau) + lam_x(tau)
                        = g + H*lam(0) + A'*lam_g(0) + lam_x(0)
                        + tau*H*dz + tau*A'*dlam_g + tau*dlam_x
                        = glag(0) + lam_x(0) + tau*(H*dz + A'*dlam_g + dlam_x)
                        = infeas(0) + tau*tinfeas
            The function is continuous in tau, but tinfeas makes a stepwise
            change when tau=tau1.
          Let us find the largest possible tau, while keeping maximum
          dual infeasibility below e.
        */
      // Tangent of the dual infeasibility at tau=0
      casadi_fill(tinfeas, nx_, 0.);
      casadi_mv(h, H_, dz, tinfeas, 0); // A'*dlam_g + dlam_x==0 by definition
      casadi_mv(a, A_, dlam+nx_, tinfeas, 1);
      casadi_axpy(nx_, 1., dlam, tinfeas);
      // How long step can we take without exceeding e?
      double tau_k = 0.;
      for (casadi_int j=0; j<n_tau; ++j) {
        // Constraint that we're watching
        i = iw[j];
        // Distance to the next tau (may be zero)
        double dtau = w[i] - tau_k;
        // Check if maximum dual infeasibilty gets exceeded
        bool found_tau = false;
        for (casadi_int k=0; k<nx_ && !found_tau; ++k) {
          if (fabs(infeas[k]+dtau*tinfeas[k])>e) {
            double tau1 = fmax(tau_k - dtau*(infeas[k]/tinfeas[k]), 0.);
            if (tau1<tau) {
              // Smallest tau found so far
              found_tau = true;
              tau = tau1;
              index = -1;
              changed_active_set = true;
            }
          }
        }
        // To not allow the active set change if e gets exceeded
        if (found_tau) break;
        // Continue to the next tau
        tau_k = w[i];
        // Update infeasibility
        casadi_axpy(nx_, dtau, tinfeas, infeas);
        // Update the infeasibility tangent for next iteration
        if (i<nx_) {
          // Set a lam_x to zero
          tinfeas[i] -= lam[i];
        } else {
          // Set a lam_a to zero
          for (casadi_int k=at_colind[i-nx_]; k<at_colind[i-nx_+1]; ++k) {
            tinfeas[at_row[k]] -= trans_a[k]*lam[i];
          }
        }
        // Accept the tau, set multiplier to zero or flip sign if equality
        if (i!=index) { // ignore if already taken care of
          changed_active_set = true;
          lam[i] = lbz[i]!=ubz[i] ? 0 : lam[i]<0 ? DMIN : -DMIN;
          sprint(msg, sizeof(msg), "Removed %lld", i);
          dlam[i] = 0.;
        }
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
        if (i==index && s!=sign) {
          sprint(msg, sizeof(msg), "Added %lld (%lld->%lld)", i, s, sign);
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
