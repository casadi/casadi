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
    kkt_ = Sparsity::kkt(H_, A_, false, false);

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
    alloc_w(nx_, true); // glag
    alloc_w(nx_, true); // infeas
    alloc_w(nx_, true); // tinfeas
    alloc_iw(nx_+na_, true); // neverzero
    alloc_iw(nx_+na_, true); // neverupper
    alloc_iw(nx_+na_, true); // neverlower
    alloc_iw(nx_+na_, true); // newsign
    alloc_iw(nx_+na_, true); // backupset
    alloc_iw(nx_+na_); // allzero

    // Memory for numerical solution
    alloc_w(max(sp_v_.nnz()+sp_r_.nnz(), kktd_.nnz()), true); // either v & r or trans(kktd)
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
    vector<double> nz(sp.nnz(), 0.);
    if (x!=0) casadi_copy(x, nz.size(), get_ptr(nz));
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
    casadi_int i, k, flag;
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
           *glag, *trans_a, *infeas, *tinfeas, *vr;
    kkt = w; w += kkt_.nnz();
    kktd = w; w += kktd_.nnz();
    z = w; w += nx_+na_;
    lbz = w; w += nx_+na_;
    ubz = w; w += nx_+na_;
    lam = w; w += nx_+na_;
    dz = w; w += nx_+na_;
    dlam = w; w += nx_+na_;
    vr = w; w += max(sp_v_.nnz()+sp_r_.nnz(), kktd_.nnz());
    v = vr;
    r = vr + sp_v_.nnz();
    beta = w; w += nx_+na_;
    glag = w; w += nx_;
    trans_a = w; w += AT_.nnz();
    infeas = w; w += nx_;
    tinfeas = w; w += nx_;
    casadi_int *neverzero, *neverupper, *neverlower, *newsign, *backupset;
    neverzero = iw; iw += nx_+na_;
    neverupper = iw; iw += nx_+na_;
    neverlower = iw; iw += nx_+na_;
    newsign = iw; iw += nx_+na_;
    backupset = iw; iw += nx_+na_;

    // Smallest strictly positive number
    const double DMIN = std::numeric_limits<double>::min();

    // Bounds on z
    casadi_copy(lbx, nx_, lbz);
    casadi_copy(lba, na_, lbz+nx_);
    casadi_copy(ubx, nx_, ubz);
    casadi_copy(uba, na_, ubz+nx_);

    if (verbose_) {
      print_vector("lbz", lbz, nx_+na_);
      print_vector("ubz", ubz, nx_+na_);
      print_matrix("H", h, H_);
      print_matrix("A", a, A_);
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

    // Look for all-zero rows in kkt
    const casadi_int* kkt_colind = kkt_.colind();
    const casadi_int* kkt_row = kkt_.row();
    for (casadi_int c=0; c<nx_+na_; ++c) iw[c] = 1;
    for (casadi_int c=0; c<nx_+na_; ++c) {
      for (casadi_int k=kkt_colind[c]; k<kkt_colind[c+1]; ++k) {
        if (fabs(kkt[k])>1e-16) iw[kkt_row[k]] = 0;
      }
    }

    // Permitted signs for lam
    for (casadi_int c=0; c<nx_+na_; ++c) {
      neverzero[c] = lbz[c]==ubz[c];
      neverupper[c] = isinf(ubz[c]);
      neverlower[c] = isinf(lbz[c]);
      if (iw[c]) {
        // All-zero row
        if (c<nx_) {
          // Inactive constraint would lead to singular KKT
          neverzero[c] = 1;
        } else {
          // Active constraint would lead to singular KKT
          neverupper[c] = neverlower[c] = 1;
        }
      }
    }

    // Calculate g
    casadi_fill(z+nx_, na_, 0.);
    casadi_mv(a, A_, z, z+nx_, 0);

    // Correct initial active set
    for (i=0; i<nx_+na_; ++i) {
      casadi_assert(!neverzero[i] || !neverupper[i] || !neverlower[i],
                    "No sign possible for " + str(i));
      // Use provided active set when possible
      if (neverzero[i] && lam[i]==0.) {
        lam[i] = neverupper[i] || z[i]-lbz[i] <= ubz[i]-z[i] ? -DMIN : DMIN;
      } else if (neverupper[i] && lam[i]>0.) {
        lam[i] = neverzero[i] ? -DMIN : 0.;
      } else if (neverlower[i] && lam[i]<0.) {
        lam[i] = neverzero[i] ? DMIN : 0.;
      }
    }

    // kktd sparsity
    const casadi_int* kktd_colind = kktd_.colind();
    const casadi_int* kktd_row = kktd_.row();

    // AT sparsity
    const casadi_int* at_colind = AT_.colind();
    const casadi_int* at_row = AT_.row();

    // Message buffer
    char msg[40] = "";

    // No change so far
    bool new_active_set = true;

    // Stepsize
    double tau = 0.;

    // Smallest diagonal value for the QR factorization
    double mina = -1;
    casadi_int imina = -1;

    // Primal error
    double prerr=inf, old_prerr;

    // Singularity in the last iteration
    casadi_int sing = 0; // set to avoid false positive warning

    // Backup active set is available
    bool has_backupset = false;

    // Constraint to be flipped, if any
    casadi_int sign, index=-1;

    // QP iterations
    casadi_int iter = 0;
    while (true) {
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
      old_prerr = prerr;
      prerr = 0.;
      casadi_int iprerr = -1;
      bool prerr_pos = false; // NOTE(jaendersson): suppress used unset warning
      for (i=0; i<nx_+na_; ++i) {
        if (z[i] > ubz[i]+prerr) {
          prerr = z[i]-ubz[i];
          iprerr = i;
          prerr_pos = true;
        } else if (z[i] < lbz[i]-prerr) {
          prerr = lbz[i]-z[i];
          iprerr = i;
          prerr_pos = false;
        }
      }

      // Calculate dual infeasibility
      double duerr = 0.;
      casadi_int iduerr = -1;
      bool duerr_pos = false; // set to avoid false positive warning
      for (i=0; i<nx_; ++i) {
        infeas[i] = glag[i]+lam[i];
        double duerr_trial = fabs(infeas[i]);
        if (duerr_trial>duerr) {
          duerr = duerr_trial;
          iduerr = i;
          duerr_pos = glag[i]+lam[i]>0;
        }
      }

      // Improve primal or dual feasibility
      if (!new_active_set && index<0 && (iprerr>=0 || iduerr>=0)) {
        if (prerr>=duerr) {
          // Try to improve primal feasibility by adding a constraint
          if (lam[iprerr]==0.) {
            // Add the most violating constraint
            index = iprerr;
            sign = z[iprerr]<lbz[iprerr] ? -1 : 1;
            sprint(msg, sizeof(msg), "Added %lld to reduce |pr|", iprerr);
          } else {
            // After a full-step, lam[iprerr] should be zero
            if (prerr < 0.5*old_prerr) {
              // Keep iterating while error is decreasing at a fast-linear rate
              new_active_set = true;
              sprint(msg, sizeof(msg), "|pr| refinement. Rate: %g", prerr/old_prerr);
            }
          }
        } else {
          // Try to improve dual feasibility by removing a constraint
          // We need to increase or decrease infeas[iduerr]. Sensitivity:
          casadi_fill(w, nx_+na_, 0.);
          w[iduerr] = duerr_pos ? -1. : 1.;
          casadi_mv(a, A_, w, w+nx_, 0);
          // Find the best lam[i] to make zero
          casadi_int best_ind = -1;
          double best_w = 0.;
          for (i=0; i<nx_+na_; ++i) {
            // Make sure variable influences duerr
            if (w[i]==0.) continue;
            // Make sure removing the constraint decreases dual infeasibility
            if (w[i]>0. ? lam[i]>=0. : lam[i]<=0.) continue;
            // Maximum infeasibility from setting from setting lam[i]=0
            double new_duerr;
            if (i<nx_) {
              new_duerr = fabs(glag[i]);
            } else {
              new_duerr = 0.;
              for (k=at_colind[i-nx_]; k<at_colind[i-nx_+1]; ++k) {
                new_duerr = fmax(new_duerr, fabs(infeas[at_row[k]]-trans_a[k]*lam[i]));
              }
            }
            // Skip if duerr increases
            if (new_duerr>duerr) continue;
            // Check if best so far
            if (fabs(w[i])>best_w) {
              best_w = fabs(w[i]);
              best_ind = i;
            }
          }
          // Accept, if any
          if (best_ind>=0) {
            index = best_ind;
            sign = 0;
            sprint(msg, sizeof(msg), "Removed %lld to reduce |du|", best_ind);
          }
        }
      }

      // If a constraint was added
      if (index>=0) {
        // Make sure we maintain non-singularity
        if (!sing) {
          // New column that we're trying to add
          casadi_fill(dz, nx_+na_, 0.);
          if (index<nx_ ? sign!=0 : sign==0) {
            dz[index] = 1; // sign is irrelevant
          } else {
            for (k=kkt_colind[index]; k<kkt_colind[index+1]; ++k) {
              dz[kkt_row[k]] = kkt[k];
            }
          }
          // Express it using the other columns
          casadi_qr_solve(dz, 1, 0, sp_v_, v, sp_r_, r, beta,
                          get_ptr(prinv_), get_ptr(pc_), w);
          // If dz[index] is zero, columns are linearly dependent
          if (fabs(dz[index])<1e-12) {
            // Column that we're removing
            casadi_fill(w, nx_+na_, 0.);
            if (index<nx_ ? sign==0 : sign!=0) {
              w[index] = 1; // sign is irrelevant
            } else {
              for (k=kkt_colind[index]; k<kkt_colind[index+1]; ++k) {
                w[kkt_row[k]] = kkt[k];
              }
            }
            // Find best constraint we can flip, if any
            casadi_int best_ind=-1, best_sign;
            double best_slack = -inf;
            for (i=0; i<nx_+na_; ++i) {
              // Can't be the same
              if (i==index) continue;
              // Make sure constraint is flippable
              if (lam[i]==0 ? neverlower[i] && neverupper[i] : neverzero[i]) continue;
              // If dz[i]!=0, column i is redundant
              if (fabs(dz[i])<1e-12) continue;
              // We want to make sure that the new, flipped column i is linearly
              // independent with other columns. We have:
              // flipped_column[index] = dz[0]*column[0] + ... + dz[N-1]*column[N-1]
              // We also require that flipped_column[i] isn't othogonal to
              // (old) column[index], as this will surely lead to singularity
              // This will not cover all cases of singularity, but many important
              // ones. General singularity handling is done below.
              if (i<nx_ ? lam[i]==0. : lam[i]!=0.) {
                // Flipped column is a unit vector
                if (fabs(w[i])<1e-12) continue;
              } else {
                // Flipped column is kkt(:,i)
                double d = 0;
                for (k=kkt_colind[i]; k<kkt_colind[i+1]; ++k) d += kkt[k]*w[kkt_row[k]];
                if (fabs(d)<1e-12) continue;
              }
              // Dual infeasibility
              double new_slack;
              casadi_int new_sign;
              // Check if the best so far
              if (lam[i]==0.) {
                // Which bound is closer?
                new_sign = lbz[i]-z[i] >= z[i]-ubz[i] ? -1 : 1;
                // Better than negative slack, worse than positive slack
                new_slack = 0;
              } else {
                // Slack to the bound
                new_slack = lam[i]>0 ? ubz[i]-z[i] : z[i]-lbz[i];
                new_sign = 0;
              }
              // Discarded?
              casadi_int skip_ind = -1, skip_sign;
              double skip_slack;
              // Best so far?
              if (new_slack > best_slack) {
                if (best_ind>=0) {
                  skip_ind=best_ind;
                  skip_sign=best_sign;
                  skip_slack=best_slack;
                }
                best_slack = new_slack;
                best_ind = i;
                best_sign = new_sign;
              } else {
                skip_ind = i;
                skip_sign = new_sign;
                skip_slack = new_slack;
              }
              // Logic can be improved, issue a warning
              if (skip_ind>=0) {
                print("Note: Discarded %lld to resolve singularity: "
                      "lam=%g, z=%g, lbz=%g, ubz=%g, dz=%g, slack=%g, sign=%lld\n",
                      skip_ind, lam[skip_ind], z[skip_ind],
                      lbz[skip_ind], ubz[skip_ind], dz[skip_ind], skip_slack, skip_sign);
              }
            }

            // Accept, if any
            if (best_ind>=0) {
              lam[best_ind] = best_sign==0 ? 0 : best_sign>0 ? DMIN : -DMIN;
              sprint(msg, sizeof(msg), "%lld->%lld, %lld->%lld",
                     index, sign, best_ind, best_sign);
            } else {
              print("Note: Singularity about to happen\n");
            }
          }
        }

        // Accept active set change
        lam[index] = sign==0 ? 0 : sign>0 ? DMIN : -DMIN;
        new_active_set = true;
        index = -1;
      }

      // Debugging
      if (verbose_) {
        print_vector("z", z, nx_+na_);
        print_vector("lam", lam, nx_+na_);
        print_signs("sign(lam)", lam, nx_+na_);
      }

      // Copy kkt to kktd
      casadi_project(kkt, kkt_, kktd, kktd_, w);

      // Loop over kktd entries (left two blocks of the transposed KKT)
      for (casadi_int c=0; c<nx_; ++c) {
        if (lam[c]!=0) {
          // Zero out column, set diagonal entry to 1
          for (k=kktd_colind[c]; k<kktd_colind[c+1]; ++k) {
            kktd[k] = kktd_row[k]==c ? 1. : 0.;
          }
        }
      }

      // Loop over kktd entries (right two blocks of the transposed KKT)
      for (casadi_int c=0; c<na_; ++c) {
        if (lam[nx_+c]==0) {
          // Zero out column, set diagonal entry to -1
          for (k=kktd_colind[nx_+c]; k<kktd_colind[nx_+c+1]; ++k) {
            kktd[k] = kktd_row[k]==nx_+c ? -1. : 0.;
          }
        }
      }

      if (verbose_) {
        print_matrix("KKT", kktd, kktd_);
      }

      // QR factorization
      casadi_qr(kktd_, kktd, w, sp_v_, v, sp_r_, r, beta, get_ptr(prinv_), get_ptr(pc_));
      if (verbose_) {
        print_matrix("QR(R)", r, sp_r_);
      }

      // Check singularity
      sing = casadi_qr_singular(&mina, &imina, r, sp_r_, get_ptr(pc_), 1e-12);

      // Save or revert to known activeset
      if (!sing) {
        // Remember current active set
        for (i=0; i<nx_+na_; ++i) backupset[i] = lam[i]>0. ? 1 : lam[i]<0 ? -1 : 0;
        has_backupset = true;
      } else if (has_backupset) {
        // Revert to nonsingular active set
        for (i=0; i<nx_+na_; ++i) {
          switch (backupset[i]) {
            case -1: lam[i] = fmin(lam[i], -DMIN); break;
            case  1: lam[i] = fmax(lam[i],  DMIN); break;
            case  0: lam[i] = 0.; break;
          }
        }
        has_backupset = false;
        new_active_set = false;
        sing = 0;
        continue;
      }

      if (iter % 10 == 0) {
        // Print header
        print("%5s %5s %10s %10s %6s %10s %6s %10s %6s %10s %40s\n",
              "Iter", "Null", "fk", "|pr|", "con", "|du|", "var",
              "mindiag(R)", "con", "last tau", "Note");
      }

      // Print iteration progress:
      print("%5d %5d %10.2g %10.2g %6d %10.2g %6d %10.2g %6d %10.2g %40s\n",
            iter, sing, fk, prerr, iprerr, duerr, iduerr,
            mina, imina, tau, msg);

      // Successful return if still no change
      if (!new_active_set) {
        flag = 0;
        break;
      }

      // Break if close enough to optimum
      if (!sing && prerr<1e-12 && duerr<1e-12) {
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
      new_active_set = false;
      sign=0;
      index=-1;

      // Calculate search direction
      if (!sing) {
        // KKT residual
        for (i=0; i<nx_+na_; ++i) {
          if (lam[i]>0.) {
            dz[i] = ubz[i]-z[i];
          } else if (lam[i]<0.) {
            dz[i] = lbz[i]-z[i];
          } else if (i<nx_) {
            dz[i] = -glag[i];
          } else {
            dz[i] = lam[i];
          }
        }

        // Print search direction
        if (verbose_) {
          print_vector("neg kkt residual", dz, nx_+na_);
        }

        // Solve to get primal-dual step
        casadi_qr_solve(dz, 1, 1, sp_v_, v, sp_r_, r, beta,
                        get_ptr(prinv_), get_ptr(pc_), w);
      } else {
        // Get a linear combination of the columns in kktd
        casadi_qr_colcomb(dz, r, sp_r_, get_ptr(pc_), imina, 0);
      }

      // Calculate change in Lagrangian gradient
      casadi_fill(dlam, nx_, 0.);
      casadi_mv(h, H_, dz, dlam, 0); // gradient of the objective
      casadi_mv(a, A_, dz+nx_, dlam, 1); // gradient of the Lagrangian

      // Step in lam(x)
      casadi_scal(nx_, -1., dlam);

      // For inactive constraints, lam(x) step is zero
      for (i=0; i<nx_; ++i) if (lam[i]==0.) dlam[i] = 0.;

      // Step in lam(g)
      casadi_copy(dz+nx_, na_, dlam+nx_);

      // Step in z(g)
      casadi_fill(dz+nx_, na_, 0.);
      casadi_mv(a, A_, dz, dz+nx_, 0);

      // Avoid steps that are nonzero due to numerics
      for (i=0; i<nx_+na_; ++i) if (fabs(dz[i])<1e-14) dz[i] = 0.;

      // Print search direction
      if (verbose_) {
        print_vector("dz", dz, nx_+na_);
        print_vector("dlam", dlam, nx_+na_);
      }

      // Tangent of the dual infeasibility at tau=0
      casadi_fill(tinfeas, nx_, 0.);
      casadi_mv(h, H_, dz, tinfeas, 0); // A'*dlam_g + dlam_x==0 by definition
      casadi_mv(a, A_, dlam+nx_, tinfeas, 1);
      casadi_axpy(nx_, 1., dlam, tinfeas);

      // Handle singularity
      if (sing) {
        // Change in err in the search direction
        double prtau = iprerr<0 ? 0. : prerr_pos ? dz[iprerr]/prerr : -dz[iprerr]/prerr;
        double dutau = iduerr<0 ? 0. : tinfeas[iduerr]/infeas[iduerr];
        double derr = prerr>=duerr ? prtau : dutau;

        // QR factorization of the transpose
        casadi_trans(kktd, kktd_, vr, kktd_, iw);
        casadi_copy(vr, kktd_.nnz(), kktd);
        casadi_qr(kktd_, kktd, w, sp_v_, v, sp_r_, r, beta, get_ptr(prinv_), get_ptr(pc_));

        // Best flip
        double tau_test, best_tau = inf;

        // For all nullspace vectors
        casadi_int nullity_tr, nulli, imina_tr;
        double minat_tr;
        nullity_tr = casadi_qr_singular(&minat_tr, &imina_tr, r, sp_r_,
                                        get_ptr(pc_), 1e-12);
        for (nulli=0; nulli<nullity_tr; ++nulli) {
          // Get a linear combination of the rows in kktd
          casadi_qr_colcomb(w, r, sp_r_, get_ptr(pc_), imina_tr, nulli);
          if (verbose_) {
            print_vector("normal", w, nx_+na_);
          }
          // Look for the best constraint for increasing rank
          for (i=0; i<nx_+na_; ++i) {
            // Check if old column can be removed without decreasing rank
            if (fabs(i<nx_ ? dz[i] : dlam[i])<1e-12) continue;
            // If dot(w, kktd(:,i)-kktd_flipped(:,i))==0, rank won't increase
            double d = i<nx_ ? w[i] : -w[i];
            for (k=kkt_colind[i]; k<kkt_colind[i+1]; ++k) d -= kkt[k]*w[kkt_row[k]];
            if (fabs(d)<1e-12) {
              continue;
            }
            // Is constraint active?
            if (lam[i]==0.) {
              // Make sure that step is nonzero
              if (fabs(dz[i])<1e-12) continue;
              // Step needed to bring z to lower bound
              if (!neverlower[i]) {
                tau_test = (lbz[i]-z[i])/dz[i];
                // Ensure nonincrease in max(prerr, duerr)
                if (!((derr>0. && tau_test>0.) || (derr<0. && tau_test<0.))) {
                  // Only allow removing constraints if tau_test==0
                  if (fabs(tau_test)>=1e-16) {
                    // Check if best so far
                    if (fabs(tau_test)<fabs(best_tau)) {
                      best_tau = tau_test;
                      index = i;
                      sign = -1;
                      sprint(msg, sizeof(msg), "Enforced lbz[%lld] for regularity", i);
                    }
                  }
                }
              }
              // Step needed to bring z to upper bound
              if (!neverupper[i]) {
                tau_test = (ubz[i]-z[i])/dz[i];
                // Ensure nonincrease in max(prerr, duerr)
                if (!((derr>0. && tau_test>0.) || (derr<0. && tau_test<0.))) {
                  // Only allow removing constraints if tau_test==0
                  if (fabs(tau_test)>=1e-16) {
                    // Check if best so far
                    if (fabs(tau_test)<fabs(best_tau)) {
                      best_tau = tau_test;
                      index = i;
                      sign = 1;
                      sprint(msg, sizeof(msg), "Enforced ubz[%lld] for regularity", i);
                    }
                  }
                }
              }
            } else {
              // Make sure that step is nonzero
              if (fabs(dlam[i])<1e-12) continue;
              // Step needed to bring lam to zero
              if (!neverzero[i]) {
                tau_test = -lam[i]/dlam[i];
                // Ensure nonincrease in max(prerr, duerr)
                if ((derr>0. && tau_test>0.) || (derr<0. && tau_test<0.)) continue;
                // Check if best so far
                if (fabs(tau_test)<fabs(best_tau)) {
                  best_tau = tau_test;
                  index = i;
                  sign = 0;
                  sprint(msg, sizeof(msg), "Dropped %s[%lld] for regularity",
                         lam[i]>0 ? "lbz" : "ubz", i);
                }
              }
            }
          }
        }

        // Cannot restore feasibility
        if (index<0) {
          casadi_warning("Cannot restore feasibility");
          flag = 1;
          break;
        }

        // Quick return?
        if (fabs(best_tau)<1e-12) {
          tau = 0.;
          continue;
        }

        // Scale step so that tau=1 is full step
        casadi_scal(nx_+na_, best_tau, dz);
        casadi_scal(nx_+na_, best_tau, dlam);
        casadi_scal(nx_, best_tau, tinfeas);
      }

      // Maximum step size
      tau = 1.;

      // Check if the step is nonzero
      bool zero_step = true;
      for (i=0; i<nx_+na_ && zero_step; ++i) zero_step = dz[i]==0.;
      for (i=0; i<nx_+na_ && zero_step; ++i) zero_step = dlam[i]==0.;
      if (zero_step) {
        tau = 0.;
        continue;
      }

      // Acceptable primal error (must be non-increasing)
      double e = fmax(prerr, 1e-10);

      // Check if violation with tau=0 and not improving
      double dz_max = 0.;
      for (i=0; i<nx_+na_ && tau>0.; ++i) {
        if (-dz[i]>dz_max && z[i]<=lbz[i]-e) {
          index = i;
          sign = -1;
          tau = 0.;
          sprint(msg, sizeof(msg), "Enforcing lbz[%lld]", i);
        } else if (dz[i]>dz_max && z[i]>=ubz[i]+e) {
          index = i;
          sign = 1;
          tau = 0.;
          sprint(msg, sizeof(msg), "Enforcing ubz[%lld]", i);
        }
      }

      // No step needs to be taken if tau==0
      if (tau==0.) continue;

      // Check primal feasibility in the search direction
      for (i=0; i<nx_+na_ && tau>0.; ++i) {
        double tau1 = tau;
        if (dz[i]==0.) continue; // Skip zero steps
        // Trial primal step
        double trial_z=z[i] + tau*dz[i];
        if (dz[i]<0 && trial_z<lbz[i]-e) {
          // Trial would increase maximum infeasibility
          tau = (lbz[i]-e-z[i])/dz[i];
          index = i;
          sign = -1;
          sprint(msg, sizeof(msg), "Enforcing lbz[%lld]", i);
        } else if (dz[i]>0 && trial_z>ubz[i]+e) {
          // Trial would increase maximum infeasibility
          tau = (ubz[i]+e-z[i])/dz[i];
          index = i;
          sign = 1;
          sprint(msg, sizeof(msg), "Enforcing ubz[%lld]", i);
        }
        // Consistency check
        casadi_assert(tau<=tau1, "Inconsistent step size calculation");
      }

      // Calculate and order all tau for which there is a sign change
      casadi_fill(w, nx_+na_, 1.);
      casadi_int n_tau = 0;
      if (tau>0) {
        for (i=0; i<nx_+na_; ++i) {
          if (dlam[i]==0.) continue; // Skip zero steps
          if (lam[i]==0.) continue; // Skip inactive constraints
          // Trial dual step
          double trial_lam = lam[i] + tau*dlam[i];
          // Skip if no sign change
          if (lam[i]>0 ? trial_lam>=0 : trial_lam<=0) continue;
          // Location of the sign change
          w[i] = -lam[i]/dlam[i];
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
      }

      // Get current sign
      for (i=0; i<nx_+na_; ++i) newsign[i] = lam[i]>0. ? 1 : lam[i]<0 ? -1 : 0;

      // Acceptable dual error (must be non-increasing)
      e = fmax(duerr, 1e-10);
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
      // How long step can we take without exceeding e?
      double tau_k = 0.;
      for (casadi_int j=0; j<n_tau; ++j) {
        // Constraint that we're watching
        i = iw[j];
        // Distance to the next tau (may be zero)
        double dtau = w[i] - tau_k;
        // Check if maximum dual infeasibilty gets exceeded
        bool found_tau = false;
        for (k=0; k<nx_ && !found_tau; ++k) {
          if (fabs(infeas[k]+dtau*tinfeas[k])>e) {
            double tau1 = fmax(tau_k - dtau*(infeas[k]/tinfeas[k]), 0.);
            if (tau1<tau) {
              // Smallest tau found so far
              found_tau = true;
              tau = tau1;
              sprint(msg, sizeof(msg), "Dropping %s[%lld]\n",
                     lam[i]>0 ? "ubz" : "lbz", i);
              index = i;
              sign = 0;
            }
          }
        }
        // To not allow the active set change if e gets exceeded
        if (found_tau) break;
        // Continue to the next tau
        tau_k = w[i];
        // Update infeasibility
        casadi_axpy(nx_, dtau, tinfeas, infeas);
        // Update sign or tinfeas
        if (neverzero[i]) {
          // lam changes sign, no change in tinfeas
          newsign[i] = lam[i]<0 ? 1 : -1;
        } else {
          // lam becomes zero, update the infeasibility tangent
          if (i<nx_) {
            // Set a lam_x to zero
            tinfeas[i] -= lam[i];
          } else {
            // Set a lam_a to zero
            for (k=at_colind[i-nx_]; k<at_colind[i-nx_+1]; ++k) {
              tinfeas[at_row[k]] -= trans_a[k]*lam[i];
            }
          }
        }
      }

      if (verbose_) {
        print("tau = %g\n", tau);
      }

      // Take primal step
      casadi_axpy(nx_+na_, tau, dz, z);
      casadi_axpy(nx_+na_, tau, dlam, lam);

      // Update sign
      for (i=0; i<nx_+na_; ++i) {
        switch (newsign[i]) {
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
