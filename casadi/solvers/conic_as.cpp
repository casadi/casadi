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
  }

  Options ConicAs::options_
  = {{&Conic::options_},
     {{"nlpsol",
       {OT_STRING,
        "Name of solver."}},
      {"nlpsol_options",
       {OT_DICT,
        "Options to be passed to solver."}}
     }
  };

  void ConicAs::init(const Dict& opts) {
    // Initialize the base classes
    Conic::init(opts);

    // Default options
    string nlpsol_plugin = "ipopt";
    Dict nlpsol_options;

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="nlpsol") {
        nlpsol_plugin = op.second.to_string();
      } else if (op.first=="nlpsol_options") {
        nlpsol_options = op.second;
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
    alloc_iw(na_); // casadi_trans
    alloc_w(nx_ + na_); // casadi_project, [alpha_x, -lambda_g], [lambda_x, alpha_g]
    alloc_w(nx_, true); // alpha_x
    alloc_w(na_, true); // alpha_a

    // Memory for numerical solution
    alloc_w(sp_v_.nnz(), true); // v
    alloc_w(sp_r_.nnz(), true); // r
    alloc_w(nx_+na_, true); // beta
    alloc_w(2*na_+2*nx_); // casadi_qr

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

  int ConicAs::
  eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
    // Local variables
    casadi_int i;
    double lb, ub;
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
    double *kkt, *kktd, *xk, *lam_xk, *lam_ak, *v, *r, *beta, *alpha_x, *alpha_a, *gk;
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

    // Copy kkt to kktd
    casadi_project(kkt, kkt_, kktd, kktd_, w);

    // Make sure that simple bounds are respected
    for (i=0; i<nx_; ++i) {
      lb = lbx ? lbx[i] : 0.;
      ub = ubx ? ubx[i] : 0.;
      xk[i] = std::max(std::min(xk[i], ub), lb);
    }

    // Calculate g
    casadi_fill(gk, na_, 0.);
    casadi_mv(a, A_, xk, gk, 0);

    // Make sure that violated constraints have the correct multiplier signs
    for (i=0; i<na_; ++i) {
      lb = lba ? lba[i] : 0.;
      ub = uba ? uba[i] : 0.;
      if (gk[i]<lb && lam_ak[i]>=0) lam_ak[i] = -std::numeric_limits<double>::min();
      if (gk[i]>ub && lam_ak[i]<=0) lam_ak[i] = std::numeric_limits<double>::min();
    }

    // Calculate alpha_x
    for (i=0; i<nx_; ++i) {
      alpha_x[i] = xk[i];
      if (lam_xk[i]<0 && lbx) {
        alpha_x[i] -= lbx[i];
      } else if (lam_xk[i]>0 && ubx) {
        alpha_x[i] -= ubx[i];
      }
    }

    // Adjust signs for active constraints
//    for (i=0; i<nx_; ++i)


//    casadi_fill(w, nx_, 0.);
  //  casadi_mv(h, H_, xk, w, 0);

    // Calculate alpha_a
    for (i=0; i<na_; ++i) {
      alpha_a[i] = gk[i];
      if (lam_ak[i]<0 && lba) {
        alpha_a[i] -= lba[i];
      } else if (lam_ak[i]>0 && uba) {
        alpha_a[i] -= uba[i];
      }
    }

//    cout << "kktd = " << endl;
  //  print_matrix(kktd, kktd_);
    //cout << "w = " << endl;
    //print_vector(w, nx_ + na_);
    //cout << "endl" << endl;

    // Multiply kktd with diag([alpha_x; -lambda_a]) from the left
    for (i=0; i<nx_; ++i) w[i] = alpha_a[i];
    for (i=0; i<na_; ++i) w[nx_+i] = -lam_ak[i];
    casadi_row_scal(kktd, kktd_, w);

    // Subtract diag([lambda_x; alpha_x])
    for (i=0; i<nx_; ++i) w[i] = -lam_xk[i];
    for (i=0; i<na_; ++i) w[nx_+i] = -alpha_a[i];
    casadi_add_diag(kktd, kktd_, w);

    // QR factorization
    casadi_qr(kktd_, kktd, w, sp_v_, v, sp_r_, r, beta, get_ptr(prinv_), get_ptr(pc_));

    // Calculate residual
    casadi_fill(w, nx_, 0.);
    casadi_mv(h, H_, xk, w, 0);
    for (i=0; i<nx_; ++i) w[i] *= alpha_x[i];
    for (i=0; i<na_; ++i) w[nx_+i] *= -lam_ak[i]*alpha_a[i];

    cout << "residual = " << endl;
    print_vector(w, nx_ + na_);


    cout << "kktd scaled, shifted = " << endl;
    print_matrix(kktd, kktd_);
    cout << "beta" << endl;
    print_vector(beta, nx_ + na_);
    cout << "v = " << endl;
    print_matrix(v, sp_v_);
    cout << "r = " << endl;
    print_matrix(r, sp_r_);

    // Get solution
    casadi_copy(xk, nx_, x);
    casadi_copy(lam_xk, nx_, lam_x);
    casadi_copy(lam_ak, na_, lam_a);

    return 1;
  }

} // namespace casadi
