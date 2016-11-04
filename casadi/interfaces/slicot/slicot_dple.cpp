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


#include "slicot_dple.hpp"

#include "../../core/std_vector_tools.hpp"
#include "../../core/function/mx_function.hpp"
#include "../../core/function/sx_function.hpp"

#include <cassert>
#include <ctime>
#include <numeric>

// Need an 8-byte integer since libslicot0 is compiled with  fdefault-integer-8
typedef long long int f_int;

extern "C" {
  int mb03vd_(f_int* n, f_int* p, f_int* ilo, f_int* ihi, double *a, f_int* lda1, f_int* lda2,
              double* tau, f_int* ldtau, double* dwork, f_int *info);
  int mb03vy_(f_int* n, f_int* p, f_int* ilo, f_int* ihi, double *a, f_int* lda1, f_int* lda2,
              const double* tau, f_int* ldtau, double* dwork, f_int *ld_work, f_int *info);
  int mb03wd_(char* job, char* compz, f_int* n, f_int* p, f_int* ilo, f_int* ihi, f_int* iloz,
              f_int* ihiz, double *h, f_int* ldh1, f_int* ldh2, double* z, f_int* ldz1,
              f_int* ldz2, double* wr, double *wi, double* dwork, f_int *ld_work, f_int *info);
}

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_DPLE_SLICOT_EXPORT
  casadi_register_dple_slicot(Dple::Plugin* plugin) {
    plugin->creator = SlicotDple::creator;
    plugin->name = "slicot";
    plugin->doc = SlicotDple::meta_doc.c_str();
    plugin->version = 23;
    plugin->exposed.periodic_shur = slicot_periodic_schur;
    return 0;
  }

  extern "C"
  void CASADI_DPLE_SLICOT_EXPORT casadi_load_dple_slicot() {
    Dple::registerPlugin(casadi_register_dple_slicot);
  }

  Options SlicotDple::options_
  = {{&FunctionInternal::options_},
     {{"linear_solver",
       {OT_STRING,
        "User-defined linear solver class. Needed for sensitivities."}},
      {"linear_solver_options",
        {OT_DICT,
         "Options to be passed to the linear solver."}},
      {"psd_num_zero",
        {OT_DOUBLE,
          "Numerical zero used in Periodic Schur decomposition with slicot."
          "This option is needed when your systems has Floquet multipliers"
          "zero or close to zero"}}
     }
  };


  SlicotDple::
  SlicotDple(const std::string& name, const SpDict & st,
                       int nrhs, bool transp) : Dple(name, st, nrhs, transp) {

  }

  SlicotDple::~SlicotDple() {

  }

  void SlicotDple::init(const Dict& opts) {

    Dple::init(opts);

    linear_solver_ = "lapacklu";
    psd_num_zero_ = 1e-12;

    // Read user options
    for (auto&& op : opts) {
      if (op.first=="linear_solver") {
        linear_solver_ = op.second.as_string();
      } else if (op.first=="linear_solver_options") {
        linear_solver_options_ = op.second;
      } else if (op.first=="psd_num_zero") {
        psd_num_zero_ = op.second;
      }
    }

    casadi_assert_message(!pos_def_,
                          "pos_def option set to True: Solver only handles the indefinite case.");
    casadi_assert_message(const_dim_,
                          "const_dim option set to False: Solver only handles the True case.");

    //for (int k=0;k<K_;k++) {
    //  casadi_assert_message(A_[k].isdense(), "Solver requires arguments to be dense.");
    //  casadi_assert_message(V_[k].isdense(), "Solver requires arguments to be dense.");
    //}


    n_ = V_.colind()[1];
    K_ = V_.size1()/n_;

    alloc_w(n_*n_*K_, true); // VZ_
    alloc_w(n_*n_*K_, true); // T_
    alloc_w(n_*n_*K_, true); // Z_
    alloc_w(n_*n_*K_, true); // X_

    alloc_w(n_*n_*K_, true); // Xbar_

    alloc_w(n_*n_*K_, true); // nnKa_
    alloc_w(n_*n_*K_, true); // nnKb_

    alloc_w(n_, true); // eig_real_
    alloc_w(n_, true); // eig_imag_

    alloc_w(2*2*n_*K_, true); // F_
    alloc_w(2*2*K_, true); // FF_


    // There can be at most n partitions
    alloc_iw(n_+1, true); // partition_

    alloc_w(std::max(n_+K_-2, 4*n_)+(n_-1)*K_+2*n_); // dwork_
    alloc_w(4*K_*4, true); // A_
    alloc_w(4*K_, true); // B_
  }


  void SlicotDple::set_work(void* mem, const double**& arg, double**& res,
                                int*& iw, double*& w) const {
    auto m = static_cast<SlicotDpleMemory*>(mem);

    // Set work in base classes
    Dple::set_work(mem, arg, res, iw, w);

    // Lagrange multipliers of the NLP
    m->VZ = w; w += n_*n_*K_;
    m->T = w; w += n_*n_*K_;
    m->Z = w; w += n_*n_*K_;
    m->X = w; w += n_*n_*K_;

    m->Xbar = w; w += n_*n_*K_;
    m->nnKa = w; w += n_*n_*K_;
    m->nnKb = w; w += n_*n_*K_;

    m->eig_real = w; w += n_;
    m->eig_imag = w; w += n_;

    m->F = w; w += 2*2*n_*K_;
    m->FF = w; w += 2*2*K_;

    m->A = w; w += 4*K_*4;
    m->B = w; w += 4*K_;
    m->dwork = w;
    m->partition = iw;

  }


  /** \brief Initalize memory block */
  void SlicotDple::init_memory(void* mem) const {
    Dple::init_memory(mem);
    auto m = static_cast<SlicotDpleMemory*>(mem);

    // Construct linear solvers for low-order Discrete Periodic Sylvester Equations
    // I00X
    // XI00
    // 0XI0
    // 00XI
    //  Special case K=1
    // I+X
    // Solver complexity:  K
    m->dpse_solvers.resize(3);
    for (int i=0;i<3;++i) {
      int np = std::pow(2, i);

      Sparsity sp;
      if (K_==1) {
        sp = Sparsity::dense(np, np);
      } else {
        std::vector<int> row_ind = range(0, np*(np+1)*K_+np+1, np+1);
        std::vector<int> col(np*(np+1)*K_);

        int k = 0;
        for (int l=0;l<np;++l) {
          col[np*(np+1)*k+l*(np+1)] = l;
          for (int m=0;m<np;++m) {
            col[np*(np+1)*k+l*(np+1)+m+1] = (K_-1)*np+m;
          }
        }

        for (k=1;k<K_;++k) {
          for (int l=0;l<np;++l) {
            for (int m=0;m<np;++m) {
              col[np*(np+1)*k+l*(np+1)+m] = (k-1)*np+m;
            }
            col[np*(np+1)*k+l*(np+1)+np] = k*np +l;
          }

        }

        sp = Sparsity(np*K_, np*K_, row_ind, col);
      }

      m->dpse_solvers[i].reserve(n_*(n_+1)/2);
      for (int k=0;k<n_*(n_+1)/2;++k) {
        m->dpse_solvers[i].push_back(Linsol("solver", linear_solver_));
        m->dpse_solvers[i][k].reset(sp);
      }
    }
  }

  /// \cond INTERNAL
  inline int SlicotDple::partindex(const SlicotDpleMemory* m, int i, int j, int k, int r, int c) const {
    return k*n_*n_+(m->partition[i]+r)*n_ + m->partition[j]+c;
  }

  //  A : n-by-l   B: m-by-l
  //  C = A*B'
  void dense_mul_nt(int n, int m, int l, const double *A, const double *B, double *C) {
    for (int i=0;i<n;++i) {
      for (int j=0;j<m;++j) {
        for (int k=0;k<l;++k) {
          C[n*i + j] += A[n*i + k]*B[m*j + k];
        }
      }
    }
  }

  //  A : n-by-l   B: l-by-m
  //  C = A*B
  void dense_mul_nn(int n, int m, int l, const double *A, const double *B, double *C) {
    for (int i=0;i<n;++i) {
      for (int j=0;j<m;++j) {
        for (int k=0;k<l;++k) {
          C[n*i + j] += A[n*i + k]*B[l*k + j];
        }
      }
    }
  }

  //  A : l-by-n   B: l-by-m
  //  C = A'*B
  void dense_mul_tn(int n, int m, int l, const double *A, const double *B, double *C) {
    for (int i=0;i<n;++i) {
      for (int j=0;j<m;++j) {
        for (int k=0;k<l;++k) {
          C[n*i + j] += A[l*k + i]*B[l*k + j];
        }
      }
    }
  }

  /// \endcond

  void SlicotDple::eval(void* mem, const double** arg, double** res, int* iw, double* w) const {
    auto m = static_cast<SlicotDpleMemory*>(mem);

    setup(mem, arg, res, iw, w);

    std::cout << n_ << ":" << K_ << std::endl;
    // Transpose operation (after #554)
    for (int k=0;k<K_;++k) {
      for (int i=0;i<n_;++i) {
        for (int j=0;j<n_;++j) {
          m->X[k*n_*n_+i*n_+j] = arg[DPLE_A][k*n_*n_+n_*j+i];
        }
      }
    }

    slicot_periodic_schur(n_, K_, m->X, m->T, m->Z, m->dwork, m->eig_real, m->eig_imag, psd_num_zero_);

    std::cout << std::vector<double>(m->T, m->T+n_*n_*K_) << std::endl;

    if (error_unstable_) {
      for (int i=0;i<n_;++i) {
        double modulus = sqrt(m->eig_real[i]*m->eig_real[i]+m->eig_imag[i]*m->eig_imag[i]);
        casadi_assert_message(modulus+eps_unstable_ <= 1,
          "SlicotDple: system is unstable."
          "Found an eigenvalue " << m->eig_real[i] << " + " <<
          m->eig_imag[i] << "j, with modulus " << modulus <<
          " (corresponding eps= " << 1-modulus << ")." <<
          std::endl << "Use options and 'error_unstable'"
          "and 'eps_unstable' to influence this message.");
      }
    }

    std::fill(m->A, m->A+4*K_*4, 1);

    int* partition = m->partition;

    // Find a block partition of the T hessenberg form
    partition[0] = 0;int partition_i=1;
    int i = 0;
    int j = 0;
    while (j<n_) {
      while (i<n_ && m->T[i+n_*j]!=0) {
        i+=1;
      }
      j = i;
      partition[partition_i++] = i;
      i += 1;
    }

    // Main loops to loop over blocks of the block-upper triangular A
    // Outer main loop
    for (int l=0;l<partition_i-1;++l) {

      // Inner main loop
      for (int r=0;r<l+1;++r) {

        int n1 = partition[r+1]-partition[r];
        int n2 = partition[l+1]-partition[l];
        int np = n1*n2;

        casadi_assert(n1-1+n2-1>=0);

        Linsol & solver = m->dpse_solvers[n1-1+n2-1][((l+1)*l)/2+r];

        // ********** START ***************
        double * A = m->A;
        double * T = m->T;
        // Special case if K==1
        if (K_==1) {
          for (int ll=0;ll<np;++ll) {
            for (int mm=0;mm<np;++mm) {
              A[ll*np+mm] = -T[partindex(m, r, r, 0, ll/n2, mm/n2)]*T[partindex(m, l, l, 0, ll%n2, mm%n2)];
              if (ll==mm) {
                A[ll*np+mm]+= 1;
              }
            }
          }
        } else {
          std::cout << K_ << np << n1 << n2 << std::endl;
          // Other cases
          int k;
          for (k=0;k<K_-1;++k) {
            for (int ll=0;ll<np;++ll) {
              for (int mm=0;mm<np;++mm) {
                A[np*(np+1)*((k+1)%K_)+ll*(np+1)+mm] =
                    -T[partindex(m, r, r, k, ll/n2, mm/n2)]*T[partindex(m, l, l, k, ll%n2, mm%n2)];
              }
            }
          }
          std::cout << "Apart" << std::vector<double>(A, A+4*4) << std::endl;

          for (int ll=0;ll<np;++ll) {
            for (int mm=0;mm<np;++mm) {
              A[np*(np+1)*((k+1)%K_)+ll*(np+1)+mm+1] =
                -T[partindex(m, r, r, k, ll/n2, mm/n2)]*T[partindex(m, l, l, k, ll%n2, mm%n2)];
            }
          }
        }
        // ********** STOP ***************
        // Solve Discrete Periodic Sylvester Equation Solver

        std::cout << "A" << std::vector<double>(A, A+4*4) << std::endl;

        solver.pivoting(m->A);
        solver.factorize(m->A);

      }
    }


    if (!transp_) {
      for (int d=0;d<nrhs_;++d) {
        std::cout << "nrhs" << nrhs_ << std::endl;

        // ********** START ***************
        // V = blocks([mul([sZ[k].T, V[k], sZ[k]]) for k in range(p)])

        for (int k=0;k<K_;++k) { // K
          double * nnKa = m->nnKa+k*n_*n_;
          double * nnKb = m->nnKb+k*n_*n_;
          // n^2 K

          std::fill(nnKa, nnKa+n_*n_, 0);
          // nnKa[k] <- V[k]*Z[k+1]
          // n^3 K
          dense_mul_nt(n_, n_, n_, arg[DPLE_V]+d*n_*n_*K_ + k*n_*n_, m->Z+((k+1) % K_)*n_*n_, nnKa);
          std::fill(nnKb, nnKb+n_*n_, 0);
          // nnKb[k] <- Z[k+1]'*V[k]*Z[k+1]
          dense_mul_nn(n_, n_, n_, m->Z + ((k+1) % K_)*n_*n_, nnKa, nnKb);
        }

        // ********** STOP ****************

        std::fill(m->X, m->X+n_*n_*K_, 0);

        // Main loops to loop over blocks of the block-upper triangular A
        // Outer main loop
        for (int l=0;l<partition_i-1;++l) { // n
          // F serves an an accumulator for intermediate summation results
          // n^2 K
          std::fill(m->F, m->F+2*2*n_*K_, 0);

          // ********** START ***************
          //for i in range(l):
          //  F[i] = [sum(mul(X[i][j][k], A[l][j][k].T) for j in range(l)) for k in range(p) ]

          for (int i=0;i<l;++i) { // n^2
            int na1 = partition[i+1]-partition[i];
            int nb2 = partition[l+1]-partition[l];
            for (int j=0;j<l;++j) { // n^3
              int na2 = partition[j+1]-partition[j];
              for (int k=0;k<K_;++k) {
                for (int ii=0;ii<na1;++ii) {
                  for (int jj=0;jj<nb2;++jj) {
                    for (int kk=0;kk<na2;++kk) {
                      // n^3 K
                      m->F[k*4*n_+4*i+2*ii+jj] +=
                          m->X[partindex(m, i, j, k, ii, kk)]*m->T[partindex(m, l, j, k, jj, kk)];
                    }
                  }
                }
              }
            }
          }
          // ********** STOP ***************




          // Inner main loop
          for (int r=0;r<l+1;++r) { // n^2

            // ********** START ***************
            // F[r] = [sum(mul(X[r][j][k], A[l][j][k].T) for j in range(l)) for k in range(p) ]
            int na1 = partition[r+1]-partition[r];
            int nb2 = partition[l+1]-partition[l];

            if (r==l) {
              for (int j=0;j<l;++j) { // n^3
                int na2 = partition[j+1]-partition[j];
                for (int k=0;k<K_;++k) { // n^3 K
                  for (int ii=0;ii<na1;++ii) {
                    for (int jj=0;jj<nb2;++jj) {
                      for (int kk=0;kk<na2;++kk) {
                        m->F[k*4*n_+4*r+2*ii+jj] +=
                          m->X[partindex(m, r, j, k, ii, kk)]*m->T[partindex(m, l, j, k, jj, kk)];
                      }
                    }
                  }
                }
              }
            }
            // ********** STOP ***************


            // ********** START ***************
            // FF =   [sum(mul(A[r][i][k], X[i][l][k]) for i in range(r)) for k in range(p)]
            // Each entry of FF is na1-by-na2
            // n^2 K
            std::fill(m->FF, m->FF+2*2*K_, 0);
            {
              int na1 = partition[r+1]-partition[r];
              for (int i=0;i<r;++i) { // n^3
                int nb2 = partition[l+1]-partition[l];
                int na2 = partition[i+1]-partition[i];
                for (int k=0;k<K_;++k) { // n^3 K
                  for (int ii=0;ii<na1;++ii) {
                    for (int jj=0;jj<nb2;++jj) {
                      for (int kk=0;kk<na2;++kk) {
                        m->FF[k*4+2*ii+jj] +=
                            m->T[partindex(m, r, i, k, ii, kk)]*m->X[partindex(m, i, l, k, kk, jj)];
                      }
                    }
                  }
                }
              }
            }


            // ********** STOP ***************

            int n1 = partition[r+1]-partition[r];
            int n2 = partition[l+1]-partition[l];
            int np = n1*n2;

            Linsol & solver = m->dpse_solvers[n1-1+n2-1][((l+1)*l)/2+r];

            // M <- V
            for (int k=0;k<K_;++k) {
              for (int ll=0;ll<n1;++ll) {
                for (int mm=0;mm<n2;++mm) {
                  // n^2 K
                  m->B[np*((k+1)%K_)+ll*n2+mm] =
                    m->nnKb[k*n_*n_+ (partition[r]+ll)*n_ + partition[l]+mm];
                }
              }
            }



            // ********** START ***************
            // M+= [sum(mul(A[r][i][k], F[i][k])  for i in range(r+1)) for k in rang(p)]
            {
              int na1 = partition[r+1]-partition[r];
              int nb2 = partition[l+1]-partition[l];

              for (int i=0;i<r+1;++i) { // n^3
                int na2 = partition[i+1]-partition[i];
                for (int k=0;k<K_;++k) { // n^3 K
                  for (int ii=0;ii<na1;++ii) {
                    for (int jj=0;jj<nb2;++jj) {
                      for (int kk=0;kk<na2;++kk) {
                        m->B[np*((k+1)%K_)+ii*n2+jj] +=
                            m->T[partindex(m, r, i, k, ii, kk)]*m->F[k*4*n_+4*i+2*kk+jj];
                      }
                    }
                  }
                }
              }
            }
            // ********** STOP ***************

            // ********** START ***************
            // M+= [mul(FF[k], A[l][l][k].T) for k in rang(p)]
            {
              int na1 = partition[r+1]-partition[r];
              int na2 = partition[l+1]-partition[l];
              int nb2 = partition[l+1]-partition[l];
              for (int k=0;k<K_;++k) { // n^2 K
                for (int ii=0;ii<na1;++ii) {
                  for (int jj=0;jj<nb2;++jj) {
                    for (int kk=0;kk<na2;++kk) {
                      m->B[np*((k+1)%K_)+ii*n2+jj] += m->FF[k*4+2*ii+kk]*m->T[partindex(m, l, l, k, jj, kk)];
                    }
                  }
                }
              }
            }
            // ********** STOP ***************
            std::cout << "B" << std::vector<double>(m->B, m->B+4) << std::endl;
            // Critical observation: Prepare step is not needed
            // n^2 K
            solver.solve(m->B, 1, true);

            // Extract solution and store it in X
            double * sol = m->B;
            std::cout << std::vector<double>(sol, sol+4) << std::endl;
            // ********** START ***************
            for (int ii=0;ii<partition[r+1]-partition[r];++ii) {
              for (int jj=0;jj<partition[l+1]-partition[l];++jj) {
                for (int k=0;k<K_;++k) { // n^2 K
                  m->X[partindex(m, r, l, k, ii, jj)] = sol[n1*n2*k+n2*ii+jj];
                }
              }
            }

            for (int ii=0;ii<partition[r+1]-partition[r];++ii) {
              for (int jj=0;jj<partition[l+1]-partition[l];++jj) {
                for (int k=0;k<K_;++k) { // n^2 K
                  m->X[partindex(m, l, r, k, jj, ii)] = sol[n1*n2*k+n2*ii+jj];
                }
              }
            }
            // ********** STOP ***************

          }

          // n^3 K
          std::fill(res[DPLE_P]+d*n_*n_*K_, res[DPLE_P]+(d+1)*n_*n_*K_, 0);
        }

        std::fill(m->nnKa, m->nnKa+n_*n_, 0);
        for (int k=0;k<K_;++k) {
          // nnKa[k] <- V[k]*Z[k]'
          // n^3 K
          dense_mul_nn(n_, n_, n_, m->X + k*n_*n_, m->Z+ k*n_*n_, m->nnKa+ k*n_*n_);
          // output <- Z[k]*V[k]*Z[k]'
          dense_mul_tn(n_, n_, n_, m->Z + k*n_*n_, m->nnKa+ k*n_*n_,
                       res[DPLE_P]+d*n_*n_*K_+ k*n_*n_);
        }

      }
    } else { // Transposed

      for (int d=0;d<nrhs_;++d) {


        const double *P_bar = arg[DPLE_V]+d*n_*n_*K_;
        double *Vbar = res[DPLE_P]+d*n_*n_*K_;
        std::fill(Vbar, Vbar+n_*n_*K_, 0);
        std::fill(m->Xbar, m->Xbar+n_*n_*K_, 0);

        // X_bar = [mul([Z[k].T, X_bar[k] , Z[k]]) for k in range(p)]
        for (int k=0;k<K_;++k) {
          double* nnKa = m->nnKa + k*n_*n_;
          std::fill(nnKa, nnKa + n_*n_, 0);

          // n^3 K
          // nnKa[k] <- nnKb*Z[k]
          dense_mul_nt(n_, n_, n_, P_bar+n_*n_*k, m->Z+k*n_*n_, nnKa);
          // Xbar <- Z[k]*V[k]*Z[k]'
          dense_mul_nn(n_, n_, n_, m->Z+k*n_*n_, nnKa, m->Xbar+k*n_*n_);
        }

        // Main loops to loop over blocks of the block-upper triangular A
        // Outer main loop
        for (int l=partition_i-2;l>=0;--l) { // n

          std::fill(m->F, m->F+n_*n_*K_, 0);

          // Inner main loop
          for (int r=l;r>=0;--r) { // n^2

            int n1 = partition[r+1]-partition[r];
            int n2 = partition[l+1]-partition[l];
            int np = n1*n2;

            Linsol & solver = m->dpse_solvers[n1-1+n2-1][((l+1)*l)/2+r];


            // ********** START ***************
            for (int ii=0;ii<partition[r+1]-partition[r];++ii) {
              for (int jj=0;jj<partition[l+1]-partition[l];++jj) {
                for (int k=0;k<K_;++k) {
                  // n^2 K
                  m->B[n1*n2*k+n2*ii+jj] = m->Xbar[partindex(m, r, l, k, ii, jj)];
                }
              }
            }
            if (r!=l) {
              for (int ii=0;ii<partition[r+1]-partition[r];++ii) {
                for (int jj=0;jj<partition[l+1]-partition[l];++jj) {
                  for (int k=0;k<K_;++k) {
                    m->B[n1*n2*k+n2*ii+jj]+= m->Xbar[partindex(m, l, r, k, jj, ii)];
                  }
                }
              }
            }
            // ********** STOP ***************

            // n^2 K
            solver.solve(m->B, 1, false);

            // for k in range(p): V_bar[r][l][k]+=M_bar[k]
            double *Mbar = m->B;

            for (int k=0;k<K_;++k) {
              for (int ll=0;ll<n1;++ll) {
                for (int mm=0;mm<n2;++mm) {
                   // n^2 K
                   Vbar[partindex(m, r, l, k, ll, mm)]+= Mbar[np*((k+1)%K_)+ll*n2+mm];
                }
              }
            }

            // ********** START ***************
            //     for k in range(p):
            //       FF_bar[k] += mul(M_bar[k], A[l][l][k])
            // n^2 K
            std::fill(m->FF, m->FF+n_*n_*K_, 0);
            {
              int na1 = partition[r+1]-partition[r];
              int na2 = partition[l+1]-partition[l];
              int nb2 = partition[l+1]-partition[l];
              for (int k=0;k<K_;++k) {
                for (int ii=0;ii<na1;++ii) {
                  for (int jj=0;jj<nb2;++jj) {
                    for (int kk=0;kk<na2;++kk) {
                      // n^2 K
                      m->FF[k*4+2*ii+kk] +=
                        Mbar[np*((k+1)%K_)+ii*n2+jj]*m->T[partindex(m, l, l, k, jj, kk)];
                    }
                  }
                }
              }
            }

            // ********** STOP ***************


            // for k in range(p):
            //  for i in range(r):
            //    X_bar[i][l][k] +=mul(A[r][i][k].T, FF_bar[k])
            // ********** START ***************
            {
              int na1 = partition[r+1]-partition[r];
              for (int i=0;i<r;++i) { // n^3
                int nb2 = partition[l+1]-partition[l];
                int na2 = partition[i+1]-partition[i];
                for (int k=0;k<K_;++k) {
                  for (int ii=0;ii<na1;++ii) {
                    for (int jj=0;jj<nb2;++jj) {
                      for (int kk=0;kk<na2;++kk) {
                        // n^3 K
                        m->Xbar[partindex(m, i, l, k, kk, jj)] +=
                            m->T[partindex(m, r, i, k, ii, kk)]*m->FF[k*4+2*ii+jj];
                      }
                    }
                  }
                }
              }
            }
            // ********** STOP ***************

            //for k in range(p):
            //  for i in range(r+1):
            //    F_bar[i][k]+= mul(A[r][i][k].T, M_bar[k])
            // ********** START ***************
            {
              int na1 = partition[r+1]-partition[r];
              int nb2 = partition[l+1]-partition[l];

              for (int i=0;i<r+1;++i) { // n^3
                int na2 = partition[i+1]-partition[i];
                for (int k=0;k<K_;++k) {
                  for (int ii=0;ii<na1;++ii) {
                    for (int jj=0;jj<nb2;++jj) {
                      for (int kk=0;kk<na2;++kk) {
                        // n^3 K
                        m->F[k*4*n_+4*i+2*kk+jj] +=
                            m->T[partindex(m, r, i, k, ii, kk)]*Mbar[np*((k+1)%K_)+ii*n2+jj];
                      }
                    }
                  }
                }
              }
            }
            // ********** STOP ***************

            // if r>=l:
            //   for k in range(p):
            //     for j in range(l):
            //       X_bar[r][j][k] +=mul(F_bar[r][k], A[l][j][k])
            // ********** START ***************
            if (r>=l) {
              int na1 = partition[r+1]-partition[r];
              int nb2 = partition[l+1]-partition[l];

              for (int j=0;j<l;++j) { // n^3
                int na2 = partition[j+1]-partition[j];
                for (int k=0;k<K_;++k) {
                  for (int ii=0;ii<na1;++ii) {
                    for (int jj=0;jj<nb2;++jj) {
                      for (int kk=0;kk<na2;++kk) {
                        // n^3 K
                        m->Xbar[partindex(m, r, j, k, ii, kk)] +=
                            m->F[k*4*n_+4*r+2*ii+jj]*m->T[partindex(m, l, j, k, jj, kk)];
                      }
                    }
                  }
                }
              }
            }
            // ********** STOP ***************

          }

          //for i in range(l):
          // for k in range(p):
          //   for j in range(l):
          //     X_bar[i][j][k]+=mul(F_bar[i][k], A[l][j][k])
          // ********** START ***************

          for (int i=0;i<l;++i) { // n^2
            int na1 = partition[i+1]-partition[i];
            int nb2 = partition[l+1]-partition[l];
            for (int j=0;j<l;++j) { // n^3
              int na2 = partition[j+1]-partition[j];
              for (int k=0;k<K_;++k) {
                for (int ii=0;ii<na1;++ii) {
                  for (int jj=0;jj<nb2;++jj) {
                    for (int kk=0;kk<na2;++kk) {
                      // n^3 K
                      m->Xbar[partindex(m, i, j, k, ii, kk)] +=
                          m->F[k*4*n_+4*i+2*ii+jj]*m->T[partindex(m, l, j, k, jj, kk)];
                    }
                  }
                }
              }
            }
          }
          // ********** STOP ***************

        }


        // V_bar = [mul([sZ[k], V_bar[k] , sZ[k].T]) for k in range(p)]
        for (int k=0;k<K_;++k) {
          double* nnKa = m->nnKa+k*n_*n_;
          std::fill(nnKa, nnKa+n_*n_, 0);
          // n^3 K
          dense_mul_nn(n_, n_, n_, Vbar + n_*n_*k, m->Z + ((k+1)%K_)*n_*n_, nnKa);
          std::fill(Vbar + n_*n_*k, Vbar + n_*n_*(k+1), 0);
          dense_mul_tn(n_, n_, n_, m->Z + ((k+1)%K_)*n_*n_, nnKa, Vbar+n_*n_*k);
        }

        // Force symmetry
        for (int k=0;k<K_;++k) {
          for (int r=0;r<n_;++r) {
            for (int l=0;l<r;++l) {
              // n^2 K
              double s = (Vbar[n_*n_*k+r+l*n_]+Vbar[n_*n_*k+r*n_+l])/2;
              Vbar[n_*n_*k+r+l*n_] = s;
              Vbar[n_*n_*k+r*n_+l] = s;
            }
          }
        }

        //std::fill(P_bar.data().begin(), P_bar.data().end(), 0);
      }
    }

  }

  void slicot_mb03vd(int n, int p, int ilo, int ihi, double * a, int lda1, int lda2, double * tau,
                     int ldtau, double * dwork) {
     if (dwork==0) {
       std::vector<double> work = std::vector<double>(n);
       slicot_mb03vd(n, p, ilo, ihi, a, lda1, lda2, tau, ldtau, &work[0]);
       return;
     }
     f_int n_ = n;
     f_int p_ = p;
     f_int ilo_ = ilo;
     f_int ihi_ = ihi;
     f_int lda1_ = lda1;
     f_int lda2_ = lda2;
     f_int ldtau_ = ldtau;
     f_int ret_ = 0;

     mb03vd_(&n_, &p_, &ilo_, &ihi_, a, &lda1_, &lda2_, tau, &ldtau_, dwork, &ret_);

     if (ret_<0) {
       casadi_error("mb03vd wrong arguments:" << ret_);
     } else if (ret_>0) {
       casadi_error("mb03vd error code:" << ret_);
     }


  }

  void slicot_mb03vy(int n, int p, int ilo, int ihi, double * a, int lda1, int lda2,
                     const double * tau, int ldtau, double * dwork, int ldwork) {
     if (dwork==0) {
       std::vector<double> work = std::vector<double>(4*n);
       slicot_mb03vy(n, p, ilo, ihi, a, lda1, lda2, tau, ldtau, &work[0], 4*n);
       return;
     }
     f_int n_ = n;
     f_int p_ = p;
     f_int ilo_ = ilo;
     f_int ihi_ = ihi;
     f_int lda1_ = lda1;
     f_int lda2_ = lda2;
     f_int ldtau_ = ldtau;
     f_int ldwork_ = ldwork;
     f_int ret_=0;
     mb03vy_(&n_, &p_, &ilo_, &ihi_, a, &lda1_, &lda2_, tau, &ldtau_, dwork, &ldwork_, &ret_);

     if (ret_ < 0) {
       casadi_error("mb03vy wrong arguments:" << ret_);
     } else if (ret_>0) {
       casadi_error("mb03vy error code:" << ret_);
     }


  }

  void slicot_mb03wd(char job, char compz, int n, int p, int ilo, int ihi, int iloz, int ihiz,
                     double *h, int ldh1, int ldh2, double* z, int ldz1, int ldz2, double* wr,
                     double *wi, double * dwork, int ldwork) {
      if (dwork==0) {
       std::vector<double> work = std::vector<double>(ihi-ilo+p-1);
       slicot_mb03wd(job, compz, n, p, ilo, ihi, iloz, ihiz, h, ldh1, ldh2, z, ldz1, ldz2, wr, wi,
                     &work[0], ihi-ilo+p-1);
       return;
     }
     f_int n_ = n;
     f_int p_ = p;
     f_int ilo_ = ilo;
     f_int ihi_ = ihi;
     f_int iloz_ = ilo;
     f_int ihiz_ = ihi;
     f_int ldh1_ = ldh1;
     f_int ldh2_ = ldh2;
     f_int ldz1_ = ldz1;
     f_int ldz2_ = ldz2;
     f_int ldwork_ = ldwork;
     f_int ret_ = 0;
     mb03wd_(&job, &compz, &n_, &p_, &ilo_, &ihi_, &iloz_, &ihiz_, h, &ldh1_, &ldh2_,
             z, &ldz1_, &ldz2_, wr, wi, dwork, &ldwork_, &ret_);

     if (ret_<0) {
       casadi_error("mb03wd wrong arguments:" << ret_);
     } else if (ret_>0) {
       casadi_error("mb03wd error code:" << ret_);
     }
  }

  void slicot_periodic_schur(int n, int K, const std::vector< double > & a,
                             std::vector< double > & t,  std::vector< double > & z,
                             std::vector<double> &eig_real, std::vector<double> &eig_imag,
                             double num_zero) {
    t.resize(n*n*K);
    std::vector<double> dwork(std::max(n+K-2, 4*n)+(n-1)*K);
    if (eig_real.size()!=n) {
      eig_real.resize(n);
    }

    if (eig_imag.size()!=n) {
      eig_imag.resize(n);
    }
    slicot_periodic_schur(n, K, get_ptr(a), get_ptr(t), get_ptr(z), get_ptr(dwork), get_ptr(eig_real), get_ptr(eig_imag), num_zero);
  }


  void slicot_periodic_schur(int n, int K, const double* a,
                             double* t,  double * z,
                             double* dwork, double* eig_real,
                             double *eig_imag, double num_zero) {
    int mem_base = std::max(n+K-2, 4*n);
    int mem_needed = mem_base+(n-1)*K;

    // a is immutable, we need a mutable pointer, so we use available buffer
    std::copy(a, a+n*n*K, z);

    slicot_mb03vd(n, K, 1, n, z, n, n, dwork+mem_base, n-1, dwork);
    std::copy(z, z+n*n*K, t);

    slicot_mb03vy(n, K, 1, n, z, n, n, dwork+mem_base, n-1, dwork, mem_needed);

    // Set numerical zeros to zero
    if (num_zero>0) {
      for (int k = 0;k<n*n*K;++k) {
        double &r = t[k];
        if (fabs(r)<num_zero) r = 0.0;
      }
    }

    slicot_mb03wd('S', 'V', n, K, 1, n, 1, n, t, n, n, z, n, n,
                  eig_real, eig_imag, dwork, mem_needed);
  }

  void slicot_periodic_schur(const std::vector< Matrix<double> > & a,
                             std::vector< Matrix<double> > & t, std::vector< Matrix<double> > & z,
                             std::vector< double > & eig_real, std::vector< double > & eig_imag,
                             double num_zero) {
    int K = a.size();
    int n = a[0].size1();
    for (int k=0;k<K;++k) {
      casadi_assert_message(a[k].is_square(), "a must be square");
      casadi_assert_message(a[k].size1()==n, "a must be n-by-n");
      casadi_assert_message(a[k].is_dense(), "a must be dense");
    }


    std::vector<double> a_data(n*n*K);
    // Copy data into consecutive structure
    for (int k=0;k<K;++k) {
      std::copy(a[k].nonzeros().begin(), a[k].nonzeros().end(), a_data.begin()+k*n*n);
    }

    std::vector<double> t_data(n*n*K);
    std::vector<double> z_data(n*n*K);

    slicot_periodic_schur(n, K, a_data, t_data, z_data, eig_real, eig_imag, num_zero);

    t.resize(K);
    z.resize(K);
    for (int k=0;k<K;++k) {
      t[k] = DM::zeros(n, n);
      std::copy(t_data.begin()+k*n*n, t_data.begin()+(k+1)*n*n, t[k].nonzeros().begin());
      z[k] = DM::zeros(n, n);
      std::copy(z_data.begin()+k*n*n, z_data.begin()+(k+1)*n*n, z[k].nonzeros().begin());
    }
  }

} // namespace casadi
