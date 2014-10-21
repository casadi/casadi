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


#include "psd_indef_dple_internal.hpp"
#include "slicot_tools.hpp"
#include <cassert>
#include "../../core/std_vector_tools.hpp"
#include "../../core/matrix/matrix_tools.hpp"
#include "../../core/mx/mx_tools.hpp"
#include "../../core/sx/sx_tools.hpp"
#include "../../core/function/mx_function.hpp"
#include "../../core/function/sx_function.hpp"

#include "../../core/profiling.hpp"
#include "../../core/casadi_options.hpp"
#include <ctime>

#include <numeric>

INPUTSCHEME(DPLEInput)
OUTPUTSCHEME(DPLEOutput)

using namespace std;
namespace casadi {

  extern "C"
  int CASADI_DPLESOLVER_SLICOT_EXPORT
  casadi_register_dplesolver_slicot(DpleInternal::Plugin* plugin) {
    plugin->creator = PsdIndefDpleInternal::creator;
    plugin->name = "slicot";
    plugin->doc = PsdIndefDpleInternal::meta_doc.c_str();;
    plugin->version = 21;
    return 0;
  }

  extern "C"
  void CASADI_DPLESOLVER_SLICOT_EXPORT casadi_load_dplesolver_slicot() {
    DpleInternal::registerPlugin(casadi_register_dplesolver_slicot);
  }

  PsdIndefDpleInternal::PsdIndefDpleInternal(const DpleStructure & st,
                                             int nrhs,
                                             bool transp)
      : DpleInternal(st, nrhs, transp) {

    // set default options
    setOption("name", "unnamed_psd_indef_dple_solver"); // name of the function

    setOption("pos_def", false);
    setOption("const_dim", true);

    addOption("linear_solver",            OT_STRING, GenericType(),
              "User-defined linear solver class. Needed for sensitivities.");
    addOption("linear_solver_options",    OT_DICTIONARY,   GenericType(),
              "Options to be passed to the linear solver.");
    addOption("psd_num_zero",             OT_REAL,         1e-12,
              "Numerical zero used in Periodic Schur decomposition with slicot."
              "This option is needed when your systems has Floquet multipliers"
              "zero or close to zero");
  }

  PsdIndefDpleInternal::~PsdIndefDpleInternal() {

  }

  void PsdIndefDpleInternal::init() {

    DpleInternal::init();

    casadi_assert_message(!pos_def_,
                          "pos_def option set to True: Solver only handles the indefinite case.");
    casadi_assert_message(const_dim_,
                          "const_dim option set to False: Solver only handles the True case.");

    DenseIO::init();

    //for (int k=0;k<K_;k++) {
    //  casadi_assert_message(A_[k].isDense(), "Solver requires arguments to be dense.");
    //  casadi_assert_message(V_[k].isDense(), "Solver requires arguments to be dense.");
    //}

    n_ = A_[0].size1();

    // Allocate data structures
    VZ_.resize(n_*n_*K_);
    T_.resize(n_*n_*K_);
    Z_.resize(n_*n_*K_);
    X_.resize(n_*n_*K_);

    Xbar_.resize(n_*n_*K_);

    nnKa_.resize(K_, DMatrix::zeros(n_, n_));
    nnKb_.resize(K_, DMatrix::zeros(n_, n_));

    eig_real_.resize(n_);
    eig_imag_.resize(n_);


    F_.resize(2*2*n_*K_);

    FF_.resize(2*2*K_);

    dwork_.resize(std::max(n_+K_-2, 4*n_)+(n_-1)*K_+2*n_);

    // There can be at most n partitions
    partition_.reserve(n_);

    if (hasSetOption("linear_solver")) {
      std::string linear_solver_name = getOption("linear_solver");

      // Construct linear solvers for low-order Discrete Periodic Sylvester Equations
      // I00X
      // XI00
      // 0XI0
      // 00XI
      //  Special case K=1
      // I+X
      // Solver complexity:  K
      dpse_solvers_.resize(3);
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
        LinearSolver solver(linear_solver_name, sp, 1);
        solver.init();

        dpse_solvers_[i].resize(n_*(n_+1)/2, solver);
        for (int k=0;k<n_*(n_+1)/2;++k) {
          dpse_solvers_[i][k].makeUnique();
          dpse_solvers_[i][k].setInput(1, LINSOL_A);
        }

      }

    } else {
      casadi_error("Must set linear_solver option.");
    }

    if (CasadiOptions::profiling && CasadiOptions::profilingBinary) {
      profileWriteName(CasadiOptions::profilingLog, this, "PsdIndefSolver",
                       ProfilingData_FunctionType_Other, 4);

      profileWriteSourceLine(CasadiOptions::profilingLog, this, 0, "periodic schur form", -1);
      profileWriteSourceLine(CasadiOptions::profilingLog, this, 1, "nominal", -1);
      profileWriteSourceLine(CasadiOptions::profilingLog, this, 2, "forward", -1);
      profileWriteSourceLine(CasadiOptions::profilingLog, this, 3, "adjoint", -1);
    }

    psd_num_zero_ = getOption("psd_num_zero");

  }

  /// \cond INTERNAL
  inline int PsdIndefDpleInternal::partindex(int i, int j, int k, int r, int c) {
    return k*n_*n_+(partition_[i]+r)*n_ + partition_[j]+c;
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


  void PsdIndefDpleInternal::evaluate() {
    // Obtain a periodic Schur form

    DenseIO::readInputs();

    double t_linear_solve_ = 0; // Time spent in linear solve
    double t_psd_ = 0;          // Time spent in Periodic Schur decomposition
    double t_total_ = 0;        // Time spent in total

    double time_total_start = clock();

    // Set up timers for profiling
    double time_zero=0;
    double time_start=0;
    double time_stop=0;
    if (CasadiOptions::profiling && CasadiOptions::profilingBinary) {
      time_zero = getRealTime();
      profileWriteEntry(CasadiOptions::profilingLog, this);
    }

    if (CasadiOptions::profiling) {
      time_start = getRealTime(); // Start timer
    }

    // Transpose operation (after #554)
    for (int k=0;k<K_;++k) {
      for (int i=0;i<n_;++i) {
        for (int j=0;j<n_;++j) {
          X_[k*n_*n_+i*n_+j] = inputD(DPLE_A).data()[k*n_*n_+n_*j+i];
        }
      }
    }

    double time_psd_start = clock();
    slicot_periodic_schur(n_, K_, X_, T_, Z_, dwork_, eig_real_, eig_imag_, psd_num_zero_);
    t_psd_+=(clock()-time_psd_start)/CLOCKS_PER_SEC;

    if (CasadiOptions::profiling && CasadiOptions::profilingBinary) {
      time_stop = getRealTime(); // Stop timer
      profileWriteTime(CasadiOptions::profilingLog, this, 0, time_stop-time_start,
                       time_stop-time_zero);
    }

    if (CasadiOptions::profiling) {
      time_start = getRealTime(); // Start timer
    }

    if (error_unstable_) {
      for (int i=0;i<n_;++i) {
        double modulus = sqrt(eig_real_[i]*eig_real_[i]+eig_imag_[i]*eig_imag_[i]);
        casadi_assert_message(modulus+eps_unstable_ <= 1,
          "PsdIndefDpleInternal: system is unstable."
          "Found an eigenvalue " << eig_real_[i] << " + " <<
          eig_imag_[i] << "j, with modulus " << modulus <<
          " (corresponding eps= " << 1-modulus << ")." <<
          std::endl << "Use options and 'error_unstable'"
          "and 'eps_unstable' to influence this message.");
      }
    }

    // Find a block partition of the T hessenberg form
    partition_.resize(1, 0);
    int i = 0;
    int j = 0;
    while (j<n_) {
      while (i<n_ && T_[i+n_*j]!=0) {
        i+=1;
      }
      j = i;
      partition_.push_back(i);
      i += 1;
    }

    // Main loops to loop over blocks of the block-upper triangular A
    // Outer main loop
    for (int l=0;l<partition_.size()-1;++l) {

      // Inner main loop
      for (int r=0;r<l+1;++r) {

        int n1 = partition_[r+1]-partition_[r];
        int n2 = partition_[l+1]-partition_[l];
        int np = n1*n2;

        casadi_assert(n1-1+n2-1>=0);

        LinearSolver & solver = dpse_solvers_[n1-1+n2-1][((l+1)*l)/2+r];

        // ********** START ***************
        // Populate the appropriate solver with kron(Arr, All)
        std::vector<double> &A = solver.input(LINSOL_A).data();

        // Special case if K==1
        if (K_==1) {
          for (int ll=0;ll<np;++ll) {
            for (int m=0;m<np;++m) {
              A[ll*np+m] = -T_[partindex(r, r, 0, ll/n2, m/n2)]*T_[partindex(l, l, 0, ll%n2, m%n2)];
              if (ll==m) {
                A[ll*np+m]+= 1;
              }
            }
          }
        } else {
          // Other cases
          int k;
          for (k=0;k<K_-1;++k) {
            for (int ll=0;ll<np;++ll) {
              for (int m=0;m<np;++m) {
                A[np*(np+1)*((k+1)%K_)+ll*(np+1)+m] =
                    -T_[partindex(r, r, k, ll/n2, m/n2)]*T_[partindex(l, l, k, ll%n2, m%n2)];
              }
            }
          }

          for (int ll=0;ll<np;++ll) {
            for (int m=0;m<np;++m) {
              A[np*(np+1)*((k+1)%K_)+ll*(np+1)+m+1] =
                -T_[partindex(r, r, k, ll/n2, m/n2)]*T_[partindex(l, l, k, ll%n2, m%n2)];
            }
          }
        }
        // ********** STOP ***************
        // Solve Discrete Periodic Sylvester Equation Solver

        double time_linear_solve_start = clock();
        solver.prepare();
        t_linear_solve_ += (clock()-time_linear_solve_start)/CLOCKS_PER_SEC;

      }
    }

    if (CasadiOptions::profiling && CasadiOptions::profilingBinary) {
      time_stop = getRealTime(); // Stop timer
      profileWriteTime(CasadiOptions::profilingLog, this, 1, time_stop-time_start,
                       time_stop-time_zero);
    }

    if (!transp_) {
      for (int d=0;d<nrhs_;++d) {

        if (CasadiOptions::profiling) {
          time_start = getRealTime(); // Start timer
        }

        // ********** START ***************
        // V = blocks([mul([sZ[k].T, V[k], sZ[k]]) for k in range(p)])

        for (int k=0;k<K_;++k) { // K
          // n^2 K
          nnKa_[k].set(0.0);
          // nnKa[k] <- V[k]*Z[k+1]
          // n^3 K
          dense_mul_nt(n_, n_, n_,
            &inputD(1+d).data()[k*n_*n_], &Z_[((k+1) % K_)*n_*n_], &nnKa_[k].data()[0]);
          nnKb_[k].set(0.0);
          // nnKb[k] <- Z[k+1]'*V[k]*Z[k+1]
          dense_mul_nn(n_, n_, n_,
            &Z_[((k+1) % K_)*n_*n_], &nnKa_[k].data()[0], &nnKb_[k].data()[0]);
        }

        // ********** STOP ****************

        std::fill(X_.begin(), X_.end(), 0);

        // Main loops to loop over blocks of the block-upper triangular A
        // Outer main loop
        for (int l=0;l<partition_.size()-1;++l) { // n
          // F serves an an accumulator for intermediate summation results
          // n^2 K
          std::fill(F_.begin(), F_.end(), 0);

          // ********** START ***************
          //for i in range(l):
          //  F[i] = [sum(mul(X[i][j][k], A[l][j][k].T) for j in range(l)) for k in range(p) ]

          for (int i=0;i<l;++i) { // n^2
            int na1 = partition_[i+1]-partition_[i];
            int nb2 = partition_[l+1]-partition_[l];
            for (int j=0;j<l;++j) { // n^3
              int na2 = partition_[j+1]-partition_[j];
              for (int k=0;k<K_;++k) {
                for (int ii=0;ii<na1;++ii) {
                  for (int jj=0;jj<nb2;++jj) {
                    for (int kk=0;kk<na2;++kk) {
                      // n^3 K
                      F_[k*4*n_+4*i+2*ii+jj] +=
                          X_[partindex(i, j, k, ii, kk)]*T_[partindex(l, j, k, jj, kk)];
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
            int na1 = partition_[r+1]-partition_[r];
            int nb2 = partition_[l+1]-partition_[l];

            if (r==l) {
              for (int j=0;j<l;++j) { // n^3
                int na2 = partition_[j+1]-partition_[j];
                for (int k=0;k<K_;++k) { // n^3 K
                  for (int ii=0;ii<na1;++ii) {
                    for (int jj=0;jj<nb2;++jj) {
                      for (int kk=0;kk<na2;++kk) {
                        F_[k*4*n_+4*r+2*ii+jj] +=
                            X_[partindex(r, j, k, ii, kk)]*T_[partindex(l, j, k, jj, kk)];
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
            std::fill(FF_.begin(), FF_.end(), 0);
            {
              int na1 = partition_[r+1]-partition_[r];
              for (int i=0;i<r;++i) { // n^3
                int nb2 = partition_[l+1]-partition_[l];
                int na2 = partition_[i+1]-partition_[i];
                for (int k=0;k<K_;++k) { // n^3 K
                  for (int ii=0;ii<na1;++ii) {
                    for (int jj=0;jj<nb2;++jj) {
                      for (int kk=0;kk<na2;++kk) {
                        FF_[k*4+2*ii+jj] +=
                            T_[partindex(r, i, k, ii, kk)]*X_[partindex(i, l, k, kk, jj)];
                      }
                    }
                  }
                }
              }
            }
            // ********** STOP ***************

            int n1 = partition_[r+1]-partition_[r];
            int n2 = partition_[l+1]-partition_[l];
            int np = n1*n2;

            LinearSolver & solver = dpse_solvers_[n1-1+n2-1][((l+1)*l)/2+r];

            // M <- V
            std::vector<double> &M = solver.input(LINSOL_B).data();
            for (int k=0;k<K_;++k) {
              for (int ll=0;ll<n1;++ll) {
                for (int m=0;m<n2;++m) {
                  // n^2 K
                  M[np*((k+1)%K_)+ll*n2+m] =
                    nnKb_[k].data()[(partition_[r]+ll)*n_ + partition_[l]+m];
                }
              }
            }

            // ********** START ***************
            // M+= [sum(mul(A[r][i][k], F[i][k])  for i in range(r+1)) for k in rang(p)]
            {
              int na1 = partition_[r+1]-partition_[r];
              int nb2 = partition_[l+1]-partition_[l];

              for (int i=0;i<r+1;++i) { // n^3
                int na2 = partition_[i+1]-partition_[i];
                for (int k=0;k<K_;++k) { // n^3 K
                  for (int ii=0;ii<na1;++ii) {
                    for (int jj=0;jj<nb2;++jj) {
                      for (int kk=0;kk<na2;++kk) {
                        M[np*((k+1)%K_)+ii*n2+jj] +=
                            T_[partindex(r, i, k, ii, kk)]*F_[k*4*n_+4*i+2*kk+jj];
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
              int na1 = partition_[r+1]-partition_[r];
              int na2 = partition_[l+1]-partition_[l];
              int nb2 = partition_[l+1]-partition_[l];
              for (int k=0;k<K_;++k) { // n^2 K
                for (int ii=0;ii<na1;++ii) {
                  for (int jj=0;jj<nb2;++jj) {
                    for (int kk=0;kk<na2;++kk) {
                      M[np*((k+1)%K_)+ii*n2+jj] += FF_[k*4+2*ii+kk]*T_[partindex(l, l, k, jj, kk)];
                    }
                  }
                }
              }
            }
            // ********** STOP ***************

            // Critical observation: Prepare step is not needed
            double time_linear_solve_start = clock();
            // n^2 K
            solver.solve(true);
            t_linear_solve_ += (clock()-time_linear_solve_start)/CLOCKS_PER_SEC;

            // Extract solution and store it in X
            std::vector<double> & sol = solver.output().data();

            // ********** START ***************
            for (int ii=0;ii<partition_[r+1]-partition_[r];++ii) {
              for (int jj=0;jj<partition_[l+1]-partition_[l];++jj) {
                for (int k=0;k<K_;++k) { // n^2 K
                  X_[partindex(r, l, k, ii, jj)] = sol[n1*n2*k+n2*ii+jj];
                }
              }
            }

            for (int ii=0;ii<partition_[r+1]-partition_[r];++ii) {
              for (int jj=0;jj<partition_[l+1]-partition_[l];++jj) {
                for (int k=0;k<K_;++k) { // n^2 K
                  X_[partindex(l, r, k, jj, ii)] = sol[n1*n2*k+n2*ii+jj];
                }
              }
            }
            // ********** STOP ***************

          }

          // n^3 K
          outputD(d).set(0.0);
        }

        for (int k=0;k<K_;++k) {
          nnKa_[k].set(0.0);

          // nnKa[k] <- V[k]*Z[k]'
          // n^3 K
          dense_mul_nn(n_, n_, n_, &X_[k*n_*n_], &Z_[k*n_*n_], &nnKa_[k].data()[0]);
          // output <- Z[k]*V[k]*Z[k]'
          dense_mul_tn(n_, n_, n_, &Z_[k*n_*n_], &nnKa_[k].data()[0],
                       &outputD(d).data()[k*n_*n_]);
        }

        if (CasadiOptions::profiling && CasadiOptions::profilingBinary) {
          time_stop = getRealTime(); // Stop timer
          profileWriteTime(CasadiOptions::profilingLog, this, 2, time_stop-time_start,
                           time_stop-time_zero);
        }

      }
    } else { // Transposed

      for (int d=0;d<nrhs_;++d) {

        if (CasadiOptions::profiling) {
          time_start = getRealTime(); // Start timer
        }

        DMatrix &P_bar = inputD(1+d);
        std::vector<double> &Vbar = outputD(d).data();
        std::fill(Vbar.begin(), Vbar.end(), 0);
        std::fill(Xbar_.begin(), Xbar_.end(), 0);

        // X_bar = [mul([Z[k].T, X_bar[k] , Z[k]]) for k in range(p)]
        for (int k=0;k<K_;++k) {
          nnKa_[k].set(0.0);

          // n^3 K
          // nnKa[k] <- nnKb*Z[k]
          dense_mul_nt(n_, n_, n_, &P_bar.data()[n_*n_*k], &Z_[k*n_*n_], &nnKa_[k].data()[0]);
          // Xbar <- Z[k]*V[k]*Z[k]'
          dense_mul_nn(n_, n_, n_, &Z_[k*n_*n_], &nnKa_[k].data()[0], &Xbar_[k*n_*n_]);
        }




        // Main loops to loop over blocks of the block-upper triangular A
        // Outer main loop
        for (int l=partition_.size()-2;l>=0;--l) { // n

          std::fill(F_.begin(), F_.end(), 0);

          // Inner main loop
          for (int r=l;r>=0;--r) { // n^2

            int n1 = partition_[r+1]-partition_[r];
            int n2 = partition_[l+1]-partition_[l];
            int np = n1*n2;

            LinearSolver & solver = dpse_solvers_[n1-1+n2-1][((l+1)*l)/2+r];

            std::vector<double> & B = solver.input(LINSOL_B).data();


            // ********** START ***************
            for (int ii=0;ii<partition_[r+1]-partition_[r];++ii) {
              for (int jj=0;jj<partition_[l+1]-partition_[l];++jj) {
                for (int k=0;k<K_;++k) {
                  // n^2 K
                  B[n1*n2*k+n2*ii+jj] = Xbar_[partindex(r, l, k, ii, jj)];
                }
              }
            }
            if (r!=l) {
              for (int ii=0;ii<partition_[r+1]-partition_[r];++ii) {
                for (int jj=0;jj<partition_[l+1]-partition_[l];++jj) {
                  for (int k=0;k<K_;++k) {
                    B[n1*n2*k+n2*ii+jj]+=Xbar_[partindex(l, r, k, jj, ii)];
                  }
                }
              }
            }
            // ********** STOP ***************

            double time_linear_solve_start = clock();
            // n^2 K
            solver.solve(false);
            t_linear_solve_ += (clock()-time_linear_solve_start)/CLOCKS_PER_SEC;

            // for k in range(p): V_bar[r][l][k]+=M_bar[k]
            std::vector<double> &Mbar = solver.output().data();

            for (int k=0;k<K_;++k) {
              for (int ll=0;ll<n1;++ll) {
                for (int m=0;m<n2;++m) {
                   // n^2 K
                   Vbar[partindex(r, l, k, ll, m)]+= Mbar[np*((k+1)%K_)+ll*n2+m];
                }
              }
            }

            // ********** START ***************
            //     for k in range(p):
            //       FF_bar[k] += mul(M_bar[k], A[l][l][k])
            // n^2 K
            std::fill(FF_.begin(), FF_.end(), 0);
            {
              int na1 = partition_[r+1]-partition_[r];
              int na2 = partition_[l+1]-partition_[l];
              int nb2 = partition_[l+1]-partition_[l];
              for (int k=0;k<K_;++k) {
                for (int ii=0;ii<na1;++ii) {
                  for (int jj=0;jj<nb2;++jj) {
                    for (int kk=0;kk<na2;++kk) {
                      // n^2 K
                      FF_[k*4+2*ii+kk] +=
                        Mbar[np*((k+1)%K_)+ii*n2+jj]*T_[partindex(l, l, k, jj, kk)];
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
              int na1 = partition_[r+1]-partition_[r];
              for (int i=0;i<r;++i) { // n^3
                int nb2 = partition_[l+1]-partition_[l];
                int na2 = partition_[i+1]-partition_[i];
                for (int k=0;k<K_;++k) {
                  for (int ii=0;ii<na1;++ii) {
                    for (int jj=0;jj<nb2;++jj) {
                      for (int kk=0;kk<na2;++kk) {
                        // n^3 K
                        Xbar_[partindex(i, l, k, kk, jj)] +=
                            T_[partindex(r, i, k, ii, kk)]*FF_[k*4+2*ii+jj];
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
              int na1 = partition_[r+1]-partition_[r];
              int nb2 = partition_[l+1]-partition_[l];

              for (int i=0;i<r+1;++i) { // n^3
                int na2 = partition_[i+1]-partition_[i];
                for (int k=0;k<K_;++k) {
                  for (int ii=0;ii<na1;++ii) {
                    for (int jj=0;jj<nb2;++jj) {
                      for (int kk=0;kk<na2;++kk) {
                        // n^3 K
                        F_[k*4*n_+4*i+2*kk+jj] +=
                            T_[partindex(r, i, k, ii, kk)]*Mbar[np*((k+1)%K_)+ii*n2+jj];
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
              int na1 = partition_[r+1]-partition_[r];
              int nb2 = partition_[l+1]-partition_[l];

              for (int j=0;j<l;++j) { // n^3
                int na2 = partition_[j+1]-partition_[j];
                for (int k=0;k<K_;++k) {
                  for (int ii=0;ii<na1;++ii) {
                    for (int jj=0;jj<nb2;++jj) {
                      for (int kk=0;kk<na2;++kk) {
                        // n^3 K
                        Xbar_[partindex(r, j, k, ii, kk)] +=
                            F_[k*4*n_+4*r+2*ii+jj]*T_[partindex(l, j, k, jj, kk)];
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
            int na1 = partition_[i+1]-partition_[i];
            int nb2 = partition_[l+1]-partition_[l];
            for (int j=0;j<l;++j) { // n^3
              int na2 = partition_[j+1]-partition_[j];
              for (int k=0;k<K_;++k) {
                for (int ii=0;ii<na1;++ii) {
                  for (int jj=0;jj<nb2;++jj) {
                    for (int kk=0;kk<na2;++kk) {
                      // n^3 K
                      Xbar_[partindex(i, j, k, ii, kk)] +=
                          F_[k*4*n_+4*i+2*ii+jj]*T_[partindex(l, j, k, jj, kk)];
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
          nnKa_[k].set(0.0);
          // n^3 K
          dense_mul_nn(n_, n_, n_, &Vbar[n_*n_*k], &Z_[((k+1)%K_)*n_*n_], &nnKa_[k].data()[0]);
          std::fill(&Vbar[n_*n_*k], &Vbar[n_*n_*(k+1)], 0);
          dense_mul_tn(n_, n_, n_, &Z_[((k+1)%K_)*n_*n_], &nnKa_[k].data()[0], &Vbar[n_*n_*k]);
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

       if (CasadiOptions::profiling && CasadiOptions::profilingBinary) {
          time_stop = getRealTime(); // Stop timer
          profileWriteTime(CasadiOptions::profilingLog, this, 3, time_stop-time_start,
                           time_stop-time_zero);
       }

      }
    }

    if (CasadiOptions::profiling && CasadiOptions::profilingBinary) {
      time_stop = getRealTime();
      profileWriteExit(CasadiOptions::profilingLog, this, time_stop-time_zero);
    }

    t_total_ += (clock()-time_total_start)/CLOCKS_PER_SEC;

    if (gather_stats_) {
      stats_["t_psd"] = t_psd_;
      stats_["t_total"] = t_total_;
      stats_["t_linear_solve"] = t_linear_solve_;
    }

    DenseIO::writeOutputs();

  }

  Function PsdIndefDpleInternal::getDerivative(int nfwd, int nadj) {

    // Base:
    // P_0 P_1 P_2 .. P_{nrhs-1} = f( A Q_0 Q_1 Q_2 .. Q_{nrhs-1})

    /* Allocate output list for derivative
    *
    * Three parts:
    *
    * 1)  [P_0 .. P_{nrhs-1}]
    * 2)  [P_0^f0 .. P_{nrhs-1}^f0] ...
    *       [P_0^f{nfwd-1} .. P_{nrhs-1}^f{nfwd-1}]
    * 3)  [A^b0 Q_0^b0 .. Q_{nrhs-1}^b0] ...
    *       [A^b{nadj-1} Q_0^b{nadj-1} .. Q_{nrhs-1}^b{nadj-1}]
    */
    std::vector<MX> outs_new((nrhs_+1)*nadj+nrhs_*(nfwd+1), 0);

    /* Allocate input list for derivative and populate with symbolics
    * Three parts:
    *
    * 1)  [A Q_0 .. Q_{nrhs-1}]
    * 2)  [A^f0 Q_0^f0 .. Q_{nrhs-1}^f0] ...
    *       [A^f{nfwd-1} Q_0^f{nfwd-1} .. Q_{nrhs-1}^f{nfwd-1}]
    * 3)  [P_0^b0 .. P_{nrhs-1}^b0] ...
    *       [P_0^b{nadj-1} .. P_{nrhs-1}^b{nadj-1}]
    */
    std::vector<MX> ins_new((nrhs_+1)*(nfwd+1)+nrhs_*nadj);

    // Part 1
    ins_new[0] = MX::sym("A", input(DPLE_A).sparsity());
    for (int i=0; i<nrhs_; ++i) {
      ins_new[i+1] = MX::sym("Q", input(DPLE_V).sparsity());
    }

    // Part 2
    for (int q=0; q<nrhs_; ++q) {
      for (int k=0;k<nfwd;++k) {
        MX& Qf  = ins_new.at((nrhs_+1)*(k+1)+q+1);
        MX& Af  = ins_new.at((nrhs_+1)*(k+1));

        Qf = MX::sym("Qf", input(DPLE_V).sparsity());
        Af = MX::sym("Af", input(DPLE_A).sparsity());
      }
    }

    // Part 3
    for (int q=0; q<nrhs_; ++q) {
      for (int k=0;k<nadj;++k) {
          MX& Pb  = ins_new.at((nrhs_+1)*(nfwd+1)+nrhs_*k+q);

          Pb = MX::sym("Pb", output(DPLE_P).sparsity());
      }
    }
    // Create a call to yourself
    //
    // P_0 .. P_{nrhs-1} = f( A Q_0 .. Q_{nrhs-1})
    std::vector<MX> ins_P;
    ins_P.insert(ins_P.begin(), ins_new.begin(), ins_new.begin()+(nrhs_+1));
    std::vector<MX> Ps = shared_from_this<Function>().call(ins_P);

    for (int i=0;i<nrhs_;++i) {
      outs_new[i] = Ps[i];
    }

    // Prepare a solver for forward seeds
    PsdIndefDpleInternal* node = new PsdIndefDpleInternal(st_, nfwd, transp_);
    node->setOption(dictionary());

    DpleSolver f;
    f.assignNode(node);
    f.init();

    // Prepare a solver for adjoint seeds
    PsdIndefDpleInternal* node2 = new PsdIndefDpleInternal(st_, nadj, !transp_);
    node2->setOption(dictionary());

    DpleSolver b;
    b.assignNode(node2);
    b.init();

    for (int q=0; q<nrhs_; ++q) {

      // Forward
      /* P_q^f0 .. P_q^f{nfwd-1} = f(A,
      *          Q_q^f0 + A P_q (A^f0)^T + A^f0 P_q A^T,
      *          ...
      *          Q_q^f{nfwd-1} + A P_q (A^f{nfwd-1})^T + A^f{nfwd-1} P_q A^T)
      */
      std::vector<MX> ins_f;
      const MX& A  = ins_new[0];
      std::vector<MX> As = horzsplit(A, n_);
      const MX& P  = Ps[q];
      std::vector<MX> Ps_ = horzsplit(P, n_);

      ins_f.push_back(A);
      for (int k=0;k<nfwd;++k) {
        const MX& Qf  = ins_new.at((nrhs_+1)*(k+1)+q+1);
        const MX& Af  = ins_new.at((nrhs_+1)*(k+1));

        std::vector<MX> Qfs = horzsplit(Qf, n_);
        std::vector<MX> Afs = horzsplit(Af, n_);

        std::vector<MX> sum(K_, 0);
        for (int i=0;i<K_;++i) {
          // Note: Qf is symmetrised here
          MX temp;
          if (transp_) {
            temp = mul(As[i].T(), mul(Ps_[i], Afs[i])) + Qfs[i]/2;
          } else {
            temp = mul(As[i], mul(Ps_[i], Afs[i].T())) + Qfs[i]/2;
          }
          sum[i] = temp + temp.T();
        }
        ins_f.push_back(horzcat(sum));
      }

      std::vector<MX> outs = f.call(ins_f);
      for (int i=0;i<nfwd;++i) {
        outs_new.at(nrhs_+q+nrhs_*i) = outs[i];
      }

      // Adjoint
      /* rev(Q_b^b0) .. rev(Q_b^b{nadj-1}) += f(rev(A),
      *          rev(P_q^b0) ... rev(P_q^b{nadj-1})
      *         )
      *
      *  A^b0 += 2 Q_q^b0 A P_q
      *   ....
      *  A^b{nadj-1} += 2 Q_q^b{nadj-1} A P_q
      */
      std::vector<MX> ins_b;
      ins_b.push_back(A);
      for (int k=0;k<nadj;++k) {
        const MX& Pb  = ins_new.at((nrhs_+1)*(nfwd+1)+nrhs_*k+q);
        // Symmetrise P_q^bk
        std::vector<MX> Pbs = horzsplit(Pb, n_);
        for (int i=0;i<K_;++i) {
          Pbs[i]+= Pbs[i].T();
        }

        ins_b.push_back(horzcat(Pbs)/2);
      }

      outs = b.call(ins_b);
      for (int i=0;i<nadj;++i) {
        MX& Qb = outs_new.at(nrhs_*(nfwd+1)+(nrhs_+1)*i+q+1);

        Qb += outs[i];
        std::vector<MX> Qbs = horzsplit(Qb, n_);
        MX& Ab = outs_new.at(nrhs_*(nfwd+1)+(nrhs_+1)*i);

        std::vector<MX> sum(K_, 0);
        for (int j=0;j<K_;++j) {
          if (transp_) {
            sum[j]+= (2*mul(Ps_[j], mul(As[j], Qbs[j])));
          } else {
            sum[j]+= 2*mul(Qbs[j], mul(As[j], Ps_[j]));
          }
        }
        Ab += horzcat(sum);
      }

    }

    MXFunction ret(ins_new, outs_new);
    ret.init();

    return ret;

  }

  void PsdIndefDpleInternal::deepCopyMembers(
      std::map<SharedObjectNode*, SharedObject>& already_copied) {
    DpleInternal::deepCopyMembers(already_copied);
  }

  PsdIndefDpleInternal* PsdIndefDpleInternal::clone() const {
    // Return a deep copy
    PsdIndefDpleInternal* node = new PsdIndefDpleInternal(st_, nrhs_, transp_);
    node->setOption(dictionary());
    return node;
  }

} // namespace casadi
