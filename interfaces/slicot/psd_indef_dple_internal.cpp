/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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
#include "../../symbolic/stl_vector_tools.hpp"
#include "../../symbolic/matrix/matrix_tools.hpp"
#include "../../symbolic/mx/mx_tools.hpp"
#include "../../symbolic/sx/sx_tools.hpp"
#include "../../symbolic/fx/mx_function.hpp"
#include "../../symbolic/fx/sx_function.hpp"

#include <numeric>

INPUTSCHEME(DPLEInput)
OUTPUTSCHEME(DPLEOutput)

using namespace std;
namespace CasADi{

  PsdIndefDpleInternal::PsdIndefDpleInternal(const std::vector< CRSSparsity > & A, const std::vector< CRSSparsity > &V, int nfwd, int nadj) : DpleInternal(A,V, nfwd, nadj) {
  
    // set default options
    setOption("name","unnamed_psd_indef_dple_solver"); // name of the function 

    setOption("pos_def",false);
    setOption("const_dim",true);
    
    addOption("linear_solver",            OT_LINEARSOLVER, GenericType(), "User-defined linear solver class. Needed for sensitivities.");
    addOption("linear_solver_options",    OT_DICTIONARY,   GenericType(), "Options to be passed to the linear solver.");
    
  }

  PsdIndefDpleInternal::~PsdIndefDpleInternal(){ 

  }

  void PsdIndefDpleInternal::init(){
  
    DpleInternal::init();

    casadi_assert_message(!pos_def_,"pos_def option set to True: Solver only handles the indefinite case.");
    casadi_assert_message(const_dim_,"const_dim option set to False: Solver only handles the True case.");
    
    n_ = A_[0].size1();
    
    // Allocate data structures
    VZ_.resize(n_*n_*K_);
    T_.resize(n_*n_*K_);
    Z_.resize(n_*n_*K_);
    X_.resize(n_*n_*K_);
    dX_.resize(n_*n_*K_);
    
    Xbar_.resize(n_*n_*K_);
    
    nnKa_.resize(K_,DMatrix::zeros(n_,n_));
    nnKb_.resize(K_,DMatrix::zeros(n_,n_));
    
    eig_real_.resize(n_);
    eig_imag_.resize(n_);
    
    
    F_.resize(2*2*n_*K_);
    
    FF_.resize(2*2*K_);
    
    dwork_.resize(std::max(n_+K_-2,4*n_)+(n_-1)*K_+2*n_);
    
    // There can be at most n partitions
    partition_.reserve(n_);

    if(hasSetOption("linear_solver")){
      linearSolverCreator linear_solver_creator = getOption("linear_solver");
      
      // Construct linear solvers for low-order Discrete Periodic Sylvester Equations
      // I00X
      // XI00 
      // 0XI0
      // 00XI
      //  Special case K=1
      // I+X
      dpse_solvers_.resize(3);
      for (int i=0;i<3;++i) {
        int np = std::pow(2,i);
          
        CRSSparsity sp;
        if (K_==1) {
          sp = sp_dense(np,np);
        } else {
          std::vector<int> row_ind = range(0,np*(np+1)*K_+np+1,np+1);
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
 
          sp = CRSSparsity(np*K_,np*K_,col,row_ind);
        }
        LinearSolver solver = linear_solver_creator(sp,1);
        solver.init();
        
        dpse_solvers_[i].resize(n_*(n_+1)/2,solver);
        for (int k=0;k<n_*(n_+1)/2;++k) {
          dpse_solvers_[i][k].makeUnique();
          dpse_solvers_[i][k].setInput(1,LINSOL_A);
        }
        
      }
        
    } else {
      casadi_error("Must set linear_solver option.");
    }

  }
  
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
  
  
  void PsdIndefDpleInternal::evaluate(){
    // Obtain a periodic Schur form
    slicot_periodic_schur(n_,K_,input(DPLE_A).data(),T_,Z_,dwork_,eig_real_,eig_imag_);
    
    if (error_unstable_) {
      for (int i=0;i<n_;++i) {
        double modulus = sqrt(eig_real_[i]*eig_real_[i]+eig_imag_[i]*eig_imag_[i]);
        casadi_assert_message(modulus+eps_unstable_ <= 1,"PsdIndefDpleInternal: system is unstable. Found an eigenvalue " << eig_real_[i] << " + " << eig_imag_[i] << "j, with modulus " << modulus << " (corresponding eps= " << 1-modulus << ")." << std::endl << "Use options and 'error_unstable' and 'eps_unstable' to influence this message.");
      }
    }
    
    // Find a block partition of the T hessenberg form
    partition_.resize(1,0);
    int i = 0;
    int j = 0;
    while(j<n_) {
      while(i<n_ && T_[i+n_*j]!=0) {
        i+=1;
      }
      j = i;
      partition_.push_back(i);
    }
    
    // ********** START ***************
    // V = blocks([mul([sZ[k].T,Vs[k],sZ[k]]) for k in range(p)])
    
    for(int k=0;k<K_;++k) {
      nnKa_[k].set(0.0);
      nnKb_[k].set(0.0);
      // nnKa[k] <- V[k]*Z[k+1]
      dense_mul_nt(n_,n_,n_,&input(DPLE_V).data()[k*n_*n_],&Z_[((k+1) % K_)*n_*n_],&nnKa_[k].data()[0]);
      // nnKb[k] <- Z[k+1]'*V[k]*Z[k+1]
      dense_mul_nn(n_,n_,n_,&Z_[((k+1) % K_)*n_*n_],&nnKa_[k].data()[0],&nnKb_[k].data()[0]);
    }
    
    // ********** STOP ****************
    
    // Main loops to loop over blocks of the block-upper triangular A
    // Outer main loop
    for (int l=0;l<partition_.size()-1;++l) {
      // F serves an an accumulator for intermediate summation results
      std::fill(F_.begin(),F_.end(),0);
      
      // ********** START ***************
      //for i in range(l):
      //  F[i] = [sum(mul(X[i][j][k],A[l][j][k].T) for j in range(l)) for k in range(p) ]

      for (int i=0;i<l;++i) {
        int na1 = partition_[i+1]-partition_[i];
        int nb2 = partition_[l+1]-partition_[l];
        for (int j=0;j<l;++j) {
          int na2 = partition_[j+1]-partition_[j];
          for (int k=0;k<K_;++k) {
            for (int ii=0;ii<na1;++ii) {
              for (int jj=0;jj<nb2;++jj) {
                for (int kk=0;kk<na2;++kk) {
                  F_[k*4*n_+4*i+2*ii+jj] += X_[partindex(i,j,k,ii,kk)]*T_[partindex(l,j,k,jj,kk)];
                }
              }
            }
          }
        }
      }
      // ********** STOP ***************

      // Inner main loop
      for (int r=0;r<l+1;++r) {
      
        // ********** START ***************
        // F[r] = [sum(mul(X[r][j][k],A[l][j][k].T) for j in range(l)) for k in range(p) ]
        int na1 = partition_[r+1]-partition_[r];
        int nb2 = partition_[l+1]-partition_[l];

        for (int k=0;k<K_;++k) {
          std::fill(F_.begin()+k*4*n_+4*r,F_.begin()+k*4*n_+4*r+4,0);
        }
          
        for (int j=0;j<l;++j) {
          int na2 = partition_[j+1]-partition_[j];
          for (int k=0;k<K_;++k) {
            for (int ii=0;ii<na1;++ii) {
              for (int jj=0;jj<nb2;++jj) {
                for (int kk=0;kk<na2;++kk) {
                  F_[k*4*n_+4*r+2*ii+jj] += X_[partindex(r,j,k,ii,kk)]*T_[partindex(l,j,k,jj,kk)];
                }
              }
            }
          }
        }
        // ********** STOP ***************
        
        // ********** START ***************
        // FF =   [sum(mul(A[r][i][k],X[i][l][k]) for i in range(r)) for k in range(p)]
        // Each entry of FF is na1-by-na2
        std::fill(FF_.begin(),FF_.end(),0);
        { 
          int na1 = partition_[r+1]-partition_[r];
          for (int i=0;i<r;++i) {
            int nb2 = partition_[l+1]-partition_[l];
            int na2 = partition_[i+1]-partition_[i];
            for (int k=0;k<K_;++k) {
              for (int ii=0;ii<na1;++ii) {
                for (int jj=0;jj<nb2;++jj) {
                  for (int kk=0;kk<na2;++kk) {
                    FF_[k*4+2*ii+jj] += T_[partindex(r,i,k,ii,kk)]*X_[partindex(i,l,k,kk,jj)];
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
              M[np*((k+1)%K_)+ll*n2+m] = nnKb_[k].data()[(partition_[r]+ll)*n_ + partition_[l]+m];
            }
          }
        }

        // ********** START ***************
        // M+= [sum(mul(A[r][i][k],F[i][k])  for i in range(r+1)) for k in rang(p)]
        {
          int na1 = partition_[r+1]-partition_[r];
          int nb2 = partition_[l+1]-partition_[l];
          
          for (int i=0;i<r+1;++i) {
            int na2 = partition_[i+1]-partition_[i];
            for (int k=0;k<K_;++k) {
              for (int ii=0;ii<na1;++ii) {
                for (int jj=0;jj<nb2;++jj) {
                  for (int kk=0;kk<na2;++kk) {
                    M[np*((k+1)%K_)+ii*n2+jj] += T_[partindex(r,i,k,ii,kk)]*F_[k*4*n_+4*i+2*kk+jj];
                  }
                }
              }
            }
          }
        }
        // ********** STOP ***************
        
        // ********** START ***************
        // M+= [mul(FF[k],A[l][l][k].T) for k in rang(p)]
        {
          int na1 = partition_[r+1]-partition_[r];
          int na2 = partition_[l+1]-partition_[l];
          int nb2 = partition_[l+1]-partition_[l];
          for (int k=0;k<K_;++k) {
            for (int ii=0;ii<na1;++ii) {
              for (int jj=0;jj<nb2;++jj) {
                for (int kk=0;kk<na2;++kk) {
                  M[np*((k+1)%K_)+ii*n2+jj] += FF_[k*4+2*ii+kk]*T_[partindex(l,l,k,jj,kk)];
                }
              }
            }
          }
        }
        // ********** STOP ***************

        // ********** START ***************
        // Populate the appropriate solver with kron(Arr,All)
        std::vector<double> &A = solver.input(LINSOL_A).data();
        
        // Special case if K==1
        if (K_==1) {
          for (int ll=0;ll<np;++ll) {
            for (int m=0;m<np;++m) {
              A[ll*np+m] = -T_[partindex(r,r,0,ll/n2,m/n2)]*T_[partindex(l,l,0,ll%n2,m%n2)];
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
                A[np*(np+1)*((k+1)%K_)+ll*(np+1)+m] = -T_[partindex(r,r,k,ll/n2,m/n2)]*T_[partindex(l,l,k,ll%n2,m%n2)];
              }
            }
          }

          for (int ll=0;ll<np;++ll) {
            for (int m=0;m<np;++m) {
              A[np*(np+1)*((k+1)%K_)+ll*(np+1)+m+1] = -T_[partindex(r,r,k,ll/n2,m/n2)]*T_[partindex(l,l,k,ll%n2,m%n2)];
            }
          }
        }
        // ********** STOP ***************
        
        // Solve Discrete Periodic Sylvester Equation Solver
        solver.prepare();
        solver.solve(true);

        // Extract solution and store it in X
        std::vector<double> & sol = solver.output().data();

        // ********** START ***************
        for (int ii=0;ii<partition_[r+1]-partition_[r];++ii) {
          for (int jj=0;jj<partition_[l+1]-partition_[l];++jj) {
            for (int k=0;k<K_;++k) {
              X_[partindex(r,l,k,ii,jj)] = sol[n1*n2*k+n2*ii+jj];
            }
          }
        }
        
        if (r!=l) {
          for (int ii=0;ii<partition_[r+1]-partition_[r];++ii) {
            for (int jj=0;jj<partition_[l+1]-partition_[l];++jj) {
              for (int k=0;k<K_;++k) {
                X_[partindex(l,r,k,jj,ii)] = sol[n1*n2*k+n2*ii+jj];
              }
            }
          }
        }
        // ********** STOP ***************
        
      }
    }
    
    output(DPLE_P).set(0.0);

    for(int k=0;k<K_;++k) {
      nnKa_[k].set(0.0);
      
      // nnKa[k] <- V[k]*Z[k]'
      dense_mul_nn(n_,n_,n_,&X_[k*n_*n_],&Z_[k*n_*n_],&nnKa_[k].data()[0]);
      // output <- Z[k]*V[k]*Z[k]'
      dense_mul_tn(n_,n_,n_,&Z_[k*n_*n_],&nnKa_[k].data()[0],&output(DPLE_P).data()[k*n_*n_]);
    }
    
    if (nfwd_==0 && nadj_==0) return;
    // Forward sensitivitites
    

    for (int d=0;d<nfwd_;++d) {
      
      // dV2 = [dV+mul([a_dot,x,a.T])+mul([a,x,a_dot.T]) for vp,a,a_dot,x in zip(Vp,As,Ap,X) ]  
      for(int k=0;k<K_;++k) {
        std::fill(nnKb_[k].begin(),nnKb_[k].end(),0);    
        std::fill(nnKa_[k].begin(),nnKa_[k].end(),0);
        
        // nnKb <- dV
        std::copy(&input(DPLE_NUM_IN*(d+1)+DPLE_V).data()[n_*n_*k],&input(DPLE_NUM_IN*(d+1)+DPLE_V).data()[n_*n_*(k+1)], &nnKb_[k].data()[0]);
        
        // Force seed to be symmetric
        for (int r=0;r<n_;++r) {
          for (int l=0;l<r;++l) {
            double s = (nnKb_[k].data()[r+l*n_]+nnKb_[k].data()[r*n_+l])/2;
            nnKb_[k].data()[r+l*n_] = s;
            nnKb_[k].data()[r*n_+l] = s;
          }
        }
        
        // nnKa <- x*a'
        dense_mul_nt(n_,n_,n_,&output(DPLE_P).data()[n_*n_*k],&input(DPLE_A).data()[n_*n_*k],&nnKa_[k].data()[0]);
        
        // nnKb += da*nnKa
        dense_mul_nn(n_,n_,n_,&input(DPLE_NUM_IN*(d+1)+DPLE_A).data()[n_*n_*k],&nnKa_[k].data()[0],&nnKb_[k].data()[0]);
        
        std::fill(nnKa_[k].begin(),nnKa_[k].end(),0);
        
        // nnKa <- x*da'
        dense_mul_nt(n_,n_,n_,&output(DPLE_P).data()[n_*n_*k],&input(DPLE_NUM_IN*(d+1)+DPLE_A).data()[n_*n_*k],&nnKa_[k].data()[0]);
        
        // nnKb += a*nnKa
        dense_mul_nn(n_,n_,n_,&input(DPLE_A).data()[n_*n_*k],&nnKa_[k].data()[0],&nnKb_[k].data()[0]);
        
      }
      
      // ********** START ***************
      // V = blocks([mul([sZ[k].T,dV2[k],sZ[k]]) for k in range(p)])
      
      for(int k=0;k<K_;++k) {
        nnKa_[k].set(0.0);
        // nnKa[k] <- dV2[k]*Z[k+1]
        dense_mul_nt(n_,n_,n_,&nnKb_[k].data()[0],&Z_[((k+1) % K_)*n_*n_],&nnKa_[k].data()[0]);
        nnKb_[k].set(0.0);
        // nnKb[k] <- Z[k+1]'*dV2[k]*Z[k+1]
        dense_mul_nn(n_,n_,n_,&Z_[((k+1) % K_)*n_*n_],&nnKa_[k].data()[0],&nnKb_[k].data()[0]);
      }
      
      // ********** STOP ****************
        
      std::fill(dX_.begin(),dX_.end(),0);

      // Main loops to loop over blocks of the block-upper triangular A
      // Outer main loop
      for (int l=0;l<partition_.size()-1;++l) {
        // F serves an an accumulator for intermediate summation results
        std::fill(F_.begin(),F_.end(),0);
        
        // ********** START ***************
        //for i in range(l):
        //  F[i] = [sum(mul(dX[i][j][k],A[l][j][k].T) for j in range(l)) for k in range(p) ]

        for (int i=0;i<l;++i) {
          int na1 = partition_[i+1]-partition_[i];
          int nb2 = partition_[l+1]-partition_[l];
          for (int j=0;j<l;++j) {
            int na2 = partition_[j+1]-partition_[j];
            for (int k=0;k<K_;++k) {
              for (int ii=0;ii<na1;++ii) {
                for (int jj=0;jj<nb2;++jj) {
                  for (int kk=0;kk<na2;++kk) {
                    F_[k*4*n_+4*i+2*ii+jj] += dX_[partindex(i,j,k,ii,kk)]*T_[partindex(l,j,k,jj,kk)];
                  }
                }
              }
            }
          }
        }
        // ********** STOP ***************

        // Inner main loop
        for (int r=0;r<l+1;++r) {
        
          // ********** START ***************
          // F[r] = [sum(mul(dX[r][j][k],A[l][j][k].T) for j in range(l)) for k in range(p) ]
          int na1 = partition_[r+1]-partition_[r];
          int nb2 = partition_[l+1]-partition_[l];

          for (int k=0;k<K_;++k) {
            std::fill(F_.begin()+k*4*n_+4*r,F_.begin()+k*4*n_+4*r+4,0);
          }
            
          for (int j=0;j<l;++j) {
            int na2 = partition_[j+1]-partition_[j];
            for (int k=0;k<K_;++k) {
              for (int ii=0;ii<na1;++ii) {
                for (int jj=0;jj<nb2;++jj) {
                  for (int kk=0;kk<na2;++kk) {
                    F_[k*4*n_+4*r+2*ii+jj] += dX_[partindex(r,j,k,ii,kk)]*T_[partindex(l,j,k,jj,kk)];
                  }
                }
              }
            }
          }
          // ********** STOP ***************
          
  
          // ********** START ***************
          // FF =   [sum(mul(A[r][i][k],dX[i][l][k]) for i in range(r)) for k in range(p)]
          // Each entry of FF is na1-by-na2
          std::fill(FF_.begin(),FF_.end(),0);
          { 
            int na1 = partition_[r+1]-partition_[r];
            for (int i=0;i<r;++i) {
              int nb2 = partition_[l+1]-partition_[l];
              int na2 = partition_[i+1]-partition_[i];
              for (int k=0;k<K_;++k) {
                for (int ii=0;ii<na1;++ii) {
                  for (int jj=0;jj<nb2;++jj) {
                    for (int kk=0;kk<na2;++kk) {
                      FF_[k*4+2*ii+jj] += T_[partindex(r,i,k,ii,kk)]*dX_[partindex(i,l,k,kk,jj)];
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
                M[np*((k+1)%K_)+ll*n2+m] = nnKb_[k].data()[(partition_[r]+ll)*n_ + partition_[l]+m];
              }
            }
          }

          // ********** START ***************
          // M+= [sum(mul(A[r][i][k],F[i][k])  for i in range(r+1)) for k in rang(p)]
          {
            int na1 = partition_[r+1]-partition_[r];
            int nb2 = partition_[l+1]-partition_[l];
            
            for (int i=0;i<r+1;++i) {
              int na2 = partition_[i+1]-partition_[i];
              for (int k=0;k<K_;++k) {
                for (int ii=0;ii<na1;++ii) {
                  for (int jj=0;jj<nb2;++jj) {
                    for (int kk=0;kk<na2;++kk) {
                      M[np*((k+1)%K_)+ii*n2+jj] += T_[partindex(r,i,k,ii,kk)]*F_[k*4*n_+4*i+2*kk+jj];
                    }
                  }
                }
              }
            }
          }
          // ********** STOP ***************
          
          // ********** START ***************
          // M+= [mul(FF[k],A[l][l][k].T) for k in rang(p)]
          {
            int na1 = partition_[r+1]-partition_[r];
            int na2 = partition_[l+1]-partition_[l];
            int nb2 = partition_[l+1]-partition_[l];
            for (int k=0;k<K_;++k) {
              for (int ii=0;ii<na1;++ii) {
                for (int jj=0;jj<nb2;++jj) {
                  for (int kk=0;kk<na2;++kk) {
                    M[np*((k+1)%K_)+ii*n2+jj] += FF_[k*4+2*ii+kk]*T_[partindex(l,l,k,jj,kk)];
                  }
                }
              }
            }
          }
          // ********** STOP ***************
          
          // Critical observation: Prepare step is not needed
          solver.solve(true);
          
          // Extract solution and store it in dX
          std::vector<double> & sol = solver.output().data();

          // ********** START ***************
          for (int ii=0;ii<partition_[r+1]-partition_[r];++ii) {
            for (int jj=0;jj<partition_[l+1]-partition_[l];++jj) {
              for (int k=0;k<K_;++k) {
                dX_[partindex(r,l,k,ii,jj)] = sol[n1*n2*k+n2*ii+jj];
              }
            }
          }
          
          for (int ii=0;ii<partition_[r+1]-partition_[r];++ii) {
            for (int jj=0;jj<partition_[l+1]-partition_[l];++jj) {
              for (int k=0;k<K_;++k) {
                dX_[partindex(l,r,k,jj,ii)] = sol[n1*n2*k+n2*ii+jj];
              }
            }
          }
          // ********** STOP ***************
          
        }
        
        output(DPLE_NUM_OUT*(d+1)+DPLE_P).set(0.0);
        
        for(int k=0;k<K_;++k) {
          nnKa_[k].set(0.0);
          
          // nnKa[k] <- V[k]*Z[k]'
          dense_mul_nn(n_,n_,n_,&dX_[k*n_*n_],&Z_[k*n_*n_],&nnKa_[k].data()[0]);
          // output <- Z[k]*V[k]*Z[k]'
          dense_mul_tn(n_,n_,n_,&Z_[k*n_*n_],&nnKa_[k].data()[0],&output(DPLE_NUM_OUT*(d+1)+DPLE_P).data()[k*n_*n_]);
        }
      }
  
    }
    
    for (int d=0;d<nadj_;++d) {
    
      DMatrix &P_bar = input(DPLE_NUM_IN*(nfwd_+1)+DPLE_NUM_OUT*d+DPLE_P);
      std::vector<double> &Vbar = output(DPLE_NUM_OUT*(nfwd_+1)+DPLE_NUM_IN*d+DPLE_V).data();
      std::fill(Vbar.begin(),Vbar.end(),0);
      
      std::vector<double> &Abar = output(DPLE_NUM_OUT*(nfwd_+1)+DPLE_NUM_IN*d+DPLE_A).data();
      std::fill(Abar.begin(),Abar.end(),0);
      
      std::fill(Xbar_.begin(),Xbar_.end(),0);
      
    
      // X_bar = [mul([Z[k].T, X_bar[k] , Z[k]]) for k in range(p)]
      for(int k=0;k<K_;++k) {
        nnKa_[k].set(0.0);
        
        // nnKa[k] <- nnKb*Z[k]
        dense_mul_nt(n_,n_,n_,&P_bar.data()[n_*n_*k],&Z_[k*n_*n_],&nnKa_[k].data()[0]);
        // Xbar <- Z[k]*V[k]*Z[k]'
        dense_mul_nn(n_,n_,n_,&Z_[k*n_*n_],&nnKa_[k].data()[0],&Xbar_[k*n_*n_]);
      }
      
      
      
      
      // Main loops to loop over blocks of the block-upper triangular A
      // Outer main loop
      for (int l=partition_.size()-2;l>=0;--l) {
      
        for (int k=0;k<K_;++k) {
          std::fill(F_.begin(),F_.end(),0);
        }
          
        // Inner main loop
        for (int r=l;r>=0;--r) {
        
          int n1 = partition_[r+1]-partition_[r];
          int n2 = partition_[l+1]-partition_[l];
          int np = n1*n2;
          
          LinearSolver & solver = dpse_solvers_[n1-1+n2-1][((l+1)*l)/2+r];
          
          std::vector<double> & B = solver.input(LINSOL_B).data();
          
          
          // ********** START ***************
          for (int ii=0;ii<partition_[r+1]-partition_[r];++ii) {
            for (int jj=0;jj<partition_[l+1]-partition_[l];++jj) {
              for (int k=0;k<K_;++k) {
                B[n1*n2*k+n2*ii+jj] = Xbar_[partindex(r,l,k,ii,jj)];
              }
            }
          }
          if (r!=l) {
            for (int ii=0;ii<partition_[r+1]-partition_[r];++ii) {
              for (int jj=0;jj<partition_[l+1]-partition_[l];++jj) {
                for (int k=0;k<K_;++k) {
                  B[n1*n2*k+n2*ii+jj]+=Xbar_[partindex(l,r,k,jj,ii)];
                }
              }
            }
          }
          // ********** STOP ***************
          
          solver.solve(false);
          
          
          // for k in range(p): V_bar[r][l][k]+=M_bar[k]
          std::vector<double> &Mbar = solver.output().data();
          
          for (int k=0;k<K_;++k) {
            for (int ll=0;ll<n1;++ll) {
              for (int m=0;m<n2;++m) {
                 Vbar[partindex(r,l,k,ll,m)]+= Mbar[np*((k+1)%K_)+ll*n2+m];
              }
            }
          }
          
          // ********** START ***************
          //     for k in range(p):
          //       FF_bar[k] += mul(M_bar[k],A[l][l][k])
          std::fill(FF_.begin(),FF_.end(),0);
          { 
            int na1 = partition_[r+1]-partition_[r];
            int na2 = partition_[l+1]-partition_[l];
            int nb2 = partition_[l+1]-partition_[l];
            for (int k=0;k<K_;++k) {
              for (int ii=0;ii<na1;++ii) {
                for (int jj=0;jj<nb2;++jj) {
                  for (int kk=0;kk<na2;++kk) {
                    FF_[k*4+2*ii+kk] += Mbar[np*((k+1)%K_)+ii*n2+jj]*T_[partindex(l,l,k,jj,kk)];
                  }
                }
              }
            }
          }
          
          // ********** STOP ***************
          
          
          // for k in range(p):
          //  for i in range(r):
          //    X_bar[i][l][k] +=mul(A[r][i][k].T,FF_bar[k])
          // ********** START ***************
          {
            int na1 = partition_[r+1]-partition_[r];
            for (int i=0;i<r;++i) {
              int nb2 = partition_[l+1]-partition_[l];
              int na2 = partition_[i+1]-partition_[i];
              for (int k=0;k<K_;++k) {
                for (int ii=0;ii<na1;++ii) {
                  for (int jj=0;jj<nb2;++jj) {
                    for (int kk=0;kk<na2;++kk) {
                      Xbar_[partindex(i,l,k,kk,jj)] += T_[partindex(r,i,k,ii,kk)]*FF_[k*4+2*ii+jj];
                    }
                  }
                }
              }
            }
          }
          // ********** STOP ***************

          //for k in range(p):
          //  for i in range(r+1):
          //    F_bar[i][k]+= mul(A[r][i][k].T,M_bar[k])
          // ********** START ***************
          {
            int na1 = partition_[r+1]-partition_[r];
            int nb2 = partition_[l+1]-partition_[l];
            
            for (int i=0;i<r+1;++i) {
              int na2 = partition_[i+1]-partition_[i];
              for (int k=0;k<K_;++k) {
                for (int ii=0;ii<na1;++ii) {
                  for (int jj=0;jj<nb2;++jj) {
                    for (int kk=0;kk<na2;++kk) {
                      F_[k*4*n_+4*i+2*kk+jj] += T_[partindex(r,i,k,ii,kk)]*Mbar[np*((k+1)%K_)+ii*n2+jj];
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
          //       X_bar[r][j][k] +=mul(F_bar[r][k],A[l][j][k]) 
          // ********** START ***************
          if (r>=l) {
            int na1 = partition_[r+1]-partition_[r];
            int nb2 = partition_[l+1]-partition_[l];
              
            for (int j=0;j<l;++j) {
              int na2 = partition_[j+1]-partition_[j];
              for (int k=0;k<K_;++k) {
                for (int ii=0;ii<na1;++ii) {
                  for (int jj=0;jj<nb2;++jj) {
                    for (int kk=0;kk<na2;++kk) {
                      Xbar_[partindex(r,j,k,ii,kk)] += F_[k*4*n_+4*r+2*ii+jj]*T_[partindex(l,j,k,jj,kk)];
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
        //     X_bar[i][j][k]+=mul(F_bar[i][k],A[l][j][k])
        // ********** START ***************
        //for i in range(l):
        //  F[i] = [sum(mul(X[i][j][k],A[l][j][k].T) for j in range(l)) for k in range(p) ]

        for (int i=0;i<l;++i) {
          int na1 = partition_[i+1]-partition_[i];
          int nb2 = partition_[l+1]-partition_[l];
          for (int j=0;j<l;++j) {
            int na2 = partition_[j+1]-partition_[j];
            for (int k=0;k<K_;++k) {
              for (int ii=0;ii<na1;++ii) {
                for (int jj=0;jj<nb2;++jj) {
                  for (int kk=0;kk<na2;++kk) {
                    Xbar_[partindex(i,j,k,ii,kk)] += F_[k*4*n_+4*i+2*ii+jj]*T_[partindex(l,j,k,jj,kk)];
                  }
                }
              }
            }
          }
        }
        // ********** STOP ***************
        
      }
      
      
      // V_bar = [mul([sZ[k], V_bar[k] , sZ[k].T]) for k in range(p)]
      for(int k=0;k<K_;++k) {
        nnKa_[k].set(0.0);
        dense_mul_nn(n_,n_,n_,&Vbar[n_*n_*k],&Z_[((k+1)%K_)*n_*n_],&nnKa_[k].data()[0]);
        std::fill(&Vbar[n_*n_*k],&Vbar[n_*n_*(k+1)],0);
        dense_mul_tn(n_,n_,n_,&Z_[((k+1)%K_)*n_*n_],&nnKa_[k].data()[0],&Vbar[n_*n_*k]);
      }
      
      // Force symmetry
      for(int k=0;k<K_;++k) {
        for (int r=0;r<n_;++r) {
          for (int l=0;l<r;++l) {
            double s = (Vbar[n_*n_*k+r+l*n_]+Vbar[n_*n_*k+r*n_+l])/2;
            Vbar[n_*n_*k+r+l*n_] = s;
            Vbar[n_*n_*k+r*n_+l] = s;
          }
        }
      }
      
      // A_bar = [mul([vb+vb.T,a,x]) for vb,x,a in zip(V_bar,X,As)]
      for(int k=0;k<K_;++k) {
        std::fill(nnKa_[k].begin(),nnKa_[k].end(),0);
        dense_mul_nn(n_,n_,n_,&input(DPLE_A).data()[n_*n_*k],&output(DPLE_P).data()[n_*n_*k],&nnKa_[k].data()[0]);

        dense_mul_nn(n_,n_,n_,&Vbar[n_*n_*k],&nnKa_[k].data()[0],&Abar[n_*n_*k]);
        dense_mul_tn(n_,n_,n_,&Vbar[n_*n_*k],&nnKa_[k].data()[0],&Abar[n_*n_*k]);
      }
      
      std::fill(P_bar.data().begin(),P_bar.data().end(),0);
 
    }
    
  }
  
  FX PsdIndefDpleInternal::getDerivative(int nfwd, int nadj) {
    casadi_assert(nfwd_==0 && nadj_==0);
    
    PsdIndefDpleInternal* node = new PsdIndefDpleInternal(A_,V_, nfwd, nadj);
    node->setOption(dictionary());
    
    PsdIndefDpleSolver ret;
    ret.assignNode(node);

    return ret;
    
  }

  void PsdIndefDpleInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
    DpleInternal::deepCopyMembers(already_copied);
  }
  
  PsdIndefDpleInternal* PsdIndefDpleInternal::clone() const{
    // Return a deep copy
    PsdIndefDpleInternal* node = new PsdIndefDpleInternal(A_,V_, nfwd_, nadj_);
    node->setOption(dictionary());
    return node;
  }


} // namespace CasADi


