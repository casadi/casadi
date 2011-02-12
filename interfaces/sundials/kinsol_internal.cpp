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

#include "kinsol_internal.hpp"
#include "casadi/fx/sx_function_internal.hpp"
#include "casadi/stl_vector_tools.hpp"
#include "casadi/sx/sx_tools.hpp"
#include "casadi/fx/linear_solver_internal.hpp"

using namespace std;
namespace CasADi{
namespace Sundials{

KinsolInternal::KinsolInternal(const FX& f, int nrhs) : ImplicitFunctionInternal(f,nrhs){
  addOption("linear_solver", OT_STRING, "dense");
  addOption("upper_bandwidth", OT_INTEGER);
  addOption("lower_bandwidth", OT_INTEGER);
  addOption("max_krylov", OT_INTEGER, 0);
  addOption("exact_jacobian", OT_BOOLEAN, true);

  mem_ = 0;
  u_ = 0;
  u_scale_ = 0;
  f_scale_ = 0;
}

KinsolInternal* KinsolInternal::clone() const{
  return new KinsolInternal(*this);
}

KinsolInternal::~KinsolInternal(){
  if(u_) N_VDestroy_Serial(u_);
  if(u_scale_) N_VDestroy_Serial(u_scale_);
  if(f_scale_) N_VDestroy_Serial(f_scale_);
  if(mem_) KINFree(&mem_);
}

void KinsolInternal::init(){
  ImplicitFunctionInternal::init();
  
  // Read options
  strategy_ = KIN_LINESEARCH;
  
  // Return flag
  int flag;
  
  // Generate Jacobian if not provided
  if(J_.isNull()) J_ = f_.jacobian(0,0);
  J_.init();
  
  // Allocate N_Vectors
  u_ = N_VMake_Serial(N_,&output(0)[0]);
  u_scale_ = N_VNew_Serial(N_);
  f_scale_ = N_VNew_Serial(N_);
  
  // Set scaling factors to one
  N_VConst(1.0,u_scale_);
  N_VConst(1.0,f_scale_);
  
  // Create KINSOL memory block
  mem_ = KINCreate();
  
  // Set optional inputs
  flag = KINSetUserData(mem_, this);
  casadi_assert_message(flag==KIN_SUCCESS, "KINSetUserData");
  
  // Initialize KINSOL
  KINInit(mem_,func_wrapper, u_);

  // Use exact Jacobian
  bool exact_jacobian = getOption("exact_jacobian").toInt();
  
  // attach a linear solver
  if(getOption("linear_solver")=="dense"){
    // Dense jacobian
    flag = KINDense(mem_, N_);
    casadi_assert_message(flag==KIN_SUCCESS, "KINDense");
    
    if(exact_jacobian){
      flag = KINDlsSetDenseJacFn(mem_, djac_wrapper);
      casadi_assert_message(flag==KIN_SUCCESS, "KINDlsSetDenseJacFn");
    }
    
  } else if(getOption("linear_solver")=="banded") {
    // Banded jacobian
    flag = KINBand(mem_, N_, getOption("upper_bandwidth").toInt(), getOption("lower_bandwidth").toInt());
    casadi_assert_message(flag==KIN_SUCCESS, "KINBand");
  } else if(getOption("linear_solver")=="iterative") {
    // Sparse (iterative) solver  

    // Max dimension of the Krylov space
    int maxl = getOption("max_krylov").toInt();

    // Attach the sparse solver  
    if(getOption("iterative_solver")=="gmres"){
      flag = KINSpgmr(mem_, maxl);
      casadi_assert_message(flag==KIN_SUCCESS, "KINSpgmr");
    } else if(getOption("iterative_solver")=="bcgstab") {
      flag = KINSpbcg(mem_, maxl);
      casadi_assert_message(flag==KIN_SUCCESS, "KINSpbcg");
    } else if(getOption("iterative_solver")=="tfqmr") {
      flag = KINSptfqmr(mem_, maxl);
      casadi_assert_message(flag==KIN_SUCCESS, "KINSptfqmr");
    } else {
      throw CasadiException("KINSOL: Unknown sparse solver");
    }
  } else if(getOption("linear_solver")=="user_defined") {
    casadi_assert_message(0,"not implemented");
  } else {
    throw CasadiException("Unknown linear solver ");
  }
  
  // Sensitivities
  if(nfdir_>0 || nadir_>0){
    aug_ = jac(range(getNumInputs()),0);
    aug_.setOption("number_of_fwd_dir",0);
    aug_.setOption("number_of_adj_dir",0);
    aug_.init();
  }
  
}


void KinsolInternal::evaluate(int fsens_order, int asens_order){
  // Reset the counters
  t_func_ = 0;
  t_jac_ = 0;

  if(fsens_order==0 && asens_order==0){
    // Solve the nonlinear system of equations
    int flag = KINSol(mem_, u_, strategy_, u_scale_, f_scale_);
    casadi_assert_message(flag>=KIN_SUCCESS, "KINSol");
  } else {
    // Pass inputs to the augmented system solver
    for(int i=0; i<getNumInputs(); ++i)
      aug_.setInput(input(i),i);
    
    // Call the nonlinear equation solver for the augmented system
    aug_.evaluate();
    
    // Pointer to the augmented system data
    double *augdata = &aug_.output()[0];

    // Number of equations
    int nz = output().size();
    
    // Save the results
    output().set(augdata);

    // Loop over forward directions
    if(fsens_order){
      for(int dir=0; dir<nfdir_; ++dir){
        // Pointer to jacobian data
        double *J = augdata + nz;

        // Get the result of the matrix vector multiplication
        Matrix<double>& r = fwdSens(0,dir);
        r.setZero();
        
        for(int iind=0; iind<getNumInputs(); ++iind){
          // Get the vector to multiply the matrix with
          const Matrix<double>& v = fwdSeed(iind,dir);
        
          // Number of variables
          int nx = input(iind).numel();

          // Multiply matrix from the right
          for(int i=0; i<nz; ++i)
            for(int j=0; j<nx; ++j)
              r[i] += J[i+nz*j]*v[j];
            
          // Point to the next matrix
          J += nz*nx;
        }
      }
    }

    // Reset the adjoint sensitivities
    if(asens_order){
      for(int dir=0; dir<nadir_; ++dir){
        // Pointer to jacobian data
        double *J = augdata + nz;
        
        // Get the vector to multiply the matrix with
        const Matrix<double>& v = adjSeed(0,dir);
        
        for(int iind=0; iind<getNumInputs(); ++iind){
          // Get the result of the matrix vector multiplication
          Matrix<double>& r = adjSens(iind,dir);
          r.setZero();
          
          // Number of variables
          int nx = input(iind).numel();

          // Multiply matrix from the left
          for(int i=0; i<nz; ++i)
            for(int j=0; j<nx; ++j)
              r[j] += J[i+nz*j]*v[i];
            
          // Point to the next matrix
          J += nz*nx;
        }
      }
    }
  }
}

KinsolSolver KinsolInternal::jac(int iind, int oind){
  return jac(vector<int>(1,iind),oind);
}

KinsolSolver KinsolInternal::jac(const std::vector<int> iind, int oind){
  // Single output
  casadi_assert(oind==0);
  
  // Get the function
  SXFunction f = shared_cast<SXFunction>(f_);
  casadi_assert(!f.isNull());
  
  // Get the jacobians
  Matrix<SX> Jz = f.jac(0,0);

  // Number of equations
  int nz = f.input(0).numel();

  // All variables
  vector<Matrix<SX> > f_in(f.getNumInputs());
  f_in[0] = f.inputSX(0);
  for(int i=1; i<f.getNumInputs(); ++i)
    f_in[i] = f.inputSX(i);

  // Augmented nonlinear equation
  Matrix<SX> F_aug = f.outputSX();

  // Number of right hand sides
  int nrhs = 1;
  
  // Augment variables and equations
  for(vector<int>::const_iterator it=iind.begin(); it!=iind.end(); ++it){

    // Get the jacobian
    Matrix<SX> Jx = f.jac(*it+1,0);
    
    // Number of variables
    int nx = f.input(*it+1).numel();

    // Sensitivities
    Matrix<SX> dz_dx = symbolic("dz_dx", nz, nx);
    
    // Derivative equation
    Matrix<SX> f_der = prod(Jz,dz_dx) + Jx;

    // Append variables
    append(f_in[0],vec(dz_dx));
      
    // Augment nonlinear equation
    append(F_aug,vec(f_der));
    
    // Number of right hand sides
    nrhs += nx;
  }

  // Augmented nonlinear equation
  SXFunction f_aug(f_in, F_aug);

  // Create KINSOL instance
  KinsolSolver augsol(f_aug, nrhs);
  
  // Return the created solver
  return augsol;
}

void KinsolInternal::func(N_Vector u, N_Vector fval){
  // Get time
  time1_ = clock();

  // Pass input
  f_.setInput(NV_DATA_S(u),0);
  for(int i=0; i<getNumInputs(); ++i)
    f_.setInput(input(i),i+1);
  
  // Evaluate
  f_.evaluate();
    
    // Get results
  f_.getOutput(NV_DATA_S(fval));

  // Log time
  time2_ = clock();
  t_func_ += double(time2_-time1_)/CLOCKS_PER_SEC;
}

int KinsolInternal::func_wrapper(N_Vector u, N_Vector fval, void *user_data){
try{
    casadi_assert(user_data);
    KinsolInternal *this_ = (KinsolInternal*)user_data;
    this_->func(u,fval);
    return 0;
  } catch(exception& e){
    cerr << "func failed: " << e.what() << endl;
    return 1;
  }
}

int KinsolInternal::djac_wrapper(int N, N_Vector u, N_Vector fu, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2){
  try{
    casadi_assert(user_data);
    KinsolInternal *this_ = (KinsolInternal*)user_data;
    this_->djac(N, u, fu, J, tmp1, tmp2);
    return 0;
  } catch(exception& e){
    cerr << "djac failed: " << e.what() << endl;;
    return 1;
  }
}

void KinsolInternal::djac(int N, N_Vector u, N_Vector fu, DlsMat J, N_Vector tmp1, N_Vector tmp2){
  // Get time
  time1_ = clock();

  // Pass inputs to the jacobian function
  f_.setInput(NV_DATA_S(u),0);
  for(int i=0; i<getNumInputs(); ++i)
    f_.setInput(input(i),i+1);

  // Evaluate
  J_.evaluate();
  
  // Get sparsity and non-zero elements
  const vector<int>& rowind = J_.output().rowind();
  const vector<int>& col = J_.output().col();
  const vector<double>& val = J_.output();

  // Loop over rows
  for(int i=0; i<rowind.size()-1; ++i){
    // Loop over non-zero entries
    for(int el=rowind[i]; el<rowind[i+1]; ++el){
      // Get column
      int j = col[el];
      
      // Set the element
      DENSE_ELEM(J,i,j) = val[el];
    }
  }
  
  // Log time duration
  time2_ = clock();
  t_jac_ += double(time2_-time1_)/CLOCKS_PER_SEC;
}







} // namespace Sundials
} // namespace CasADi

