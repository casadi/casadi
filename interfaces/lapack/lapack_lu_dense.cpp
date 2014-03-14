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

#include "lapack_lu_dense.hpp"
#include "../../symbolic/std_vector_tools.hpp"

using namespace std;
namespace CasADi{

  LapackLUDense::LapackLUDense(){
  }

  LapackLUDense::LapackLUDense(const Sparsity& sparsity, int nrhs){
    assignNode(new LapackLUDenseInternal(sparsity,nrhs));
  }
 
  LapackLUDenseInternal* LapackLUDense::operator->(){
    return static_cast<LapackLUDenseInternal*>(FX::operator->());
  }

  const LapackLUDenseInternal* LapackLUDense::operator->() const{
    return static_cast<const LapackLUDenseInternal*>(FX::operator->());
  }

  LapackLUDenseInternal::LapackLUDenseInternal(const Sparsity& sparsity, int nrhs) : LinearSolverInternal(sparsity,nrhs){
    // Equilibriate the matrix
    addOption("equilibration",OT_BOOLEAN,true);
    addOption("allow_equilibration_failure",OT_BOOLEAN,false);
  }

  LapackLUDenseInternal::~LapackLUDenseInternal(){
  }

  void LapackLUDenseInternal::init(){
    // Call the base class initializer
    LinearSolverInternal::init();
  
    // Get dimensions
    ncol_ = ncol();
    nrow_ = nrow();
  
    // Currently only square matrices tested
    if(ncol_!=nrow_) throw CasadiException("LapackLUDenseInternal::LapackLUDenseInternal: currently only square matrices implemented.");

    // Allocate matrix
    mat_.resize(ncol_*ncol_);
    ipiv_.resize(ncol_);
  
    // Equilibriate?
    equilibriate_ = getOption("equilibration").toInt();
    if(equilibriate_){
      r_.resize(ncol_);
      c_.resize(nrow_);
    }
    equed_ = 'N'; // No equilibration

    // Allow equilibration failures
    allow_equilibration_failure_ = getOption("allow_equilibration_failure").toInt();
  }

  void LapackLUDenseInternal::prepare(){
    prepared_ = false;
  
    // Get the elements of the matrix, dense format
    input(0).get(mat_,DENSE);

    if(equilibriate_){
      // Calculate the col and row scaling factors
      double colcnd, rowcnd; // ratio of smallest to largest col/row scaling factor
      double amax; // absolute value of the largest matrix element
      int info = -100;
      dgeequ_(&ncol_,&nrow_,getPtr(mat_),&ncol_,getPtr(r_),getPtr(c_),&colcnd, &rowcnd, &amax, &info);
      if(info < 0) throw CasadiException("LapackQRDenseInternal::prepare: dgeequ_ failed to calculate the scaling factors");
      if(info>0){
        stringstream ss;
        ss << "LapackLUDenseInternal::prepare: ";
        if(info<=ncol_)  ss << (info-1) << "-th row (zero-based) is exactly zero";
        else             ss << (info-1-ncol_) << "-th col (zero-based) is exactly zero";

        cout << "Warning: " << ss.str() << endl;


      
        if(allow_equilibration_failure_)  cout << "Warning: " << ss.str() << endl;
        else                              casadi_error(ss.str());
      }
  
      // Equilibriate the matrix if scaling was successful
      if(info!=0)
        dlaqge_(&ncol_,&nrow_,getPtr(mat_),&ncol_,getPtr(r_),getPtr(c_),&colcnd, &rowcnd, &amax, &equed_);
      else
        equed_ = 'N';
    }
  
    // Factorize the matrix
    int info = -100;
    dgetrf_(&ncol_, &ncol_, getPtr(mat_), &ncol_, getPtr(ipiv_), &info);
    if(info != 0) throw CasadiException("LapackLUDenseInternal::prepare: dgetrf_ failed to factorize the jacobian");
  
    // Sucess if reached this point
    prepared_ = true;
  }
    
  void LapackLUDenseInternal::solve(double* x, int nrhs, bool transpose){
    // Scale the right hand side
    if(transpose){
      rowScaling(x,nrhs);
    } else {
      colScaling(x,nrhs);
    }

    // Solve the system of equations
    int info = 100;
    char trans = transpose ? 'T' : 'N';
    dgetrs_(&trans, &ncol_, &nrhs, getPtr(mat_), &ncol_, getPtr(ipiv_), x, &ncol_, &info);
    if(info != 0) throw CasadiException("LapackLUDenseInternal::solve: failed to solve the linear system");

    // Scale the solution
    if(transpose){
      colScaling(x,nrhs);
    } else {
      rowScaling(x,nrhs);
    }
  }

  void LapackLUDenseInternal::colScaling(double* x, int nrhs){
    // Scale result if this was done to the matrix
    if(equed_=='R' || equed_=='B')
      for(int rhs=0; rhs<nrhs; ++rhs)
        for(int i=0; i<ncol_; ++i)
          x[i+rhs*nrow_] *= r_[i];
  }
    
  void LapackLUDenseInternal::rowScaling(double* x, int nrhs){
    // Scale right hand side if this was done to the matrix
    if(equed_=='C' || equed_=='B')
      for(int rhs=0; rhs<nrhs; ++rhs)
        for(int i=0; i<nrow_; ++i)
          x[i+rhs*nrow_] *= c_[i];
  }

  LapackLUDenseInternal* LapackLUDenseInternal::clone() const{
    return new LapackLUDenseInternal(*this);
  }

} // namespace CasADi
