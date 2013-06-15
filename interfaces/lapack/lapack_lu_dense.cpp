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
#include "../../symbolic/stl_vector_tools.hpp"

using namespace std;
namespace CasADi{

  LapackLUDense::LapackLUDense(){
  }

  LapackLUDense::LapackLUDense(const CRSSparsity& sparsity, int nrhs){
    assignNode(new LapackLUDenseInternal(sparsity,nrhs));
  }
 
  LapackLUDenseInternal* LapackLUDense::operator->(){
    return static_cast<LapackLUDenseInternal*>(FX::operator->());
  }

  const LapackLUDenseInternal* LapackLUDense::operator->() const{
    return static_cast<const LapackLUDenseInternal*>(FX::operator->());
  }

  LapackLUDenseInternal::LapackLUDenseInternal(const CRSSparsity& sparsity, int nrhs) : LinearSolverInternal(sparsity,nrhs){
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
    nrow_ = nrow();
    ncol_ = ncol();
  
    // Currently only square matrices tested
    if(nrow_!=ncol_) throw CasadiException("LapackLUDenseInternal::LapackLUDenseInternal: currently only square matrices implemented.");

    // Allocate matrix
    mat_.resize(nrow_*nrow_);
    ipiv_.resize(nrow_);
  
    // Equilibriate?
    equilibriate_ = getOption("equilibration").toInt();
    if(equilibriate_){
      r_.resize(nrow_);
      c_.resize(ncol_);
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
      // Calculate the row and column scaling factors
      double rowcnd, colcnd; // ratio of smallest to largest row/column scaling factor
      double amax; // absolute value of the largest matrix element
      int info = -100;
      dgeequ_(&nrow_,&ncol_,getPtr(mat_),&nrow_,getPtr(r_),getPtr(c_),&rowcnd, &colcnd, &amax, &info);
      if(info < 0) throw CasadiException("LapackQRDenseInternal::prepare: dgeequ_ failed to calculate the scaling factors");
      if(info>0){
        stringstream ss;
        ss << "LapackLUDenseInternal::prepare: ";
        if(info<=nrow_)  ss << (info-1) << "-th column (zero-based) is exactly zero";
        else             ss << (info-1-nrow_) << "-th row (zero-based) is exactly zero";

        cout << "Warning: " << ss.str() << endl;


      
        if(allow_equilibration_failure_)  cout << "Warning: " << ss.str() << endl;
        else                              casadi_error(ss.str());
      }
  
      // Equilibriate the matrix if scaling was successful
      if(info!=0)
        dlaqge_(&nrow_,&ncol_,getPtr(mat_),&nrow_,getPtr(r_),getPtr(c_),&rowcnd, &colcnd, &amax, &equed_);
      else
        equed_ = 'N';
    }
  
    // Factorize the matrix
    int info = -100;
    dgetrf_(&nrow_, &nrow_, getPtr(mat_), &nrow_, getPtr(ipiv_), &info);
    if(info != 0) throw CasadiException("LapackLUDenseInternal::prepare: dgetrf_ failed to factorize the jacobian");
  
    // Sucess if reached this point
    prepared_ = true;
  }
    
  void LapackLUDenseInternal::solve(double* x, int nrhs, bool transpose){
    // Scale the right hand side
    if(transpose){
      columnScaling(x,nrhs);
    } else {
      rowScaling(x,nrhs);
    }

    // Solve the system of equations
    int info = 100;
    char trans = transpose ? 'T' : 'N';
    dgetrs_(&trans, &nrow_, &nrhs, getPtr(mat_), &nrow_, getPtr(ipiv_), x, &nrow_, &info);
    if(info != 0) throw CasadiException("LapackLUDenseInternal::solve: failed to solve the linear system");

    // Scale the solution
    if(transpose){
      rowScaling(x,nrhs);
    } else {
      columnScaling(x,nrhs);
    }
  }

  void LapackLUDenseInternal::rowScaling(double* x, int nrhs){
    // Scale result if this was done to the matrix
    if(equed_=='R' || equed_=='B')
      for(int rhs=0; rhs<nrhs; ++rhs)
        for(int i=0; i<nrow_; ++i)
          x[i+rhs*ncol_] *= r_[i];
  }
    
  void LapackLUDenseInternal::columnScaling(double* x, int nrhs){
    // Scale right hand side if this was done to the matrix
    if(equed_=='C' || equed_=='B')
      for(int rhs=0; rhs<nrhs; ++rhs)
        for(int i=0; i<ncol_; ++i)
          x[i+rhs*ncol_] *= c_[i];
  }

  LapackLUDenseInternal* LapackLUDenseInternal::clone() const{
    return new LapackLUDenseInternal(*this);
  }

} // namespace CasADi
