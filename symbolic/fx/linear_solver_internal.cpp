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

#include "linear_solver_internal.hpp"
#include "../stl_vector_tools.hpp"

using namespace std;
namespace CasADi{

LinearSolverInternal::LinearSolverInternal(const CRSSparsity& sparsity) : sparsity_(sparsity){
  casadi_assert(!sparsity.isNull());
  addOption("trans", OT_BOOLEAN, false);
}

void LinearSolverInternal::init(){
  // Transpose?
  transpose_ = getOption("trans");
  
  casadi_assert_message(sparsity_.size1()==sparsity_.size2(),"LinearSolverInternal::init: the matrix must be square but got " << sparsity_.dimString());
  
  casadi_assert_message(!isSingular(sparsity_),"LinearSolverInternal::init: singularity - the matrix is structurally rank-deficient. sprank(J)=" << rank(sparsity_) << " (in stead of "<< sparsity_.size1() << ")");
  
  // Allocate space for inputs
  input_.resize(2);
  input(0) = DMatrix(sparsity_);
  input(1) = DMatrix(sparsity_.size1(),1,0);
  
  // Allocate space for outputs
  output_.resize(1);
  output(0) = input(1);
  
  // Not prepared
  prepared_ = false;

  // Call the base class initializer
  FXInternal::init();
}

LinearSolverInternal::~LinearSolverInternal(){
}
 
void LinearSolverInternal::evaluate(int nfdir, int nadir){
/*  Factorization fact;
  if(called_once){
    // Check if any element has changed
    bool any_change = false;
    const vector<double>& val = input(0).data();
    for(int i=0; i<val.size(); ++i){
      if(val[i] != a[i]){
        any_change = true;
        break;
      }
    }
    
    // Reuse factored matrix if matrix hasn't changed
    fact = any_change ? SAMEPATTERN : FACTORED;
  } else {
    fact = DOFACT;
    called_once = true;
  }*/
  
  // Call the solve routine
  prepare();
  
  // Make sure preparation successful
  if(!prepared_) 
    throw CasadiException("LinearSolverInternal::evaluate: Preparation failed");
  
  // Solve the factorized system
  solve();
}
 
void LinearSolverInternal::solve(){
  // Get input and output vector
  const vector<double>& b = input(1).data();
  vector<double>& x = output().data();
  
  // Copy input to output
  copy(b.begin(),b.end(),x.begin());
  
  // Solve the factorized system
  solve(getPtr(x),1,transpose_);
}
 
} // namespace CasADi

  


