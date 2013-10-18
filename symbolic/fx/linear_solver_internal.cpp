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

INPUTSCHEME(LinsolInput)
OUTPUTSCHEME(LinsolOutput)

using namespace std;
namespace CasADi{

  LinearSolverInternal::LinearSolverInternal(const CRSSparsity& sparsity, int nrhs){
    // No OO derivatives supported/needed
    setOption("number_of_fwd_dir",0);
    setOption("number_of_adj_dir",0);
    setOption("max_number_of_fwd_dir",0);
    setOption("max_number_of_adj_dir",0);

    // Make sure arguments are consistent
    casadi_assert(!sparsity.isNull());
    casadi_assert_message(sparsity.size1()==sparsity.size2(),"LinearSolverInternal::init: the matrix must be square but got " << sparsity.dimString());  
    casadi_assert_message(!isSingular(sparsity),"LinearSolverInternal::init: singularity - the matrix is structurally rank-deficient. sprank(J)=" << rank(sparsity) << " (in stead of "<< sparsity.size1() << ")");

    // Allocate inputs
    input_.resize(LINSOL_NUM_IN);
    input(LINSOL_A) = DMatrix(sparsity);
    input(LINSOL_B) = DMatrix(nrhs,sparsity.size1(),0);
    input(LINSOL_T) = 0.;
  
    // Allocate outputs
    output_.resize(LINSOL_NUM_OUT);
    output(LINSOL_X) = input(LINSOL_B);

    inputScheme_ = SCHEME_LinsolInput;
    outputScheme_ = SCHEME_LinsolOutput;
  }

  void LinearSolverInternal::init(){
    // Call the base class initializer
    FXInternal::init();
    
    // Not prepared
    prepared_ = false;
  }

  LinearSolverInternal::~LinearSolverInternal(){
  }
 
  void LinearSolverInternal::evaluate(int nfdir, int nadir){
    casadi_assert_message(nfdir==0 && nadir==0,"Directional derivatives for LinearSolver not supported. Reformulate or wrap in an MXFunction instance.");

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
    const vector<double>& b = input(LINSOL_B).data();
    vector<double>& x = output(LINSOL_X).data();
    bool transpose = input(LINSOL_T).toScalar()!=0.;
    int nrhs = input(LINSOL_B).size1();

    // Copy input to output
    copy(b.begin(),b.end(),x.begin());
  
    // Solve the factorized system in-place
    solve(getPtr(x),nrhs,transpose);
  }
 
} // namespace CasADi
 
