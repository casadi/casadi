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

#include "implicit_function_internal.hpp"

using namespace std;
namespace CasADi{

  ImplicitFunctionInternal::ImplicitFunctionInternal(const FX& f) : f_(f){
    addOption("linear_solver",            OT_LINEARSOLVER, GenericType(), "User-defined linear solver class. Needed for sensitivities.");
    addOption("linear_solver_options",    OT_DICTIONARY,   GenericType(), "Options to be passed to the linear solver.");
  }

  void ImplicitFunctionInternal::init(){
    // Initialize the residual function
    if(!f_.isInit()) f_.init();
  
    // Allocate inputs
    setNumInputs(f_.getNumInputs()-1);
    for(int i=0; i<getNumInputs(); ++i){
      input(i) = f_.input(i+1);
    }
  
    // Allocate outputs
    setNumOutputs(f_.getNumOutputs());
    output(0) = f_.input(0);
    for(int i=1; i<getNumOutputs(); ++i){
      output(i) = f_.output(i);
    }

    // Call the base class initializer
    FXInternal::init();

    // Number of equations
    N_ = output().size();

    // Generate Jacobian if not provided
    if(J_.isNull()) J_ = f_.jacobian(0,0);
    J_.init();
  
    casadi_assert_message(J_.output().size1()==J_.output().size2(),"ImplicitFunctionInternal::init: the jacobian must be square but got " << J_.output().dimString());
  
    casadi_assert_message(!isSingular(J_.output().sparsity()),"ImplicitFunctionInternal::init: singularity - the jacobian is structurally rank-deficient. sprank(J)=" << sprank(J_.output()) << " (in stead of "<< J_.output().size1() << ")");
  
    // Get the linear solver creator function
    if(linsol_.isNull() && hasSetOption("linear_solver")){
      linearSolverCreator linear_solver_creator = getOption("linear_solver");
  
      // Allocate an NLP solver
      linsol_ = linear_solver_creator(CRSSparsity());
  
      // Pass options
      if(hasSetOption("linear_solver_options")){
	const Dictionary& linear_solver_options = getOption("linear_solver_options");
	linsol_.setOption(linear_solver_options);
      }
    }
  
    // Initialize the linear solver, if provided
    if(!linsol_.isNull()){
      linsol_.setSparsity(J_.output().sparsity());
      linsol_.init();
    }
    
    // Allocate memory for directional derivatives
    ImplicitFunctionInternal::updateNumSens(false);

  }

  void ImplicitFunctionInternal::updateNumSens(bool recursive){
    // Call the base class if needed
    if(recursive) FXInternal::updateNumSens(recursive);
  
    // Request more directional derivatives for the residual function
    f_.requestNumSens(nfdir_,nadir_);
  }

  void ImplicitFunctionInternal::evaluate_sens(int nfdir, int nadir, bool linsol_prepared) {
    // Make sure that a linear solver has been provided
    casadi_assert_message(!linsol_.isNull(),"Sensitivities of an implicit function requires a provided linear solver");
    casadi_assert_message(!J_.isNull(),"Sensitivities of an implicit function requires an exact Jacobian");
  
  
    // General scheme:  f(z,x_i) = 0
    //
    //  Forward sensitivities:
    //     dot(f(z,x_i)) = 0
    //     df/dz dot(z) + Sum_i df/dx_i dot(x_i) = 0
    //
    //     dot(z) = [df/dz]^(-1) [ Sum_i df/dx_i dot(x_i) ] 
    //
    //     dot(y_i) = dy_i/dz dot(z) + Sum_j dy_i/dx_i dot(x_i)
    //
    //  Adjoint sensitivitites:

    if (!linsol_prepared) {
      // Pass inputs
      J_.setInput(output(),0);
      for(int i=0; i<getNumInputs(); ++i)
	J_.setInput(input(i),i+1);

      // Evaluate jacobian
      J_.evaluate();

      // Pass non-zero elements, scaled by -gamma, to the linear solver
      linsol_.setInput(J_.output(),0);

      // Prepare the solution of the linear system (e.g. factorize)
      linsol_.prepare();
    }
  
    // Pass inputs to function
    f_.setInput(output(0),0);
    for(int i=0; i<getNumInputs(); ++i)
      f_.input(i+1).set(input(i));

    // Pass input seeds to function
    for(int dir=0; dir<nfdir; ++dir){
      f_.fwdSeed(0,dir).setZero();
      for(int i=0; i<getNumInputs(); ++i){
	f_.fwdSeed(i+1,dir).set(fwdSeed(i,dir));
      }
    }
  
    // Solve for the adjoint seeds
    for(int dir=0; dir<nadir; ++dir){
      // Negate adjoint seed and pass to function
      Matrix<double>& faseed = f_.adjSeed(0,dir);
      faseed.set(adjSeed(0,dir));
      for(vector<double>::iterator it=faseed.begin(); it!=faseed.end(); ++it){
	*it = -*it;
      }
    
      // Solve the transposed linear system
      linsol_.solve(&faseed.front(),1,true);

      // Set auxillary adjoint seeds
      for(int oind=1; oind<getNumOutputs(); ++oind){
	f_.adjSeed(oind,dir).set(adjSeed(oind,dir));
      }
    }
  
    // Evaluate
    f_.evaluate(nfdir,nadir);
  
    // Get the forward sensitivities
    for(int dir=0; dir<nfdir; ++dir){
      // Negate intermediate result and copy to output
      Matrix<double>& fsens = fwdSens(0,dir);
      fsens.set(f_.fwdSens(0,dir));
      for(vector<double>::iterator it=fsens.begin(); it!=fsens.end(); ++it){
	*it = -*it;
      }
    
      // Solve the linear system
      linsol_.solve(&fsens.front());
    }
  
    // Get auxillary forward sensitivities
    if(getNumOutputs()>1){
      // Pass the seeds to the implicitly defined variables
      for(int dir=0; dir<nfdir; ++dir){
	f_.fwdSeed(0,dir).set(fwdSens(0,dir));
      }
    
      // Evaluate
      f_.evaluate(nfdir);
  
      // Get the sensitivities
      for(int dir=0; dir<nfdir; ++dir){
	for(int oind=1; oind<getNumOutputs(); ++oind){
	  fwdSens(oind,dir).set(f_.fwdSens(oind,dir));
	}
      }
    }
  
    // Get the adjoint sensitivities
    for(int dir=0; dir<nadir; ++dir){
      for(int i=0; i<getNumInputs(); ++i){
	f_.adjSens(i+1,dir).get(adjSens(i,dir));
      }
    }
  }

  ImplicitFunctionInternal::~ImplicitFunctionInternal(){
  }

  void ImplicitFunctionInternal::setJacobian(FX &J) {
    J_ = J;
  }
 
 
} // namespace CasADi

  


