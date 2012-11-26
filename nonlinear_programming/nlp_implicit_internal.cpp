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

#include "nlp_implicit_internal.hpp"

#include "symbolic/mx/mx_tools.hpp"
#include "symbolic/fx/mx_function.hpp"

using namespace std;
namespace CasADi {

NLPImplicitInternal* NLPImplicitInternal::clone() const{
  // Return a deep copy
  NLPImplicitInternal* node = new NLPImplicitInternal(f_,nrhs_);
  if(!node->is_init_)
    node->init();
  return node;
}
  
NLPImplicitInternal::NLPImplicitInternal(const FX& f, int nrhs) : ImplicitFunctionInternal(f,nrhs) {

  addOption("nlp_solver",       OT_NLPSOLVER, GenericType(), "The NLPSolver used to solve the implicit system.");
  addOption("nlp_solver_options",       OT_DICTIONARY, GenericType(), "Options to be passed to the NLPSolver");
  addOption("linear_solver",    OT_LINEARSOLVER, GenericType(), "User-defined linear solver class. Needed for sensitivities.");
  addOption("linear_solver_options",    OT_DICTIONARY, GenericType(), "Options to be passed to the linear solver.");

}

NLPImplicitInternal::~NLPImplicitInternal(){ 
}

void NLPImplicitInternal::evaluate(int nfdir, int nadir) {
  // Obtain initial guess
  nlp_solver_.input(NLP_X_INIT).set(output(0));
  
  // Add other arguments
  int k = 0;
  
  for (int i=1;i<f_.getNumInputs();++i) {
    std::copy(input(i-1).data().begin(),input(i-1).data().end(),nlp_solver_.input(NLP_P).data().begin()+k); k+= input(i-1).size();
  }
  
  // Solve NLP
  nlp_solver_.evaluate();

  // Copy the outputs
  output(0).set(nlp_solver_.output(NLP_X_OPT));
  
  // Save auxillary outputs
  for(int i=1; i<getNumOutputs(); ++i){
    output(i).set(f_.output(i));
  }
  
  // End of function if no sensitivities
  if(nfdir==0 && nadir==0)
    return;
  
  // Make sure that a linear solver has been provided
  casadi_assert_message(!linsol_.isNull(),"Sensitivities of an implicit function requires a provided linear solver");
  casadi_assert_message(!J_.isNull(),"Sensitivities of an implicit function requires an exact Jacobian");

  // Pass inputs
  J_.setInput(nlp_solver_.output(NLP_X_OPT),0);
  for(int i=0; i<getNumInputs(); ++i)
    J_.setInput(input(i),i+1);

  // Evaluate jacobian
  J_.evaluate();

  // Pass non-zero elements, scaled by -gamma, to the linear solver
  linsol_.setInput(J_.output(),0);

  // Prepare the solution of the linear system (e.g. factorize)
  linsol_.prepare();
  
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

void NLPImplicitInternal::init(){

  ImplicitFunctionInternal::init();

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
  
  // Generate Jacobian if not provided
  if(J_.isNull()) J_ = f_.jacobian(0,0);
  J_.init();
  
  // Initialize the linear solver, if provided
  if(!linsol_.isNull()){
    linsol_.setSparsity(J_.output().sparsity());
    linsol_.init();
  }
    
  casadi_assert_message(f_.getNumInputs()>0,"NLPImplicitInternal: the supplied f must have at least one input.");
  
  MX V = msym("V",f_.input().sparsity());
  
  std::vector< CRSSparsity > sps;
  for (int k=1;k<f_.getNumInputs();++k)
    sps.push_back(f_.input(k).sparsity());
  
  // So that we can pass it on to createParent
  std::pair< MX, std::vector< MX > > mypair = createParent(sps);

  // V groups all parameters in an MX
  MX P(mypair.first);
  std::vector< MX > inputs(mypair.second);
  
  // We're going to use two-argument objective and constraints to allow the use of parameters
  std::vector< MX > args;
  args.push_back(V);
  args.push_back(P);
    
  MXFunction NLP_f(args,0); NLP_f.init();
  
  std::vector< MX > args_call;
  args_call.push_back(V);
  args_call.insert(args_call.end(),mypair.second.begin(),mypair.second.end());

  MXFunction NLP_g(args,f_.call(args_call)); NLP_g.init();
  
  // Create an nlpsolver instance
  NLPSolverCreator nlp_solvercreator = getOption("nlp_solver");
  nlp_solver_ = nlp_solvercreator(NLP_f,NLP_g,FX(),FX());
  if(hasSetOption("nlp_solver_options")){
    nlp_solver_.setOption(getOption("nlp_solver_options"));
  }
  nlp_solver_.setOption("parametric",true);
  
  // Initialize the NLP solver
  nlp_solver_.init();
  
  
  nlp_solver_.input(NLP_LBG).setAll(0);
  nlp_solver_.input(NLP_UBG).setAll(0);
  
}

} // namespace CasADi

