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
#include "../mx/mx_node.hpp"
#include "../mx/mx_tools.hpp"
#include <iterator>

using namespace std;
namespace CasADi{

  ImplicitFunctionInternal::ImplicitFunctionInternal(const FX& f, const FX& jac, const LinearSolver& linsol) : f_(f), jac_(jac), linsol_(linsol){
    addOption("linear_solver",            OT_LINEARSOLVER, GenericType(), "User-defined linear solver class. Needed for sensitivities.");
    addOption("linear_solver_options",    OT_DICTIONARY,   GenericType(), "Options to be passed to the linear solver.");
  }

  ImplicitFunctionInternal::~ImplicitFunctionInternal(){
  }

  void ImplicitFunctionInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
    FXInternal::deepCopyMembers(already_copied);
    f_ = deepcopy(f_,already_copied);
    jac_ = deepcopy(jac_,already_copied);
    linsol_ = deepcopy(linsol_,already_copied);
  }

  void ImplicitFunctionInternal::init(){
    // Initialize the residual function
    f_.init(false);
  
    // Get the number of equations and check consistency
    casadi_assert_message(f_.output().dense() && f_.output().size2()==1, "Residual must be a dense vector");
    casadi_assert_message(f_.input().dense() && f_.input().size2()==1, "Unknown must be a dense vector");
    n_ = f_.output().size();
    casadi_assert_message(n_ == f_.input().size(), "Dimension mismatch");
    casadi_assert_message(f_.getNumOutputs()==1, "Auxiliary outputs of ImplicitFunctions are no longer allowed, cf. #669");

    // Allocate inputs
    setNumInputs(f_.getNumInputs()-1);
    for(int i=0; i<getNumInputs(); ++i){
      input(i) = f_.input(i+1);
    }
  
    // Allocate output
    setNumOutputs(1);
    output(0) = f_.input(0);
  
    // Call the base class initializer
    FXInternal::init();

    // Generate Jacobian if not provided
    if(jac_.isNull()) jac_ = f_.jacobian(0,0);
    jac_.init(false);
  
    // Check for structural singularity in the Jacobian
    casadi_assert_message(!isSingular(jac_.output().sparsity()),"ImplicitFunctionInternal::init: singularity - the jacobian is structurally rank-deficient. sprank(J)=" << sprank(jac_.output()) << " (instead of "<< jac_.output().size1() << ")");
  
    // Get the linear solver creator function
    if(linsol_.isNull()){
      if(hasSetOption("linear_solver")){
        linearSolverCreator linear_solver_creator = getOption("linear_solver");
        
        // Allocate an NLP solver
        linsol_ = linear_solver_creator(jac_.output().sparsity());
        
        // Pass options
        if(hasSetOption("linear_solver_options")){
          const Dictionary& linear_solver_options = getOption("linear_solver_options");
          linsol_.setOption(linear_solver_options);
        }
        
        // Initialize
        linsol_.init();
      }
    } else {
      // Initialize the linear solver, if provided
      linsol_.init(false);
      casadi_assert(linsol_.input().sparsity()==jac_.output().sparsity());
    }
    
    // Allocate memory for directional derivatives
    ImplicitFunctionInternal::updateNumSens(false);
    
    // No factorization yet;
    fact_up_to_date_ = false;
  }

  void ImplicitFunctionInternal::updateNumSens(bool recursive){
    // Call the base class if needed
    if(recursive) FXInternal::updateNumSens(recursive);
  
    // Request more directional derivatives for the residual function
    f_.requestNumSens(nfdir_,nadir_);
  }

  void ImplicitFunctionInternal::evaluate(int nfdir, int nadir){
    // Mark factorization as out-of-date. TODO: make this conditional
    fact_up_to_date_ = false;

    // Solve the nonlinear system of equations
    solveNonLinear();

    // Quick return if no sensitivities
    if(nfdir==0 && nadir==0) return;

    // Make sure that a linear solver has been provided
    casadi_assert_message(!linsol_.isNull(),"Sensitivities of an implicit function requires a provided linear solver");
    casadi_assert_message(!jac_.isNull(),"Sensitivities of an implicit function requires an exact Jacobian");
  
    // Evaluate and factorize the Jacobian
    if (!fact_up_to_date_) {
      // Pass inputs
      jac_.setInput(output(),0);
      for(int i=0; i<getNumInputs(); ++i)
        jac_.setInput(input(i),i+1);

      // Evaluate jacobian
      jac_.evaluate();

      // Pass non-zero elements, scaled by -gamma, to the linear solver
      linsol_.setInput(jac_.output(),0);

      // Prepare the solution of the linear system (e.g. factorize)
      linsol_.prepare();
      fact_up_to_date_ = true;
    }

  
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
      linsol_.solve(&faseed.front(),1,false);
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
      linsol_.solve(&fsens.front(),1,true);
    }
  
    // Get the adjoint sensitivities
    for(int dir=0; dir<nadir; ++dir){
      for(int i=0; i<getNumInputs(); ++i){
        f_.adjSens(i+1,dir).get(adjSens(i,dir));
      }
    }
  }

  void ImplicitFunctionInternal::evaluateMX(MXNode* node, const MXPtrV& arg, MXPtrV& res, const MXPtrVV& fseed, MXPtrVV& fsens, const MXPtrVV& aseed, MXPtrVV& asens, bool output_given){
    // Evaluate non-differentiated
    vector<MX> argv = MXNode::getVector(arg);
    if(!output_given){
      *res[0] = callSelf(argv).front();
    }

    // Quick return if no derivatives
    int nfwd = fsens.size();
    int nadj = aseed.size();
    if(nfwd==0 && nadj==0) return;

    // Temporaries
    vector<int> row_offset(1,0);
    vector<MX> rhs;

    // Arguments when calling f/f_der
    vector<MX> v;
    int nf_in = f_.getNumInputs();
    v.reserve(nf_in*(1+nfwd) + nadj);
    v.push_back(*res[0]);
    v.insert(v.end(),argv.begin(),argv.end());

    // Get an expression for the Jacobian
    MX J = jac_.call(v).front();

    // Directional derivatives of f
    FX f_der = f_.derivative(nfwd,nadj);

    // Forward sensitivities, collect arguments for calling f_der
    for(int d=0; d<nfwd; ++d){
      v.push_back(MX::sparse(output().shape()));
      argv = MXNode::getVector(fseed[d]);
      v.insert(v.end(),argv.begin(),argv.end());
    }

    // Adjoint sensitivities, solve to get arguments for calling f_der
    if(nadj>0){
      for(int d=0; d<nadj; ++d){
        rhs.push_back(trans(*aseed[d][0]));
        row_offset.push_back(row_offset.back()+1);
        *aseed[d][0] = MX();
      }
      rhs = vertsplit(J->getSolve(vertcat(rhs),false,linsol_),row_offset);
      for(int d=0; d<nadj; ++d){
        v.push_back(trans(rhs[d]));
      }
      row_offset.resize(1);
      rhs.clear();
    }
  
    // Propagate through the implicit function
    v = f_der.call(v);
    vector<MX>::const_iterator v_it = v.begin();

    // Discard non-differentiated evaluation (change?)
    v_it++;

    // Solve for the forward sensitivities
    if(nfwd>0){
      for(int d=0; d<nfwd; ++d){
        rhs.push_back(trans(*v_it++));
        row_offset.push_back(row_offset.back()+1);        
      }
      rhs = vertsplit(J->getSolve(vertcat(rhs),true,linsol_),row_offset);
      for(int d=0; d<nfwd; ++d){
        if(fsens[d][0]!=0){
          *fsens[d][0] = -trans(rhs[d]);
        }
      }
      row_offset.resize(1);
      rhs.clear();
    }

    // Collect adjoint sensitivities
    for(int d=0; d<nadj; ++d){
      // Discard adjoint corresponding to z
      v_it++;

      for(int i=0; i<asens[d].size(); ++i, ++v_it){
        if(asens[d][i]!=0){
          *asens[d][i] = - *v_it;
        }
      }
    }
    casadi_assert(v_it==v.end());
  }

 
} // namespace CasADi

  


