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

#include "ip_internal.hpp"
#include "symbolic/stl_vector_tools.hpp"
#include "symbolic/matrix/sparsity_tools.hpp"
#include "symbolic/matrix/matrix_tools.hpp"
#include "symbolic/fx/sx_function.hpp"
#include "symbolic/sx/sx_tools.hpp"
#include "symbolic/casadi_calculus.hpp"
#include <ctime>
#include <iomanip>

using namespace std;
namespace CasADi{

IPInternal::IPInternal(const FX& F, const FX& G) : NLPSolverInternal(FX(),F,G){
  casadi_warning("The IP method is experimental and incomplete. Can be used as the basis of an IP solver in CasADi.");
  addOption("linear_solver",         OT_LINEARSOLVER,   GenericType(), "The linear solver to be used by the IP method");
  addOption("linear_solver_options", OT_DICTIONARY, GenericType(), "Options to be passed to the linear solver");
}

IPInternal::~IPInternal(){
}

void IPInternal::init(){
  // Call the init method of the base class
  NLPSolverInternal::init();
    
  // Assume SXFunction
  SXFunction FF = shared_cast<SXFunction>(F_);
  SXFunction GG = shared_cast<SXFunction>(G_);

  // Split up the problem
  SXMatrix x = FF.inputExpr(0);
  SXMatrix f = FF.outputExpr(0);
  SXMatrix g = GG.outputExpr(0);
//   cout << "x = " << x << endl;
//   cout << "f = " << f << endl;
//   cout << "g = " << g << endl;
  
  // Barrier parameter
  SXMatrix t = ssym("t");
  
  // Objective of the equality constraint problem
  SXMatrix f_eq = t*f - sumAll(CasADi::log(x));
//   cout << "f_eq = " << f_eq << endl;
  
  // Hessian of the objective
  SXMatrix H = CasADi::hessian(f_eq,x);
//  cout << "H = " << H << endl;
  
  // Jacobian of the constraints
  SXMatrix A = CasADi::jacobian(g,x);
//  cout << "A = " << A << endl;
  
  // Form the KKT matrix
  SXMatrix K = vertcat(horzcat(H,trans(A)),horzcat(A,SXMatrix::sparse(ng_,ng_)));
  if(verbose()){
    cout << "K = " << K << endl;
  }
  
  // Form the right hand side of the KKT system
  SXMatrix k = vertcat(-CasADi::gradient(f_eq,x),SXMatrix::sparse(ng_));
  makeDense(k);
  if(verbose()){
    cout << "k = " << k << endl;
  }
  
  // Create a function that forms the KKT system
  SXMatrix kfcn_in[] = {x,t};
  SXMatrix kfcn_out[] = {K,k};
  kfcn_ = SXFunction(vector<SXMatrix>(kfcn_in,kfcn_in+2),vector<SXMatrix>(kfcn_out,kfcn_out+2));
  kfcn_.init();
  
  // Create a linear solver for the KKT system
  linearSolverCreator linear_solver_creator = getOption("linear_solver");
  linear_solver_ = linear_solver_creator(K.sparsity());
  linear_solver_.init();  
}

void IPInternal::evaluate(int nfdir, int nadir){
  casadi_assert(nfdir==0 && nadir==0);

  // Current value of x
  DMatrix &x_curr = output(NLP_SOLVER_X);
  copy(input(NLP_SOLVER_X0).begin(),input(NLP_SOLVER_X0).end(),x_curr.begin());
  
  // Step size
  DMatrix dx = DMatrix(x_curr.sparsity(),0);
  
  // Barrier parameter
  double t = 1;
  
  // Backtracking parameter for the barrier parameter
  double t_backtrack = 0.5;
  
  // Maximum number of newton iterations
  int max_newton_iter = 10;
  
  // Maximum number of backtrackings
  int max_backtrack = 15;
  
  // Backtracking iterations
  for(int backtrack=0; backtrack<max_backtrack; ++backtrack){
  
    // Newton iterations
    for(int newton_iter=0; newton_iter<max_newton_iter; ++newton_iter){
    
      // Form the kkt system
      kfcn_.setInput(x_curr,K_x);
      kfcn_.setInput(t,K_t);
      kfcn_.evaluate();
      
      // Solve the linear solver
      linear_solver_.setInput(kfcn_.output(K_K),0);
      linear_solver_.setInput(kfcn_.output(K_k),1);
      linear_solver_.evaluate();
      const DMatrix& dx_nu = linear_solver_.output();
      copy(dx_nu.begin(),dx_nu.begin()+dx.size(),dx.begin());
      
      // Update x
      std::transform(dx.begin(),dx.end(),x_curr.begin(),x_curr.begin(),std::plus<double>());
    }
    
    // Update backtracking parameter
    t *= t_backtrack;
  }
}

} // namespace CasADi
