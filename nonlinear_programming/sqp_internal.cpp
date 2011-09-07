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

#include "sqp_internal.hpp"
#include "casadi/stl_vector_tools.hpp"
#include "casadi/matrix/sparsity_tools.hpp"
#include "casadi/fx/qp_solver.hpp"
#include "casadi/matrix/matrix_tools.hpp"
/*#include "interfaces/qpoases/qpoases_solver.hpp"*/
#include <ctime>

using namespace std;
namespace CasADi{

SQPInternal::SQPInternal(const FX& F, const FX& G, const FX& H, const FX& J) : F_(F), G_(G), H_(H), J_(J){
  casadi_warning("The SQP method is under development");
  addOption("qp_solver", OT_QPSOLVER, GenericType(), "The QP solver to be used by the SQP method");
  addOption("qp_solver_options", OT_DICTIONARY, GenericType(), "Options to be passed to the QP solver");
}


SQPInternal::~SQPInternal(){
}

void SQPInternal::init(){
  // Initialize the functions
  F_.init();
  G_.init();
  n_ = F_.input().size();
  m_ = G_.output().size();
  
  // Call the init method of the base class
  NLPSolverInternal::init();
}

void SQPInternal::evaluate(int nfdir, int nadir){
  casadi_assert(nfdir==0 && nadir==0);
  
  
  // Parameters in the algorithm
  int maxiter = 100; // maximum number of sqp iterations
  double toldx = 1e-12; // stopping criterion for the stepsize
  double tolgL = 1e-12; // stopping criterion for the lagrangian gradient
  double merit_mu = 0;  // current 'mu' in the T1 merit function
  
  // Initial guess
  DMatrix x = input(NLP_X_INIT);

  // Get dimensions
  int m = G_.output().size(); // Number of equality constraints
  int n = x.size();  // Number of variables
  int q = 0; // Number of inequality constraints

  // Initial guess for the lagrange multipliers
  DMatrix lambda_k(m,1,0);
  DMatrix lambda_x_k(n,1,0);

  // Initial guess for the Hessian
  DMatrix Bk = DMatrix::eye(n);
  makeDense(Bk);

  // Jacobian function
  FX jfcn = G_.jacobian();
  jfcn.init();

  // Allocate a QP solver
  CRSSparsity H_sparsity = sp_dense(n,n);
  CRSSparsity A_sparsity = jfcn.output().sparsity();

  QPSolverCreator qp_solver_creator = getOption("qp_solver");
  QPSolver qp_solver = qp_solver_creator(H_sparsity,A_sparsity);
  
  // Set options if provided
  if(hasSetOption("qp_solver_options")){
    Dictionary qp_solver_options = getOption("qp_solver_options");
    qp_solver.setOption(qp_solver_options);
  }

  qp_solver.init();

  // No bounds on the control
  double inf = numeric_limits<double>::infinity();
  qp_solver.input(QP_LBX).setAll(-inf);
  qp_solver.input(QP_UBX).setAll( inf);

  // Header
  cout << " k  nls | dx         gradL      eq viol    ineq viol" << endl;
  int k = 0;

  while(true){
    // Evaluate the constraint function
    G_.setInput(x);
    G_.evaluate();
    DMatrix gk = G_.output();
    
    // Evaluate the Jacobian
    jfcn.setInput(x);
    jfcn.evaluate();
    DMatrix Jgk = jfcn.output();
    
    // Evaluate the gradient of the objective function
    F_.setInput(x);
    F_.setAdjSeed(1.0);
    F_.evaluate(0,1);
    double fk = F_.output().at(0);
    DMatrix gfk = F_.adjSens();
    
    // Pass data to QP solver
    qp_solver.setInput(Bk,QP_H);
    qp_solver.setInput(Jgk,QP_A);
    qp_solver.setInput(gfk,QP_G);
    qp_solver.setInput(-gk,QP_LBA);
    qp_solver.setInput(-gk,QP_UBA);

    // Solve the QP subproblem
    qp_solver.evaluate();

    // Get the optimal solution
    DMatrix p = qp_solver.output(QP_PRIMAL);
    
    // Get the dual solution for the inequalities
    DMatrix lambda_hat = qp_solver.output(QP_DUAL_A);
    
    // Get the dual solution for the bounds
    DMatrix lambda_x_hat = qp_solver.output(QP_DUAL_X);
    
    // Get the gradient of the Lagrangian
    DMatrix gradL = F_.adjSens() - prod(trans(Jgk),lambda_hat) - lambda_x_hat;
    
    // Pass adjoint seeds to g
    //gfcn.setAdjSeed(lambda_hat);
    //gfcn.evaluate(0,1);

    // Do a line search along p
    double mu = merit_mu;
    
    // parameters in the algorithm
    double sigma = 1.;  // Bk in BDGS is always pos.def.
    double rho = 0.5;
    double mu_safety = 1.1; // safety factor for mu (see below)
    double eta = 0.0001; // text to Noc 3.4
    double tau = 0.2;
    int maxiter = 100;

    // 1-norm of the feasability violations
    double feasviol = sum(fabs(gk)).at(0);

    // Use a quadratic model of T1 to get a lower bound on mu (eq. 18.36 in Nocedal)
    double mu_lb = ((inner_prod(gfk,p) + sigma/2.0*dot(trans(p),dot(Bk,p)))/(1.-rho)/feasviol).at(0);

    // Increase mu if it is below the lower bound
    if(mu < mu_lb){
      mu = mu_lb*mu_safety;
    }

    // Calculate T1 at x (18.27 in Nocedal)
    double T1 = fk + mu*feasviol;

    // Calculate the directional derivative of T1 at x (cf. 18.29 in Nocedal)
    double DT1 = (inner_prod(gfk,p) - mu*sum(fabs(gk))).at(0);
    
    int lsiter = 0;
    double alpha = 1;
    while(true){
      // Evaluate prospective x
      DMatrix x_new = x+alpha*p;
      F_.setInput(x_new);
      F_.evaluate();
      DMatrix fk_new = DMatrix(F_.output());

      // Evaluate gk, hk and get 1-norm of the feasability violations
      G_.setInput(x_new);
      G_.evaluate();
      DMatrix gk_new = G_.output();
      DMatrix feasviol_new = sum(fabs(gk_new));

      // New T1 function
      DMatrix T1_new = fk_new + mu*feasviol_new;

      // Check Armijo condition, SQP version (18.28 in Nocedal)
      if(T1_new.at(0) <= (T1 + eta*alpha*DT1)){
        break;
      }

      // Backtrace
      alpha = alpha*tau;
      
      // Go to next iteration
      lsiter = lsiter+1;
      if(lsiter >= maxiter){
        throw CasadiException("linesearch failed!");
      }
    }

    // Step size
    double tk = alpha;

    // Calculate the new step
    DMatrix dx = p*tk;
    x = x + dx;
    lambda_k = tk*lambda_hat + (1-tk)*lambda_k;
    lambda_x_k = tk*lambda_x_hat + (1-tk)*lambda_x_k;
    k = k+1;

    // Gather and print iteration information
    DMatrix normdx = norm_2(dx); // step size
    DMatrix normgradL = norm_2(gradL); // size of the Lagrangian gradient
    DMatrix eq_viol = sum(fabs(gk)); // constraint violation
    string ineq_viol = "nan"; // sum(max(0,-hk)); % inequality constraint violation

    cout << k << " " << lsiter << " " << normdx << " " << normgradL << " " << eq_viol << " " << ineq_viol << endl;

    // Check convergence on dx
    if(normdx.at(0) < toldx){
      cout << "Convergence (small dx)" << endl;
      break;
    } else if(normgradL.at(0) < tolgL){
      cout << "Convergence (small gradL)" << endl;
      break;
    }
      
    // Evaluate the constraint function
    G_.setInput(x);
    G_.evaluate();
    gk = G_.output();
    
    // Evaluate the Jacobian
    jfcn.setInput(x);
    jfcn.evaluate();
    Jgk = jfcn.output();
      
    // Evaluate the gradient of the objective function
    F_.setInput(x);
    F_.setAdjSeed(1.0);
    F_.evaluate(0,1);
    fk = DMatrix(F_.output()).at(0);
    gfk = F_.adjSens();

    // Check if maximum number of iterations reached
    if(k >= maxiter){
      cout << "Maximum number of SQP iterations reached!" << endl;
      break;
    }

    // Complete the damped BFGS update (Procedure 18.2 in Nocedal)
    DMatrix gradL_new = gfk - prod(trans(Jgk),lambda_k) - lambda_x_k;
    DMatrix yk = gradL_new - gradL;
    DMatrix Bdx = prod(Bk,dx);
    DMatrix dxBdx = prod(trans(dx),Bdx);
    DMatrix ydx = inner_prod(dx,yk);
    DMatrix thetak;
    if(ydx.at(0) >= 0.2*dxBdx.at(0)){
      thetak = 1.;
    } else {
      thetak = 0.8*dxBdx/(dxBdx - ydx);
    }
    DMatrix rk = thetak*dx + (1-thetak)*Bdx; // rk replaces yk to assure Bk pos.def.
    Bk = Bk - outer_prod(Bdx,Bdx)/dxBdx + outer_prod(rk,rk)/ inner_prod(rk,dx);
  }
  cout << "SQP algorithm terminated after " << (k-1) << " iterations" << endl;
  
  
}

} // namespace CasADi
