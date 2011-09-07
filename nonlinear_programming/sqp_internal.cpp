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
#include "casadi/matrix/matrix_tools.hpp"
/*#include "interfaces/qpoases/qpoases_solver.hpp"*/
#include <ctime>
#include <iomanip>

using namespace std;
namespace CasADi{

SQPInternal::SQPInternal(const FX& F, const FX& G, const FX& H, const FX& J) : F_(F), G_(G), H_(H), J_(J){
  casadi_warning("The SQP method is under development");
  addOption("qp_solver",         OT_QPSOLVER,   GenericType(), "The QP solver to be used by the SQP method");
  addOption("qp_solver_options", OT_DICTIONARY, GenericType(), "Options to be passed to the QP solver");
  addOption("maxiter",           OT_INTEGER,    100,           "Maximum number of SQP iterations");
  addOption("maxiter_ls",        OT_INTEGER,    100,           "Maximum number of linesearch iterations");
  addOption("toldx",             OT_REAL   ,    1e-12,         "Stopping criterion for the stepsize");
  addOption("tolgl",             OT_REAL   ,    1e-12,         "Stopping criterion for the Lagrangian gradient");
  addOption("sigma",             OT_REAL   ,    1.0,           "Linesearch parameter");
  addOption("rho",               OT_REAL   ,    0.5,           "Linesearch parameter");
  addOption("mu_safety",         OT_REAL   ,    1.1,           "Safety factor for linesearch mu");
  addOption("eta",               OT_REAL   ,    0.0001,        "Linesearch parameter: See Nocedal 3.4");
  addOption("tau",               OT_REAL   ,    0.2,           "Linesearch parameter");
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
  
  // Create Jacobian if necessary
  if(J_.isNull()){
    J_ = G_.jacobian();
  }
  J_.init();
  
  // Read options
  maxiter_ = getOption("maxiter");
  maxiter_ls_ = getOption("maxiter_ls");
  toldx_ = getOption("toldx");
  tolgl_ = getOption("tolgl");
  sigma_ = getOption("sigma");
  rho_ = getOption("rho");
  mu_safety_ = getOption("mu_safety");
  eta_ = getOption("eta");
  tau_ = getOption("tau");
  
  // QP solver
  int n = input(NLP_X_INIT).size();
  
  // Allocate a QP solver
  CRSSparsity H_sparsity = sp_dense(n,n);
  CRSSparsity A_sparsity = J_.output().sparsity();

  QPSolverCreator qp_solver_creator = getOption("qp_solver");
  qp_solver_ = qp_solver_creator(H_sparsity,A_sparsity);
  
  // Set options if provided
  if(hasSetOption("qp_solver_options")){
    Dictionary qp_solver_options = getOption("qp_solver_options");
    qp_solver_.setOption(qp_solver_options);
  }

  qp_solver_.init();
}

void SQPInternal::evaluate(int nfdir, int nadir){
  casadi_assert(nfdir==0 && nadir==0);
  
  // Initial guess
  DMatrix x = input(NLP_X_INIT);

  // Current cost;
  double fk;
  
  // current 'mu' in the T1 merit function
  double merit_mu = 0;  

  // Get dimensions
  int m = G_.output().size(); // Number of equality constraints
  int n = x.size();  // Number of variables

  // Initial guess for the lagrange multipliers
  DMatrix lambda_k(m,1,0);
  DMatrix lambda_x_k(n,1,0);

  // Initial guess for the Hessian
  DMatrix Bk = DMatrix::eye(n);
  makeDense(Bk);

  // No bounds on the control
  double inf = numeric_limits<double>::infinity();
  qp_solver_.input(QP_LBX).setAll(-inf);
  qp_solver_.input(QP_UBX).setAll( inf);

  // Header
  cout << " iter     objective    nls           dx         gradL      eq viol" << endl;
  int k = 0;

  while(true){
    // Evaluate the constraint function
    G_.setInput(x);
    G_.evaluate();
    DMatrix gk = G_.output();
    
    // Evaluate the Jacobian
    J_.setInput(x);
    J_.evaluate();
    DMatrix Jgk = J_.output();
    
    // Evaluate the gradient of the objective function
    F_.setInput(x);
    F_.setAdjSeed(1.0);
    F_.evaluate(0,1);
    fk = F_.output().at(0);
    DMatrix gfk = F_.adjSens();
    
    // Pass data to QP solver
    qp_solver_.setInput(Bk,QP_H);
    qp_solver_.setInput(Jgk,QP_A);
    qp_solver_.setInput(gfk,QP_G);
    qp_solver_.setInput(-gk,QP_LBA);
    qp_solver_.setInput(-gk,QP_UBA);
//     qp_solver_.setInput(input(NLP_LBX),QP_LBX);
//     qp_solver_.setInput(input(NLP_UBX),QP_UBX);

    // Solve the QP subproblem
    qp_solver_.evaluate();

    // Get the optimal solution
    DMatrix p = qp_solver_.output(QP_PRIMAL);
    
    // Get the dual solution for the inequalities
    DMatrix lambda_hat = qp_solver_.output(QP_DUAL_A);
    
    // Get the dual solution for the bounds
    DMatrix lambda_x_hat = qp_solver_.output(QP_DUAL_X);
    
    // Get the gradient of the Lagrangian
    DMatrix gradL = F_.adjSens() - prod(trans(Jgk),lambda_hat) - lambda_x_hat;
    
    // Pass adjoint seeds to g
    //gfcn.setAdjSeed(lambda_hat);
    //gfcn.evaluate(0,1);

    // Do a line search along p
    double mu = merit_mu;
    
    // 1-norm of the feasability violations
    double feasviol = sum(fabs(gk)).at(0);

    // Use a quadratic model of T1 to get a lower bound on mu (eq. 18.36 in Nocedal)
    double mu_lb = ((inner_prod(gfk,p) + sigma_/2.0*prod(trans(p),prod(Bk,p)))/(1.-rho_)/feasviol).at(0);

    // Increase mu if it is below the lower bound
    if(mu < mu_lb){
      mu = mu_lb*mu_safety_;
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
      if(T1_new.at(0) <= (T1 + eta_*alpha*DT1)){
        break;
      }

      // Backtrace
      alpha = alpha*tau_;
      
      // Go to next iteration
      lsiter = lsiter+1;
      if(lsiter >= maxiter_ls_){
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

    cout << setw(5) << k << setw(15) << fk << setw(5) << lsiter << setw(15) << normdx << setw(15) << normgradL << setw(15) << eq_viol << endl;

    // Check convergence on dx
    if(normdx.at(0) < toldx_){
      cout << "Convergence (small dx)" << endl;
      break;
    } else if(normgradL.at(0) < tolgl_){
      cout << "Convergence (small gradL)" << endl;
      break;
    }
      
    // Evaluate the constraint function
    G_.setInput(x);
    G_.evaluate();
    gk = G_.output();
    
    // Evaluate the Jacobian
    J_.setInput(x);
    J_.evaluate();
    Jgk = J_.output();
      
    // Evaluate the gradient of the objective function
    F_.setInput(x);
    F_.setAdjSeed(1.0);
    F_.evaluate(0,1);
    fk = DMatrix(F_.output()).at(0);
    gfk = F_.adjSens();

    // Check if maximum number of iterations reached
    if(k >= maxiter_){
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
  
  output(NLP_COST).set(fk);
  output(NLP_X_OPT).set(x);
}

} // namespace CasADi
