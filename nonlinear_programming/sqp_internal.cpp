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
#include "symbolic/stl_vector_tools.hpp"
#include "symbolic/matrix/sparsity_tools.hpp"
#include "symbolic/matrix/matrix_tools.hpp"
#include "symbolic/fx/sx_function.hpp"
#include "symbolic/sx/sx_tools.hpp"
#include "symbolic/casadi_calculus.hpp"
#include <ctime>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cfloat>

using namespace std;
namespace CasADi{

double inf = numeric_limits<double>::infinity();

SQPInternal::SQPInternal(const FX& F, const FX& G, const FX& H, const FX& J) : NLPSolverInternal(F,G,H,J){
  casadi_warning("The SQP method is under development");
  addOption("qp_solver",         OT_QPSOLVER,   GenericType(),    "The QP solver to be used by the SQP method");
  addOption("qp_solver_options", OT_DICTIONARY, GenericType(),    "Options to be passed to the QP solver");
  addOption("hessian_approximation", OT_STRING, "limited-memory", "limited-memory|exact");
  addOption("maxiter",           OT_INTEGER,      50,             "Maximum number of SQP iterations");
  addOption("maxiter_ls",        OT_INTEGER,       3,             "Maximum number of linesearch iterations");
  addOption("tol_pr",            OT_REAL,       1e-6,             "Stopping criterion for primal infeasibility");
  addOption("tol_du",            OT_REAL,       1e-6,             "Stopping criterion for dual infeasability");
  addOption("c1",                OT_REAL,       1E-4,             "Armijo condition, coefficient of decrease in merit");
  addOption("beta",              OT_REAL,       0.8,              "Line-search parameter, restoration factor of stepsize");
  addOption("merit_memory",      OT_INTEGER,      4,              "Size of memory to store history of merit function values");
  addOption("lbfgs_memory",      OT_INTEGER,     10,              "Size of L-BFGS memory.");
  addOption("regularize",        OT_BOOLEAN,  false,              "Automatic regularization of Lagrange Hessian.");
  
  // Monitors
  addOption("monitor",      OT_STRINGVECTOR, GenericType(),  "", "eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h|qp|dx", true);
}


SQPInternal::~SQPInternal(){
}

void SQPInternal::init(){
  // Call the init method of the base class
  NLPSolverInternal::init();
    
  // Read options
  maxiter_ = getOption("maxiter");
  maxiter_ls_ = getOption("maxiter_ls");
  c1_ = getOption("c1");
  beta_ = getOption("beta");
  merit_memsize_ = getOption("merit_memory");
  lbfgs_memory_ = getOption("lbfgs_memory");
  tol_pr_ = getOption("tol_pr");
  tol_du_ = getOption("tol_du");
  regularize_ = getOption("regularize");
  if(getOption("hessian_approximation")=="exact")
    hess_mode_ = HESS_EXACT;
  else if(getOption("hessian_approximation")=="limited-memory")
    hess_mode_ = HESS_BFGS;
   
  if (hess_mode_== HESS_EXACT && H_.isNull()) {
    if (!getOption("generate_hessian")){
      casadi_error("SQPInternal::evaluate: you set option 'hessian_approximation' to 'exact', but no hessian was supplied. Try with option \"generate_hessian\".");
    }
  }
  
  // If the Hessian is generated, we use exact approximation by default
  if (bool(getOption("generate_hessian"))){
    setOption("hessian_approximation", "exact");
  }
  
  // Allocate a QP solver
  CRSSparsity H_sparsity = hess_mode_==HESS_EXACT ? H_.output().sparsity() : sp_dense(n_,n_);
  H_sparsity = H_sparsity + DMatrix::eye(n_).sparsity();
  CRSSparsity A_sparsity = J_.isNull() ? CRSSparsity(0,n_,false) : J_.output().sparsity();

  QPSolverCreator qp_solver_creator = getOption("qp_solver");
  qp_solver_ = qp_solver_creator(H_sparsity,A_sparsity);

  // Set options if provided
  if(hasSetOption("qp_solver_options")){
    Dictionary qp_solver_options = getOption("qp_solver_options");
    qp_solver_.setOption(qp_solver_options);
  }
  qp_solver_.init();
  
  // Lagrange multipliers of the NLP
  mu_.resize(m_);
  mu_x_.resize(n_);
  
  // Lagrange gradient in the next iterate
  gLag_.resize(n_);
  gLag_old_.resize(n_);

  // Current linearization point
  x_.resize(n_);
  x_cand_.resize(n_);
  x_old_.resize(n_);

  // Constraint function value
  gk_.resize(m_);
  gk_cand_.resize(m_);
  
  // Hessian approximation
  Bk_ = DMatrix(H_sparsity);
  
  // Gradient of the objective
  gf_.resize(n_);

  // Create Hessian update function
  if(hess_mode_ == HESS_BFGS){
    // Create expressions corresponding to Bk, x, x_old, gLag and gLag_old
    SXMatrix Bk = ssym("Bk",H_sparsity);
    SXMatrix x = ssym("x",input(NLP_X_INIT).sparsity());
    SXMatrix x_old = ssym("x",x.sparsity());
    SXMatrix gLag = ssym("gLag",x.sparsity());
    SXMatrix gLag_old = ssym("gLag_old",x.sparsity());
    
    SXMatrix sk = x - x_old;
    SXMatrix yk = gLag - gLag_old;
    SXMatrix qk = mul(Bk, sk);
    
    // Calculating theta
    SXMatrix skBksk = inner_prod(sk, qk);
    SXMatrix omega = if_else(inner_prod(yk, sk) < 0.2 * inner_prod(sk, qk),
                             0.8 * skBksk / (skBksk - inner_prod(sk, yk)),
                             1);
    yk = omega * yk + (1 - omega) * qk;
    SXMatrix theta = 1. / inner_prod(sk, yk);
    SXMatrix phi = 1. / inner_prod(qk, sk);
    SXMatrix Bk_new = Bk + theta * mul(yk, trans(yk)) - phi * mul(qk, trans(qk));
    
    // Inputs of the BFGS update function
    vector<SXMatrix> bfgs_in(BFGS_NUM_IN);
    bfgs_in[BFGS_BK] = Bk;
    bfgs_in[BFGS_X] = x;
    bfgs_in[BFGS_X_OLD] = x_old;
    bfgs_in[BFGS_GLAG] = gLag;
    bfgs_in[BFGS_GLAG_OLD] = gLag_old;
    bfgs_ = SXFunction(bfgs_in,Bk_new);
    bfgs_.setOption("number_of_fwd_dir",0);
    bfgs_.setOption("number_of_adj_dir",0);
    bfgs_.init();
  }
}

void SQPInternal::evaluate(int nfdir, int nadir){
  casadi_assert(nfdir==0 && nadir==0);
  
  checkInitialBounds();
  
  // Get problem data
  const vector<double>& x_init = input(NLP_X_INIT).data();
  const vector<double>& lbx = input(NLP_LBX).data();
  const vector<double>& ubx = input(NLP_UBX).data();
  const vector<double>& lbg = input(NLP_LBG).data();
  const vector<double>& ubg = input(NLP_UBG).data();
  
  // Set the static parameter
  if (parametric_) {
    const vector<double>& p = input(NLP_P).data();
    if (!F_.isNull()) F_.setInput(p,F_.getNumInputs()-1);
    if (!G_.isNull()) G_.setInput(p,G_.getNumInputs()-1);
    if (!H_.isNull()) H_.setInput(p,H_.getNumInputs()-1);
    if (!J_.isNull()) J_.setInput(p,J_.getNumInputs()-1);
  }
    
  // Set linearization point to initial guess
  copy(x_init.begin(),x_init.end(),x_.begin());
  
  // Lagrange multipliers of the NLP
  fill(mu_.begin(),mu_.end(),0);
  fill(mu_x_.begin(),mu_x_.end(),0);

  // Lagrange gradient in the next iterate
  fill(gLag_.begin(),gLag_.end(),0);

  // Reset the Hessian or Hessian approximation
  reset_h();

  // Reset
  merit_mem_.clear();
  sigma_ = 0.;

  // Print header
  printIteration(cout);
  
  // MAIN OPTIMIZATION LOOP
  int iter = 1;
  while(true){
    // Print header occasionally
    if(iter % 10 == 0) printIteration(cout);
    
    // Evaluating Hessian if needed
    eval_h(x_,mu_,1.0,Bk_);

    if(m_>0){
      // Evaluate the constraint function (NOTE: This is not needed. The constraint function should be evaluated as a byproduct of the Jacobian below)
      eval_g(x_,gk_);
      
      // Evaluate the constraint Jacobian
      eval_jac_g();
    }
    
    // Evaluate the gradient of the objective function (NOTE: This should not be needed when exact Hessian is used. The Hessian of the Lagrangian gives the gradient of the Lagrangian. The gradient of the objective function can be obtained from knowing the Jacobian of the constraints. The calculation could be desirable anyway, due to numeric cancellation
    eval_grad_f();
    const DMatrix& gfk = F_.adjSens();

    // Pass data to QP solver
    qp_solver_.setInput(Bk_, QP_H);
    qp_solver_.setInput(gfk,QP_G);
    // Hot-starting if possible
    qp_solver_.setInput(qp_solver_.output(QP_PRIMAL), QP_X_INIT);
    //TODO: Fix hot-starting of dual variables
    //qp_solver_.setInput(mu_qp, QP_LAMBDA_INIT);
      
    if(m_>0){
      qp_solver_.setInput(J_.output(),QP_A);
      transform(lbg.begin(),lbg.end(),gk_.begin(),qp_solver_.input(QP_LBA).begin(),minus<double>());
      transform(ubg.begin(),ubg.end(),gk_.begin(),qp_solver_.input(QP_UBA).begin(),minus<double>());
    }

    transform(lbx.begin(),lbx.end(),x_.begin(),qp_solver_.input(QP_LBX).begin(),minus<double>());
    transform(ubx.begin(),ubx.end(),x_.begin(),qp_solver_.input(QP_UBX).begin(),minus<double>());
    
    if (monitored("qp")) {
      cout << "(main loop) QP_H = " << endl;
      qp_solver_.input(QP_H).printDense();
      cout << "(main loop) QP_A = " << endl;
      qp_solver_.input(QP_A).printDense();
      cout << "(main loop) QP_G = " << endl;
      qp_solver_.input(QP_G).printDense();
      cout << "(main loop) QP_LBA = " << endl;
      qp_solver_.input(QP_LBA).printDense();
      cout << "(main loop) QP_UBA = " << endl;
      qp_solver_.input(QP_UBA).printDense();
      cout << "(main loop) QP_LBX = " << endl;
      qp_solver_.input(QP_LBX).printDense();
      cout << "(main loop) QP_UBX = " << endl;
      qp_solver_.input(QP_UBX).printDense();
    }

    // Solve the QP subproblem
    qp_solver_.evaluate();

    // Get the optimal solution
    const vector<double>& dx = qp_solver_.output(QP_PRIMAL).data();
    if (monitored("dx")){
      cout << "(main loop) dx = " << endl;
      cout << dx << endl;
    }
    // Detecting indefiniteness
//    if ((norm_2(p) / norm_2(x)).at(0) > 500.){
//      casadi_warning("Search direction has very large values, indefinite Hessian might have ouccured.");
//    }
    double gain = quad_form(dx,Bk_);
    if (gain < 0){
      casadi_warning("Indefinite Hessian detected...");
    }
        
    // Get the dual solution for the inequalities
    const vector<double>& mu_qp = qp_solver_.output(QP_LAMBDA_A).data();
    const vector<double>& mu_x_qp = qp_solver_.output(QP_LAMBDA_X).data();

    // Calculate penalty parameter of merit function
    for(int j=0; j<m_; ++j){
      if( fabs(mu_qp[j]) > sigma_){
        sigma_ = fabs(mu_qp[j]) * 1.01;
      }
    }
//    for(int j = 0; j < n; ++j){
//      if( fabs(mu_x_qp[j]) > sigma_){
//        sigma_ = fabs(mu_x_qp[j]) * 1.01;
//      }
//    }

    // Calculate L1-merit function in the actual iterate
    double l1_infeas = 0.;
    for(int j=0; j<m_; ++j){
      // Left-hand side violated
      if(lbg[j] - gk_[j] > 0.){
        l1_infeas += lbg[j] - gk_[j];
      }
      else if (gk_[j] - ubg[j] > 0.){
        l1_infeas += gk_[j] - ubg[j];
      }
    }

    // Right-hand side of Armijo condition
    F_.setFwdSeed(dx);
    F_.evaluate(1, 0);
    double F_sens;
    F_.getFwdSens(F_sens);
    
    double L1dir = F_sens - sigma_ * l1_infeas;
    double L1merit = fk_ + sigma_ * l1_infeas;

    // Storing the actual merit function value in a list
    merit_mem_.push_back(L1merit);
    if (merit_mem_.size() > merit_memsize_){
      merit_mem_.pop_front();
    }

    // Default stepsize
    double t = 1.0;   
    double fk_cand;
    // Merit function value in candidate
    double L1merit_cand = 0.;

    // Line-search loop
    int ls_counter = 1;
    while (true){
      for(int i=0; i<n_; ++i) x_cand_[i] = x_[i] + t * dx[i]; 
      // Evaluating objective and constraints
      F_.setInput(x_cand_);
      F_.evaluate();
      F_.getOutput(fk_cand);
      l1_infeas = 0.;
      if (!G_.isNull()){
        G_.setInput(x_cand_);
        G_.evaluate();
        G_.getOutput(gk_cand_);

        // Calculating merit-function in candidate
        for(int j=0; j<m_; ++j){
          // Left-hand side violated
          if (lbg[j] - gk_cand_[j] > 0.){
            l1_infeas += lbg[j] - gk_cand_[j];
          }
          else if (gk_cand_[j] - ubg[j] > 0.){
            l1_infeas += gk_cand_[j] - ubg[j];
          }
        }
      }
      L1merit_cand = fk_cand + sigma_ * l1_infeas;
      // Calculating maximal merit function value so far
      double meritmax = -1E20;
      for(int k = 0; k < merit_mem_.size(); ++k){
        if (merit_mem_[k] > meritmax){
          meritmax = merit_mem_[k];
        }
      }
      if (L1merit_cand <= meritmax + t * c1_ * L1dir){ 
        // Accepting candidate
        break;
      }
      else{
        // Backtracking
        t = beta_ * t; 
      }

      // Line-search not successful, but we accept it.
      if(ls_counter == maxiter_ls_){
        break;
      }
      ++ls_counter;
    }
    // Candidate accepted
    copy(x_.begin(),x_.end(),x_old_.begin());
    copy(x_cand_.begin(),x_cand_.end(),x_.begin());
    fk_ = fk_cand;
    copy(gk_cand_.begin(),gk_cand_.end(),gk_.begin());
    for(int i=0; i<m_; ++i) mu_[i] = t * mu_qp[i] + (1 - t) * mu_[i];
    for(int i=0; i<n_; ++i) mu_x_[i] = t * mu_x_qp[i] + (1 - t) * mu_x_[i];

    // Evaluating objective gradient
    F_.setInput(x_);
    F_.setAdjSeed(1.0);
    F_.evaluate(0, 1);
    F_.getAdjSens(gLag_); 

    // Adjoint derivative of constraint function
    if(m_>0){
      G_.setAdjSeed(mu_);
      G_.evaluate(0, 1);
      transform(gLag_.begin(),gLag_.end(),G_.adjSens().begin(),gLag_.begin(),plus<double>()); // gLag_ += G_.adjSens()
    }
    transform(gLag_.begin(),gLag_.end(),mu_x_.begin(),gLag_.begin(),plus<double>()); // gLag_ += mu_x_;

    F_.setInput(x_old_);
    F_.setAdjSeed(1.0);
    F_.evaluate(0, 1);
    F_.getAdjSens(gLag_old_);
    if(m_>0){
      G_.setInput(x_old_);
      G_.setAdjSeed(mu_);
      G_.evaluate(0, 1);
      transform(gLag_old_.begin(),gLag_old_.end(),G_.adjSens().begin(),gLag_old_.begin(),plus<double>()); // gLag_old_ += G_.adjSens();
    }
    transform(gLag_old_.begin(),gLag_old_.end(),mu_x_.begin(),gLag_old_.begin(),plus<double>()); // gLag_old += mu_x_;

    // Updating Lagrange Hessian if needed. (BFGS with careful updates and restarts)
    if( hess_mode_ == HESS_BFGS){
      if (iter % lbfgs_memory_ == 0){
        // Remove off-diagonal entries
        const vector<int>& rowind = Bk_.rowind();
        const vector<int>& col = Bk_.col();
        vector<double>& data = Bk_.data();
        for(int i=0; i<rowind.size()-1; ++i){
          for(int el=rowind[i]; el<rowind[i+1]; ++el){
            int j = col[el];
            if(i!=j){
              data[el] = 0;
            }
          }
        }
      }
      
      // Pass to BFGS update function
      bfgs_.setInput(Bk_,BFGS_BK);
      bfgs_.setInput(x_,BFGS_X);
      bfgs_.setInput(x_old_,BFGS_X_OLD);
      bfgs_.setInput(gLag_,BFGS_GLAG);
      bfgs_.setInput(gLag_old_,BFGS_GLAG_OLD);
      
      // Update the Hessian approximation
      bfgs_.evaluate();
      
      // Get the updated Hessian
      bfgs_.getOutput(Bk_);
    }
    // Calculating optimality criterion
    // Primal infeasability
    double pr_infG = 0.;
    double pr_infx = 0.;
    if (!G_.isNull()){
      // Nonlinear constraints
      for(int j=0; j<m_; ++j){
        // Equality
        if (ubg[j] - lbg[j] < 1E-20){
          pr_infG += fabs(gk_cand_[j] - lbg[j]);
        }
        // Inequality, left-hand side violated
        else if(lbg[j] - gk_cand_[j] > 0.){
          pr_infG += lbg[j] - gk_cand_[j];
        }
        // Inequality, right-hand side violated
        else if(gk_cand_[j] - ubg[j] > 0.){
          pr_infG += gk_cand_[j] - ubg[j];
          //cout << color << "SQP: " << mu.elem(j) << defcol << endl;
        }
      }
      // Bound constraints
      for(int j=0; j<n_; ++j){
        // Equality
        if (ubx[j] - lbx[j] < 1E-20){
          pr_infx += fabs(x_[j] - lbx[j]);
        }
        // Inequality, left-hand side violated
        else if ( lbx[j] - x_[j] > 0.){
          pr_infx += lbx[j] - x_[j];
        }
        // Inequality, right-hand side violated
        else if ( x_[j] - ubx[j] > 0.){
          pr_infx += x_[j] - ubx[j];
        }
      }
    }
    double pr_inf = pr_infG + pr_infx;
    
    // 1-norm of lagrange gradient
    double gLag_norm1 = 0;
    for(vector<double>::const_iterator it=gLag_.begin(); it!=gLag_.end(); ++it) gLag_norm1 += fabs(*it);

    // 1-norm of step
    double dx_norm1 = 0;
    for(vector<double>::const_iterator it=dx.begin(); it!=dx.end(); ++it) dx_norm1 += fabs(*it);
    
    // Printing information about the actual iterate
    printIteration(cout,iter,fk_cand,pr_inf,gLag_norm1,dx_norm1,t,ls_counter!=maxiter_ls_,ls_counter);
    
    // Call callback function if present
    if (!callback_.isNull()) {
      callback_.input(NLP_COST).set(fk_);
      callback_.input(NLP_X_OPT).set(x_);
      callback_.input(NLP_LAMBDA_G).set(mu_);
      callback_.input(NLP_LAMBDA_X).set(mu_x_);
      callback_.input(NLP_G).set(gk_);
      callback_.evaluate();
      
      if (callback_.output(0).at(0)) {
       cout << "SQP: aborted by callback...\n"; 
       break;
      }
    }

    // Checking convergence criteria
    if (pr_inf < tol_pr_ && gLag_norm1 < tol_du_){
      cout << "SQP: Convergence achieved after " << iter << " iterations.\n";
      break;
    }

    if (iter == maxiter_){
      cout << "SQP: Maximum number of iterations reached, quiting...\n"; 
      break;
    }
    ++iter;
  }
  
  // Save results to outputs
  output(NLP_COST).set(fk_);
  output(NLP_X_OPT).set(x_);
  output(NLP_LAMBDA_G).set(mu_);
  output(NLP_LAMBDA_X).set(mu_x_);
  output(NLP_G).set(gk_);
  
  // Save statistics
  stats_["iter_count"] = iter;
}

void SQPInternal::printIteration(std::ostream &stream){
  const int w=15;
  stream << setw(w) << "iter";
  stream << setw(w) << "obj";
  stream << setw(w) << "pr_inf";
  stream << setw(w) << "du_inf";
  stream << setw(w) << "corr_norm";
  stream << setw(w) << "ls_param";
  stream << ' ';
  stream << setw(w) << "ls_trials";
  stream << endl;
}
  
void SQPInternal::printIteration(std::ostream &stream, int iter, double obj, double pr_inf, double du_inf, 
                                 double corr_norm, double ls_param, bool ls_success, int ls_trials){
  const int w=15;
  stream << scientific;
  stream << setw(w) << iter;
  stream << setw(w) << obj;
  stream << setw(w) << pr_inf;
  stream << setw(w) << du_inf;
  stream << setw(w) << corr_norm;
  stream << setw(w) << ls_param;
  stream << (ls_success ? ' ' : 'F');
  stream << setw(w) << ls_trials;
  stream << endl;
}

double SQPInternal::quad_form(const std::vector<double>& x, const DMatrix& A){
  // Assert dimensions
  casadi_assert(x.size()==A.size1() && x.size()==A.size2());
  
  // Access the internal data of A
  const std::vector<int> &A_rowind = A.rowind();
  const std::vector<int> &A_col = A.col();
  const std::vector<double> &A_data = A.data();
  
  // Return value
  double ret=0;

  // Loop over the rows of A
  for(int i=0; i<x.size(); ++i){
    // Loop over the nonzeros of A
    for(int el=A_rowind[i]; el<A_rowind[i+1]; ++el){
      // Get column
      int j = A_col[el];
      
      // Add contribution
      ret += x[i]*A_data[el]*x[j];
    }
  }
  
  return ret;
}

void SQPInternal::reset_h(){
  // Initial Hessian approximation of BFGS
  if ( hess_mode_ == HESS_BFGS){
    Bk_.set(DMatrix::eye(n_));
  }

  if (monitored("eval_h")) {
    cout << "(pre) B = " << endl;
    Bk_.printSparse();
  }
}

  void SQPInternal::eval_h(const std::vector<double>& x, const std::vector<double>& lambda, double sigma, Matrix<double>& H){

  if(hess_mode_==HESS_EXACT){

    int n_hess_in = H_.getNumInputs() - (parametric_ ? 1 : 0);
    H_.setInput(x);
    if(n_hess_in>1){
      H_.setInput(lambda, n_hess_in == 4 ? 2 : 1);
      H_.setInput(sigma, n_hess_in == 4 ? 3 : 2);
    }
    H_.evaluate();
    H_.getOutput(H);
    // Determing regularization parameter with Gershgorin theorem
    if(regularize_){
      const vector<int>& rowind = H.rowind();
      const vector<int>& col = H.col();
      vector<double>& data = H.data();
      double reg_param = 0;
      for(int i=0; i<rowind.size()-1; ++i){
	double mineig = 0;
	for(int el=rowind[i]; el<rowind[i+1]; ++el){
	  int j = col[el];
	  if(i == j){
	    mineig += data[el];
	  } else {
	    mineig -= fabs(data[el]);
	  }
	  //          cout << "(" << r << "," << col[el] << "): " << data[el] << endl; 
	}
	reg_param = fmin(reg_param,mineig);
      }
      //      cout << "Regularization parameter: " << -reg_param << endl;
      if ( reg_param < 0.){
	for(int i=0; i<rowind.size()-1; ++i){
	  for(int el=rowind[i]; el<rowind[i+1]; ++el){
	    int j = col[el];
	    if(i==j){
	      data[el] += -reg_param;
	    }
	  }
	}
      }
    }
  }
  
  if (monitored("eval_h")) {
    cout << "(main loop) B = " << endl;
    H.printSparse();
  }
}

  void SQPInternal::eval_g(const std::vector<double>& x, std::vector<double>& g){
  G_.setInput(x);
  G_.evaluate();
  G_.output().get(g,DENSE);
  
  if (monitored("eval_g")) {
    cout << "(main loop) x = " << x << endl;
    cout << "(main loop) G = " << g << endl;
  } 
}

void SQPInternal::eval_jac_g(){
  J_.setInput(x_);
  J_.evaluate();
  
  if (monitored("eval_jac_g")) {
    cout << "(main loop) x = " << x_ << endl;
    cout << "(main loop) J = " << endl;
    J_.output().printSparse();
  }
}

void SQPInternal::eval_grad_f(){
  F_.setInput(x_);
  F_.setAdjSeed(1.0);
  F_.evaluate(0,1);
  F_.getOutput(fk_);
  
  // Gradient of objective
  const DMatrix& gfk = F_.adjSens();
  
  if (monitored("eval_f")){
    cout << "(main loop) x = " << x_ << endl;
    cout << "(main loop) F = " << endl;
    F_.output().printSparse();
  }
  
  if (monitored("eval_grad_f")) {
    cout << "(main loop) x = " << x_ << endl;
    cout << "(main loop) gradF = " << endl;
    gfk.printSparse();
  }
}

} // namespace CasADi
