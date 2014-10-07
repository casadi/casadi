/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


#include "lifted_sqp_internal.hpp"
#include "core/std_vector_tools.hpp"
#include "core/matrix/sparsity_tools.hpp"
#include "core/matrix/matrix_tools.hpp"
#include "core/function/sx_function.hpp"
#include "core/sx/sx_tools.hpp"
#include "core/casadi_calculus.hpp"
#include <ctime>
#include <iomanip>

using namespace std;
namespace casadi{

LiftedSQPInternal::LiftedSQPInternal(const Function& F, const Function& G) : NlpSolverInternal(Function(),F,G){
  casadi_warning("casadi::LiftedSQP has been replaced by casadi::SCPgen. This class will be deleted.");
  addOption("qp_solver",         OT_QPSOLVER,   GenericType(), "The QP solver to be used by the SQP method");
  addOption("qp_solver_options", OT_DICTIONARY, GenericType(), "Options to be passed to the QP solver");
  addOption("max_iter",           OT_INTEGER,    100,           "Maximum number of SQP iterations");
  addOption("max_iter_ls",        OT_INTEGER,    100,           "Maximum number of linesearch iterations");
  addOption("toldx",             OT_REAL   ,    1e-12,         "Stopping criterion for the stepsize");
  addOption("tolgl",             OT_REAL   ,    1e-12,         "Stopping criterion for the Lagrangian gradient");
  addOption("sigma",             OT_REAL   ,    1.0,           "Linesearch parameter");
  addOption("rho",               OT_REAL   ,    0.5,           "Linesearch parameter");
  addOption("mu_safety",         OT_REAL   ,    1.1,           "Safety factor for linesearch mu");
  addOption("eta",               OT_REAL   ,    0.0001,        "Linesearch parameter: See Nocedal 3.4");
  addOption("tau",               OT_REAL   ,    0.2,           "Linesearch parameter");
  addOption("hessian_approximation", OT_STRING, "BFGS",        "BFGS|exact");
  addOption("num_lifted",        OT_INTEGER,   0,              "Number of variables to lift");
  
  // Monitors
  addOption("monitor",      OT_STRINGVECTOR, GenericType(),  "", "eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h|qp", true);
  
  // Set default options
  setOption("generate_jacobian", false);
}


LiftedSQPInternal::~LiftedSQPInternal(){
}

void LiftedSQPInternal::init(){
  // Call the init method of the base class
  NlpSolverInternal::init();

  // Number of lifted variables
  nv = getOption("num_lifted");
  if(verbose_){
    cout << "Initializing SQP method with " << nx_ << " variables and " << ng_ << " constraints." << endl;
    cout << "Lifting " << nv << " variables." << endl;
    if(gauss_newton_){
      cout << "Gauss-Newton objective with " << F_.input().numel() << " terms." << endl;
    }
  }
  
  // Read options
  max_iter_ = getOption("max_iter");
  max_iter_ls_ = getOption("max_iter_ls");
  toldx_ = getOption("toldx");
  tolgl_ = getOption("tolgl");
  sigma_ = getOption("sigma");
  rho_ = getOption("rho");
  mu_safety_ = getOption("mu_safety");
  eta_ = getOption("eta");
  tau_ = getOption("tau");
    
  // Assume SXFunction for now
  SXFunction ffcn = shared_cast<SXFunction>(F_);
  casadi_assert(!ffcn.isNull());
  SXFunction gfcn = shared_cast<SXFunction>(G_);
  casadi_assert(!gfcn.isNull());
  
  // Extract the free variables and split into independent and dependent variables
  SX x = ffcn.inputExpr(0);
  int nx = x.size();
  nu = nx-nv;
  SX u = x[Slice(0,nu)];
  SX v = x[Slice(nu,nu+nv)];

  // Extract the constraint equations and split into constraints and definitions of dependent variables
  SX f1 = ffcn.outputExpr(0);
  int nf1 = f1.numel();
  SX g = gfcn.outputExpr(0);
  int nf2 = g.numel()-nv;
  SX v_eq = g(Slice(0,nv));
  SX f2 = g(Slice(nv,nv+nf2));
  
  // Definition of v
  SX v_def = v_eq + v;

  // Objective function
  SX f;
  
  // Multipliers
  SX lam_x, lam_g, lam_f2;
  if(gauss_newton_){
    
    // Least square objective
    f = inner_prod(f1,f1)/2;
    
  } else {
    
    // Scalar objective function
    f = f1;
    
    // Lagrange multipliers for the simple bounds on u
    SX lam_u = ssym("lam_u",nu);
    
    // Lagrange multipliers for the simple bounds on v
    SX lam_v = ssym("lam_v",nv);
    
    // Lagrange multipliers for the simple bounds on x
    lam_x = vertcat(lam_u,lam_v);

    // Lagrange multipliers corresponding to the definition of the dependent variables
    SX lam_v_eq = ssym("lam_v_eq",nv);

    // Lagrange multipliers for the nonlinear constraints that aren't eliminated
    lam_f2 = ssym("lam_f2",nf2);

    if(verbose_){
      cout << "Allocated intermediate variables." << endl;
    }
    
    // Lagrange multipliers for constraints
    lam_g = vertcat(lam_v_eq,lam_f2);
    
    // Lagrangian function
    SX lag = f + inner_prod(lam_x,x);
    if(!f2.empty()) lag += inner_prod(lam_f2,f2);
    if(!v.empty()) lag += inner_prod(lam_v_eq,v_def);
    
    // Gradient of the Lagrangian
    SX lgrad = casadi::gradient(lag,x);
    if(!v.empty()) lgrad -= vertcat(SX::zeros(nu),lam_v_eq); // Put here to ensure that lgrad is of the form "h_extended -v_extended"
    makeDense(lgrad);
    if(verbose_){
      cout << "Generated the gradient of the Lagrangian." << endl;
    }

    // Condensed gradient of the Lagrangian
    f1 = lgrad[Slice(0,nu)];
    nf1 = nu;
    
    // Gradient of h
    SX v_eq_grad = lgrad[Slice(nu,nu+nv)];
    
    // Reverse lam_v_eq and v_eq_grad
    SX v_eq_grad_reversed = v_eq_grad;
    copy(v_eq_grad.rbegin(),v_eq_grad.rend(),v_eq_grad_reversed.begin());
    SX lam_v_eq_reversed = lam_v_eq;
    copy(lam_v_eq.rbegin(),lam_v_eq.rend(),lam_v_eq_reversed.begin());
    
    // Augment h and lam_v_eq
    v_eq.append(v_eq_grad_reversed);
    v.append(lam_v_eq_reversed);
  }

  // Residual function G
  SXVector G_in(G_NUM_IN);
  G_in[G_X] = x;
  G_in[G_LAM_X] = lam_x;
  G_in[G_LAM_G] = lam_g;

  SXVector G_out(G_NUM_OUT);
  G_out[G_D] = v_eq;
  G_out[G_G] = g;
  G_out[G_F] = f;

  rfcn_ = SXFunction(G_in,G_out);
  rfcn_.setOption("number_of_fwd_dir",0);
  rfcn_.setOption("number_of_adj_dir",0);
  rfcn_.setOption("live_variables",true);
  rfcn_.init();
  if(verbose_){
    cout << "Generated residual function ( " << shared_cast<SXFunction>(rfcn_).getAlgorithmSize() << " nodes)." << endl;
  }
  
  // Difference vector d
  SX d = ssym("d",nv);
  if(!gauss_newton_){
    vector<SX> dg = ssym("dg",nv).data();
    reverse(dg.begin(),dg.end());
    d.append(dg);
  }

  // Substitute out the v from the h
  SX d_def = (v_eq + v)-d;
  SXVector ex(3);
  ex[0] = f1;
  ex[1] = f2;
  ex[2] = f;
  substituteInPlace(v, d_def, ex, false);
  SX f1_z = ex[0];
  SX f2_z = ex[1];
  SX f_z = ex[2];
  
  // Modified function Z
  enum ZIn{Z_U,Z_D,Z_LAM_X,Z_LAM_F2,Z_NUM_IN};
  SXVector zfcn_in(Z_NUM_IN);
  zfcn_in[Z_U] = u;
  zfcn_in[Z_D] = d;
  zfcn_in[Z_LAM_X] = lam_x;
  zfcn_in[Z_LAM_F2] = lam_f2;
  
  enum ZOut{Z_D_DEF,Z_F12,Z_NUM_OUT};
  SXVector zfcn_out(Z_NUM_OUT);
  zfcn_out[Z_D_DEF] = d_def;
  zfcn_out[Z_F12] = vertcat(f1_z,f2_z);
  
  SXFunction zfcn(zfcn_in,zfcn_out);
  zfcn.init();
  if(verbose_){
    cout << "Generated reconstruction function ( " << zfcn.getAlgorithmSize() << " nodes)." << endl;
  }

  // Matrix A and B in lifted Newton
  SX B = zfcn.jac(Z_U,Z_F12);
  SX B1 = B(Slice(0,nf1),Slice(0,B.size2()));
  SX B2 = B(Slice(nf1,B.size1()),Slice(0,B.size2()));
  if(verbose_){
    cout << "Formed B1 (dimension " << B1.size1() << "-by-" << B1.size2() << ", "<< B1.size() << " nonzeros) " <<
    "and B2 (dimension " << B2.size1() << "-by-" << B2.size2() << ", "<< B2.size() << " nonzeros)." << endl;
  }
  
  // Step in u
  SX du = ssym("du",nu);
  SX dlam_f2 = ssym("dlam_f2",lam_f2.sparsity());
  
  SX b1 = f1_z;
  SX b2 = f2_z;
  SX e;
  if(nv > 0){
    
    // Directional derivative of Z
    vector<vector<SX> > Z_fwdSeed(2,zfcn_in);
    vector<vector<SX> > Z_fwdSens(2,zfcn_out);
    vector<vector<SX> > Z_adjSeed;
    vector<vector<SX> > Z_adjSens;
    
    Z_fwdSeed[0][Z_U].setZero();
    Z_fwdSeed[0][Z_D] = -d;
    Z_fwdSeed[0][Z_LAM_X].setZero();
    Z_fwdSeed[0][Z_LAM_F2].setZero();
    
    Z_fwdSeed[1][Z_U] = du;
    Z_fwdSeed[1][Z_D] = -d;
    Z_fwdSeed[1][Z_LAM_X].setZero();
    Z_fwdSeed[1][Z_LAM_F2] = dlam_f2;
    
    zfcn.eval(zfcn_in,zfcn_out,Z_fwdSeed,Z_fwdSens,Z_adjSeed,Z_adjSens);
    
    b1 += Z_fwdSens[0][Z_F12](Slice(0,nf1));
    b2 += Z_fwdSens[0][Z_F12](Slice(nf1,B.size1()));
    e = Z_fwdSens[1][Z_D_DEF];
  }
  if(verbose_){
    cout << "Formed b1 (dimension " << b1.size1() << "-by-" << b1.size2() << ", "<< b1.size() << " nonzeros) " <<
    "and b2 (dimension " << b2.size1() << "-by-" << b2.size2() << ", "<< b2.size() << " nonzeros)." << endl;
  }
  
  // Generate Gauss-Newton Hessian
  if(gauss_newton_){
    b1 = mul(trans(B1),b1);
    B1 = mul(trans(B1),B1);
    if(verbose_){
      cout << "Gauss Newton Hessian (dimension " << B1.size1() << "-by-" << B1.size2() << ", "<< B1.size() << " nonzeros)." << endl;
    }
  }
  
  // Make sure b1 and b2 are dense vectors
  makeDense(b1);
  makeDense(b2);
  
  // Quadratic approximation
  SXVector lfcn_in(LIN_NUM_IN);
  lfcn_in[LIN_X] = x;
  lfcn_in[LIN_D] = d;
  lfcn_in[LIN_LAM_X] = lam_x;
  lfcn_in[LIN_LAM_G] = lam_g;
  
  SXVector lfcn_out(LIN_NUM_OUT);
  lfcn_out[LIN_F1] = b1;
  lfcn_out[LIN_J1] = B1;
  lfcn_out[LIN_F2] = b2;
  lfcn_out[LIN_J2] = B2;
  lfcn_ = SXFunction(lfcn_in,lfcn_out);
//   lfcn_.setOption("verbose",true);
  lfcn_.setOption("number_of_fwd_dir",0);
  lfcn_.setOption("number_of_adj_dir",0);
  lfcn_.setOption("live_variables",true);
  lfcn_.init();
  if(verbose_){
    cout << "Generated linearization function ( " << shared_cast<SXFunction>(lfcn_).getAlgorithmSize() << " nodes)." << endl;
  }
    
  // Step expansion
  SXVector efcn_in(EXP_NUM_IN);
  copy(lfcn_in.begin(),lfcn_in.end(),efcn_in.begin());
  efcn_in[EXP_DU] = du;
  efcn_in[EXP_DLAM_F2] = dlam_f2;
  efcn_ = SXFunction(efcn_in,e);
  efcn_.setOption("number_of_fwd_dir",0);
  efcn_.setOption("number_of_adj_dir",0);
  efcn_.setOption("live_variables",true);
  efcn_.init();
  if(verbose_){
    cout << "Generated step expansion function ( " << shared_cast<SXFunction>(efcn_).getAlgorithmSize() << " nodes)." << endl;
  }
  
  // Current guess for the primal solution
  DMatrix &x_k = output(NLP_SOLVER_X);
  
  // Current guess for the dual solution
  DMatrix &lam_x_k = output(NLP_SOLVER_LAM_X);
  DMatrix &lam_g_k = output(NLP_SOLVER_LAM_G);

  // Allocate a QP solver
  QpSolverCreator qp_solver_creator = getOption("qp_solver");
  qp_solver_ = qp_solver_creator(B1.sparsity(),B2.sparsity());
  
  // Set options if provided
  if(hasSetOption("qp_solver_options")){
    Dictionary qp_solver_options = getOption("qp_solver_options");
    qp_solver_.setOption(qp_solver_options);
  }
  
  // Initialize the QP solver
  qp_solver_.init();
  if(verbose_){
    cout << "Allocated QP solver." << endl;
  }

  // Residual
  d_k_ = DMatrix(d.sparsity(),0);
  
  // Primal step
  dx_k_ = DMatrix(x_k.sparsity());

  // Dual step
  dlam_x_k_ = DMatrix(lam_x_k.sparsity());
  dlam_g_k_ = DMatrix(lam_g_k.sparsity());
  
}

void LiftedSQPInternal::evaluate(int nfdir, int nadir){
  casadi_assert(nfdir==0 && nadir==0);
  checkInitialBounds();
     
  // Objective value
  double f_k = numeric_limits<double>::quiet_NaN();
  
  // Current guess for the primal solution
  DMatrix &x_k = output(NLP_SOLVER_X);
  const DMatrix &x_init = input(NLP_SOLVER_X0);
  copy(x_init.begin(),x_init.end(),x_k.begin());
  
  // Current guess for the dual solution
  DMatrix &lam_x_k = output(NLP_SOLVER_LAM_X);
  DMatrix &lam_g_k = output(NLP_SOLVER_LAM_G);

  // Bounds
  const DMatrix &x_min = input(NLP_SOLVER_LBX);
  const DMatrix &x_max = input(NLP_SOLVER_UBX);
  
  const DMatrix &g_min = input(NLP_SOLVER_LBG);
  const DMatrix &g_max = input(NLP_SOLVER_UBG);
  int k=0;
  
  // Does G depend on the multipliers?
  bool has_lam_x =  !rfcn_.input(G_LAM_X).empty();
  bool has_lam_g =  !rfcn_.input(G_LAM_G).empty();
  bool has_lam_f2 = !efcn_.input(EXP_DLAM_F2).empty();
  
  while(true){
    // Evaluate residual
    rfcn_.setInput(x_k,G_X);
    if(has_lam_x) rfcn_.setInput(lam_x_k,G_LAM_X);
    if(has_lam_g) rfcn_.setInput(lam_g_k,G_LAM_G);
    rfcn_.evaluate();
    rfcn_.getOutput(d_k_,G_D);
    f_k = rfcn_.output(G_F).toScalar();
    const DMatrix& g_k = rfcn_.output(G_G);
    
    // Construct the QP
    lfcn_.setInput(x_k,LIN_X);
    if(has_lam_x) lfcn_.setInput(lam_x_k,LIN_LAM_X);
    if(has_lam_g) lfcn_.setInput(lam_g_k,LIN_LAM_G);
    lfcn_.setInput(d_k_,LIN_D);
    lfcn_.evaluate();
    DMatrix& B1_k = lfcn_.output(LIN_J1);
    const DMatrix& b1_k = lfcn_.output(LIN_F1);
    const DMatrix& B2_k = lfcn_.output(LIN_J2);
    const DMatrix& b2_k = lfcn_.output(LIN_F2);

    // Regularization
    double reg = 0;
    bool regularization = true;
    
    // Check the smallest eigenvalue of the Hessian
    if(regularization && nu==2){
      double a = B1_k.elem(0,0);
      double b = B1_k.elem(0,1);
      double c = B1_k.elem(1,0);
      double d = B1_k.elem(1,1);
      
      // Make sure no not a numbers
      casadi_assert(a==a && b==b && c==c &&  d==d);
      
      // Make sure symmetric
      if(b!=c){
        casadi_assert_warning(fabs(b-c)<1e-10,"Hessian is not symmetric: " << b << " != " << c);
        B1_k.elem(1,0) = c = b;
      }
      
      double eig_smallest = (a+d)/2 - std::sqrt(4*b*c + (a-d)*(a-d))/2;
      double threshold = 1e-8;
      if(eig_smallest<threshold){
        // Regularization
        reg = threshold-eig_smallest;
        std::cerr << "Regularization with " << reg << " to ensure positive definite Hessian." << endl;
        B1_k(0,0) += reg;
        B1_k(1,1) += reg;
      }
    }
    
    
    // Solve the QP
    qp_solver_.setInput(B1_k,QP_H);
    qp_solver_.setInput(b1_k,QP_G);
    qp_solver_.setInput(B2_k,QP_A);
    std::transform(x_min.begin(),x_min.begin()+nu,x_k.begin(),qp_solver_.input(QP_LBX).begin(),std::minus<double>());
    std::transform(x_max.begin(),x_max.begin()+nu,x_k.begin(),qp_solver_.input(QP_UBX).begin(),std::minus<double>());
    std::transform(g_min.begin()+nv,g_min.end(), b2_k.begin(),qp_solver_.input(QP_LBA).begin(),std::minus<double>());
    std::transform(g_max.begin()+nv,g_max.end(), b2_k.begin(),qp_solver_.input(QP_UBA).begin(),std::minus<double>());
    qp_solver_.evaluate();
    const DMatrix& du_k = qp_solver_.output(QP_PRIMAL);
    const DMatrix& dlam_u_k = qp_solver_.output(QP_LAMBDA_X);
    const DMatrix& dlam_f2_k = qp_solver_.output(QP_LAMBDA_A);    
    
    // Expand the step
    for(int i=0; i<LIN_NUM_IN; ++i){
      efcn_.setInput(lfcn_.input(i),i);
    }
    efcn_.setInput(du_k,EXP_DU);
    if(has_lam_f2) efcn_.setInput(dlam_f2_k,EXP_DLAM_F2);
    efcn_.evaluate();
    const DMatrix& dv_k = efcn_.output();
    
    // Expanded primal step
    copy(du_k.begin(),du_k.end(),dx_k_.begin());
    copy(dv_k.begin(),dv_k.begin()+nv,dx_k_.begin()+nu);

    // Expanded dual step
    copy(dlam_u_k.begin(),dlam_u_k.end(),dlam_x_k_.begin());
    copy(dlam_f2_k.begin(),dlam_f2_k.end(),dlam_g_k_.begin()+nv);
    copy(dv_k.rbegin(),dv_k.rbegin()+nv,dlam_g_k_.begin());
    
    // Take a full step
    transform(dx_k_.begin(),dx_k_.end(),x_k.begin(),x_k.begin(),plus<double>());
    copy(dlam_x_k_.begin(),dlam_x_k_.end(),lam_x_k.begin());
    transform(dlam_g_k_.begin(),dlam_g_k_.end(),lam_g_k.begin(),lam_g_k.begin(),plus<double>());

    // Step size
    double norm_step=0;
    for(vector<double>::const_iterator it=dx_k_.begin(); it!=dx_k_.end(); ++it)  norm_step += *it**it;
    if(!gauss_newton_){
      for(vector<double>::const_iterator it=dlam_g_k_.begin(); it!=dlam_g_k_.end(); ++it) norm_step += *it**it;
    }
    norm_step = sqrt(norm_step);
    
    // Constraint violation
    double norm_viol = 0;
    for(int i=0; i<x_k.size(); ++i){
      double d = fmax(x_k.at(i)-x_max.at(i),0.) + fmax(x_min.at(i)-x_k.at(i),0.);
      norm_viol += d*d;
    }
    for(int i=0; i<g_k.size(); ++i){
      double d = fmax(g_k.at(i)-g_max.at(i),0.) + fmax(g_min.at(i)-g_k.at(i),0.);
      norm_viol += d*d;
    }
    norm_viol = sqrt(norm_viol);
    
    // Print progress (including the header every 10 rows)
    if(k % 10 == 0){
      cout << setw(4) << "iter" << setw(20) << "objective" << setw(20) << "norm_step" << setw(20) << "norm_viol" << endl;
    }
    cout   << setw(4) <<     k  << setw(20) <<  f_k        << setw(20) <<  norm_step  << setw(20) <<  norm_viol  << endl;
    
    
    // Check if stopping criteria is satisfied
    if(norm_viol + norm_step < toldx_){
      cout << "Convergence achieved!" << endl;
      break;
    }
    
    // Increase iteration count
    k = k+1;
    
    // Check if number of iterations have been reached
    if(k >= max_iter_){
      cout << "Maximum number of iterations (" << max_iter_ << ") reached" << endl;
      break;
    }
  }
  
  // Store optimal value
  output(NLP_SOLVER_F).set(f_k);
  
  // Save statistics
  stats_["iter_count"] = k;
}

} // namespace casadi
