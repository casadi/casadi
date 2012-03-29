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

#include "lifted_sqp_internal.hpp"
#include "casadi/stl_vector_tools.hpp"
#include "casadi/matrix/sparsity_tools.hpp"
#include "casadi/matrix/matrix_tools.hpp"
#include "casadi/fx/sx_function.hpp"
#include "casadi/sx/sx_tools.hpp"
#include "casadi/casadi_calculus.hpp"
#include <ctime>
#include <iomanip>

using namespace std;
namespace CasADi{

LiftedSQPInternal::LiftedSQPInternal(const FX& F, const FX& G, const FX& H, const FX& J) : NLPSolverInternal(F,G,H,J){
  casadi_warning("The lifted SQP method is under development");
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
  addOption("hessian_approximation", OT_STRING, "BFGS",        "BFGS|exact");
  addOption("num_lifted",        OT_INTEGER,   0,              "Number of variables to lift");
  
  // Monitors
  addOption("monitor",      OT_STRINGVECTOR, GenericType(),  "", "eval_f|eval_g|eval_jac_g|eval_grad_f|eval_h|qp", true);
}


LiftedSQPInternal::~LiftedSQPInternal(){
}

void LiftedSQPInternal::init(){
  // Call the init method of the base class
  NLPSolverInternal::init();
    
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
    
  // Assume SXFunction for now
  SXFunction ffcn = shared_cast<SXFunction>(F_);
  casadi_assert(!ffcn.isNull());
  SXFunction gfcn = shared_cast<SXFunction>(G_);
  casadi_assert(!gfcn.isNull());
 
  // Number of lifted variables
  nv = getOption("num_lifted");
  
  // Extract the free variables and split into independent and dependent variables
  SXMatrix x = ffcn.inputSX();
  int nx = x.size();
  nu = nx-nv;
  SXMatrix u = x[Slice(0,nu)];
  SXMatrix v = x[Slice(nu,nu+nv)];

  // Extract the constraint equations and split into constraints and definitions of dependent variables
  SXMatrix f = ffcn.outputSX();
  SXMatrix hg = gfcn.outputSX();
  int ng = hg.numel()-nv;
  SXMatrix h = hg(Slice(0,nv));
  SXMatrix g = hg(Slice(nv,nv+ng));
  
  // Definition of v
  SXMatrix v_def = h+v;

  // Lagrange multipliers for the simple bounds on u
  SXMatrix lam_u = ssym("lam_u",nu);
  
  // Lagrange multipliers for the simple bounds on v
  SXMatrix lam_v = ssym("lam_v",nv);
  
  // Lagrange multipliers for the simple bounds on x
  SXMatrix lam_x = vertcat(lam_u,lam_v);

  // Lagrange multipliers corresponding to the definition of the dependent variables
  SXMatrix lam_h = ssym("lam_h",nv);

  // Lagrange multipliers for the nonlinear constraints that aren't eliminated
  SXMatrix lam_g = ssym("lam_g",ng);

  // Lagrange multipliers for constraints
  SXMatrix lam_hg = vertcat(lam_h,lam_g);
  
  // Lagrangian function
  SXMatrix lag = f + inner_prod(lam_x,x);
  if(!g.empty()) lag += inner_prod(lam_g,g);
  if(!v.empty()) lag += inner_prod(lam_h,v_def);
  
  // Gradient of the lagrangian
  SXMatrix lgrad = gradient(lag,x);
  if(!v.empty()) lgrad -= vertcat(SXMatrix::zeros(nu),lam_h); // Put here to ensure that lgrad is of the form "h_extended -v_extended"
  makeDense(lgrad);

  // Condensed gradient of the Lagrangian
  SXMatrix lgrad_u = lgrad[Slice(0,nu)];
  
  // Gradient of h
  SXMatrix hgrad = lgrad[Slice(nu,nu+nv)];
  
  // Reverse lam_h and hgrad
  SXMatrix hgrad_reversed = hgrad;
  copy(hgrad.rbegin(),hgrad.rend(),hgrad_reversed.begin());
  SXMatrix lam_h_reversed = lam_h;
  copy(lam_h.rbegin(),lam_h.rend(),lam_h_reversed.begin());
  
  // Augment h and lam_h
  SXMatrix h_extended = vertcat(h,hgrad_reversed);
  SXMatrix v_extended = vertcat(v,lam_h_reversed);

  // Residual function G
  SXMatrixVector G_in(G_NUM_IN);
  G_in[G_X] = x;
  G_in[G_LAM_X] = lam_x;
  G_in[G_LAM_HG] = lam_hg;

  SXMatrixVector G_out(G_NUM_OUT);
  G_out[G_H] = h_extended;
  G_out[G_LGRAD] = lgrad_u;
  G_out[G_HG] = hg;
  G_out[G_F] = f;

  rfcn = SXFunction(G_in,G_out);
  rfcn.init();
  
  // Difference vector d
  SXMatrix d = ssym("d",2*nv);
  SXMatrix dg = ssym("dg",nv);
  copy(dg.begin(),dg.end(),d.rbegin());
  
  // Substitute out the v from the h
  SXMatrix d_def = (h_extended+v_extended)-d;
  SXMatrixVector ex(3);
  ex[0] = lgrad_u;
  ex[1] = f;
  ex[2] = g;
  substituteInPlace(v_extended, d_def, ex, false, false);
  SXMatrix lgrad_u_z = ex[0];
  SXMatrix f_z = ex[1];
  SXMatrix g_z = ex[2];
  
  // Modified function Z
  enum ZIn{Z_U,Z_D,Z_LAM_X,Z_LAM_G,Z_NUM_IN};
  SXMatrixVector zfcn_in(Z_NUM_IN);
  zfcn_in[Z_U] = u;
  zfcn_in[Z_D] = d;
  zfcn_in[Z_LAM_X] = lam_x;
  zfcn_in[Z_LAM_G] = lam_g;
  
  enum ZOut{Z_D_DEF,Z_LGRADG,Z_NUM_OUT};
  SXMatrixVector zfcn_out(Z_NUM_OUT);
  zfcn_out[Z_D_DEF] = d_def;
  zfcn_out[Z_LGRADG] = vertcat(lgrad_u_z,g_z);
  
  SXFunction zfcn(zfcn_in,zfcn_out);
  zfcn.init();

  // Matrix A and B in lifted Newton
  SXMatrix B = zfcn.jac(Z_U,Z_LGRADG);
  SXMatrix B1 = B(Slice(0,nu),Slice(0,B.size2()));
  SXMatrix B2 = B(Slice(nu,B.size1()),Slice(0,B.size2()));
  
  // Step in u
  SXMatrix du = ssym("du",nu);
  SXMatrix dlam_g = ssym("dlam_g",ng);
  
  SXMatrix b1 = lgrad_u_z;
  SXMatrix b2 = g_z;
  SXMatrix e;
  if(!v.empty()){
    // Directional derivative of Z
    vector<vector<SXMatrix> > Z_fwdSeed(2,zfcn_in);
    vector<vector<SXMatrix> > Z_fwdSens(2,zfcn_out);
    vector<vector<SXMatrix> > Z_adjSeed;
    vector<vector<SXMatrix> > Z_adjSens;
    Z_fwdSeed[0][Z_U].setZero();
    Z_fwdSeed[0][Z_D] = -d;
    Z_fwdSeed[0][Z_LAM_X].setZero();
    Z_fwdSeed[0][Z_LAM_G].setZero();
    
    Z_fwdSeed[1][Z_U] = du;
    Z_fwdSeed[1][Z_D] = -d;
    Z_fwdSeed[1][Z_LAM_X].setZero();
    Z_fwdSeed[1][Z_LAM_G] = dlam_g;
    
    zfcn.eval(zfcn_in,zfcn_out,Z_fwdSeed,Z_fwdSens,Z_adjSeed,Z_adjSens,true,false);
    b1 += Z_fwdSens[0][Z_LGRADG](Slice(0,nu));
    b2 += Z_fwdSens[0][Z_LGRADG](Slice(nu,B.size1()));
    e = Z_fwdSens[1][Z_D_DEF];
  }
  
  // Quadratic approximation
  SXMatrixVector lfcn_in(LIN_NUM_IN);
  lfcn_in[LIN_X] = x;
  lfcn_in[LIN_D] = d;
  lfcn_in[LIN_LAM_X] = lam_x;
  lfcn_in[LIN_LAM_HG] = lam_hg;
  
  SXMatrixVector lfcn_out(LIN_NUM_OUT);
  lfcn_out[LIN_LHESS] = B1;
  lfcn_out[LIN_LGRAD] = b1;
  lfcn_out[LIN_GJAC] = B2;
  lfcn_out[LIN_GLIN] = b2;
  lfcn = SXFunction(lfcn_in,lfcn_out);
  lfcn.init();
    
  // Step expansion
  SXMatrixVector efcn_in = lfcn_in;
  efcn_in.push_back(du);
  efcn_in.push_back(dlam_g);
  efcn = SXFunction(efcn_in,e);
  efcn.init();
  
  // Current guess for the primal solution
  DMatrix &x_k = output(NLP_X_OPT);
  
  // Current guess for the dual solution
  DMatrix &lam_x_k = output(NLP_LAMBDA_X);
  DMatrix &lam_hg_k = output(NLP_LAMBDA_G);

  // Allocate a QP solver
  QPSolverCreator qp_solver_creator = getOption("qp_solver");
  qp_solver_ = qp_solver_creator(B1.sparsity(),B2.sparsity());
  
  // Set options if provided
  if(hasSetOption("qp_solver_options")){
    Dictionary qp_solver_options = getOption("qp_solver_options");
    qp_solver_.setOption(qp_solver_options);
  }
  
  // Initialize the QP solver
  qp_solver_.init();

  // Residual
  d_k_ = DMatrix(d.sparsity(),0);
  
  // Primal step
  dx_k_ = DMatrix(x_k.sparsity());

  // Dual step
  dlam_x_k_ = DMatrix(lam_x_k.sparsity());
  dlam_hg_k_ = DMatrix(lam_hg_k.sparsity());
  
}

void LiftedSQPInternal::evaluate(int nfdir, int nadir){
  casadi_assert(nfdir==0 && nadir==0);
  checkInitialBounds();
     
  // Objective value
  double f_k = numeric_limits<double>::quiet_NaN();
  
  // Current guess for the primal solution
  DMatrix &x_k = output(NLP_X_OPT);
  const DMatrix &x_init = input(NLP_X_INIT);
  copy(x_init.begin(),x_init.end(),x_k.begin());
  
  // Current guess for the dual solution
  DMatrix &lam_x_k = output(NLP_LAMBDA_X);
  DMatrix &lam_hg_k = output(NLP_LAMBDA_G);

  // Bounds
  const DMatrix &x_min = input(NLP_LBX);
  const DMatrix &x_max = input(NLP_UBX);
  
  const DMatrix &g_min = input(NLP_LBG);
  const DMatrix &g_max = input(NLP_UBG);
  
  int k=0;
  while(true){
    // Evaluate residual
    rfcn.setInput(x_k,G_X);
    rfcn.setInput(lam_x_k,G_LAM_X);
    rfcn.setInput(lam_hg_k,G_LAM_HG);
    rfcn.evaluate();
    rfcn.getOutput(d_k_,G_H);
    f_k = rfcn.output(G_F).toScalar();
    const DMatrix& hg_k = rfcn.output(G_HG);
    
    // Construct the QP
    lfcn.setInput(x_k,LIN_X);
    lfcn.setInput(lam_x_k,LIN_LAM_X);
    lfcn.setInput(lam_hg_k,LIN_LAM_HG);
    lfcn.setInput(d_k_,LIN_D);
    lfcn.evaluate();
    const DMatrix& B1_k = lfcn.output(LIN_LHESS);
    const DMatrix& b1_k = lfcn.output(LIN_LGRAD);
    const DMatrix& B2_k = lfcn.output(LIN_GJAC);
    const DMatrix& b2_k = lfcn.output(LIN_GLIN);

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
    const DMatrix& dlam_g_k = qp_solver_.output(QP_LAMBDA_A);
            
    // Expand the step
    for(int i=0; i<LIN_NUM_IN; ++i){
      efcn.setInput(lfcn.input(i),i);
    }
    efcn.setInput(du_k,LIN_NUM_IN);
    efcn.setInput(dlam_g_k,LIN_NUM_IN+1);
    efcn.evaluate();
    const DMatrix& dv_k = efcn.output();
    
    // Expanded primal step
    copy(du_k.begin(),du_k.end(),dx_k_.begin());
    copy(dv_k.begin(),dv_k.begin()+nv,dx_k_.begin()+nu);

    // Expanded dual step
    copy(dlam_u_k.begin(),dlam_u_k.end(),dlam_x_k_.begin());
    copy(dlam_g_k.begin(),dlam_g_k.end(),dlam_hg_k_.begin()+nv);
    copy(dv_k.rbegin(),dv_k.rbegin()+nv,dlam_hg_k_.begin());
    
    // Take a full step
    transform(dx_k_.begin(),dx_k_.end(),x_k.begin(),x_k.begin(),plus<double>());
    copy(dlam_x_k_.begin(),dlam_x_k_.end(),lam_x_k.begin());
    transform(dlam_hg_k_.begin(),dlam_hg_k_.end(),lam_hg_k.begin(),lam_hg_k.begin(),plus<double>());

    // Step size
    double norm_step=0;
    for(vector<double>::const_iterator it=dx_k_.begin(); it!=dx_k_.end(); ++it)  norm_step += *it**it;
    for(vector<double>::const_iterator it=dlam_hg_k_.begin(); it!=dlam_hg_k_.end(); ++it) norm_step += *it**it;
    norm_step = sqrt(norm_step);
    
    // Constraint violation
    double norm_viol = 0;
    for(int i=0; i<x_k.size(); ++i){
      double d = fmax(x_k.at(i)-x_max.at(i),0.) + fmax(x_min.at(i)-x_k.at(i),0.);
      norm_viol += d*d;
    }
    for(int i=0; i<hg_k.size(); ++i){
      double d = fmax(hg_k.at(i)-g_max.at(i),0.) + fmax(g_min.at(i)-hg_k.at(i),0.);
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
    if(k >= maxiter_){
      cout << "Maximum number of iterations (" << maxiter_ << ") reached" << endl;
      break;
    }
  }
  
  // Store optimal value
  output(NLP_COST).set(f_k);
}

} // namespace CasADi
