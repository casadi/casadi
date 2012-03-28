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

double norm22(const DMatrix& v){
  double ret = 0;
  for(DMatrix::const_iterator it=v.begin(); it!=v.end(); ++it){
    ret += (*it) * (*it);
  }
  return ret;
}

double norm2(const DMatrix& v){
  return sqrt(norm22(v));
}

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
  
//   // Allocate a QP solver
//   CRSSparsity H_sparsity = getOption("hessian_approximation")=="exact"? H_.output().sparsity() : sp_dense(n,n);
//   CRSSparsity A_sparsity = J_.isNull() ? CRSSparsity(0,n,false) : J_.output().sparsity();
// 
//   QPSolverCreator qp_solver_creator = getOption("qp_solver");
//   qp_solver_ = qp_solver_creator(H_sparsity,A_sparsity);
//   
//   // Set options if provided
//   if(hasSetOption("qp_solver_options")){
//     Dictionary qp_solver_options = getOption("qp_solver_options");
//     qp_solver_.setOption(qp_solver_options);
//   }
// 
//   qp_solver_.init();
  
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
  SXMatrixVector G_in(Z_NUM_IN);
  G_in[Z_U] = u;
  G_in[Z_D] = v_extended;
  G_in[Z_LAM_U] = lam_u;
  G_in[Z_LAM_G] = lam_g;
  G_in[Z_LAM_V] = lam_v;

  SXMatrixVector G_out(G_NUM_OUT);
  G_out[G_H] = h_extended;
  G_out[G_LGRAD] = lgrad_u;
  G_out[G_HG] = hg;
  G_out[G_F] = f;
  G = SXFunction(G_in,G_out);
  G.init();
  
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
  SXMatrixVector zfcn_in(Z_NUM_IN);
  zfcn_in[Z_U] = u;
  zfcn_in[Z_D] = d;
  zfcn_in[Z_LAM_U] = lam_u;
  zfcn_in[Z_LAM_G] = lam_g;
  zfcn_in[Z_LAM_V] = lam_v;
  
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
    Z_fwdSeed[0][Z_LAM_U].setZero();
    Z_fwdSeed[0][Z_LAM_G].setZero();
    Z_fwdSeed[0][Z_LAM_V].setZero();
    
    Z_fwdSeed[1][Z_U] = du;
    Z_fwdSeed[1][Z_D] = -d;
    Z_fwdSeed[1][Z_LAM_U].setZero();
    Z_fwdSeed[1][Z_LAM_G] = dlam_g;
    Z_fwdSeed[1][Z_LAM_V].setZero();
    
    zfcn.eval(zfcn_in,zfcn_out,Z_fwdSeed,Z_fwdSens,Z_adjSeed,Z_adjSens,true,false);
    b1 += Z_fwdSens[0][Z_LGRADG](Slice(0,nu));
    b2 += Z_fwdSens[0][Z_LGRADG](Slice(nu,B.size1()));
    e = Z_fwdSens[1][Z_D_DEF];
  }
  
  // Quadratic approximation
  SXMatrixVector lfcn_out(LIN_NUM_OUT);
  lfcn_out[LIN_LHESS] = B1;
  lfcn_out[LIN_LGRAD] = b1;
  lfcn_out[LIN_GJAC] = B2;
  lfcn_out[LIN_GLIN] = b2;
  lfcn = SXFunction(zfcn_in,lfcn_out);
  lfcn.init();
  
  // Step expansion
  SXMatrixVector efcn_in = zfcn_in;
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
    copy(x_k.begin(),   x_k.begin()+nu,G.input(Z_U).begin());
    copy(x_k.begin()+nu,x_k.end(),     G.input(Z_D).begin());
    copy(lam_hg_k.begin(),lam_hg_k.begin()+nv,G.input(Z_D).rbegin());
    copy(lam_x_k.begin(),lam_x_k.begin()+nu,G.input(Z_LAM_U).begin());
    copy(lam_hg_k.begin()+nv,lam_hg_k.end(),G.input(Z_LAM_G).begin());
    G.evaluate();
    G.getOutput(d_k_,G_H);
    f_k = G.output(G_F).toScalar();
    const DMatrix& hg_k = G.output(G_HG);
    
    // Construct the QP
    lfcn.setInput(G.input(Z_U),Z_U);
    lfcn.setInput(d_k_,Z_D);
    lfcn.setInput(G.input(Z_LAM_U),Z_LAM_U);
    lfcn.setInput(G.input(Z_LAM_G),Z_LAM_G);
    lfcn.setInput(G.input(Z_LAM_V),Z_LAM_V);
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
    for(int i=0; i<Z_NUM_IN; ++i){
      efcn.setInput(lfcn.input(i),i);
    }
    efcn.setInput(du_k,Z_NUM_IN);
    efcn.setInput(dlam_g_k,Z_NUM_IN+1);
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
    x_k += dx_k_;
    lam_x_k += dlam_x_k_;
    lam_hg_k += dlam_hg_k_;

    double step_du_k = norm22(dx_k_);
    double step_dmug_k = norm22(dlam_hg_k_);
    double norm_step = sqrt(step_du_k + step_dmug_k); // add mux
    double norm_res = 0;
    
    // Norm of constraint violation
    double viol_umax = norm22(fmax(x_k-x_max,0));
    double viol_umin = norm22(fmax(x_min-x_k,0));
    double viol_gmax = norm22(fmax(hg_k-g_max,0));
    double viol_gmin = norm22(fmax(g_min-hg_k,0));
    double norm_viol = sqrt(viol_umax + viol_umin + viol_gmax + viol_gmin);
    
    // Print progress (including the header every 10 rows)
    if(k % 10 == 0){
      cout << setw(4) << "iter" << setw(20) << "objective" << setw(20) << "norm_res" << setw(20) << "norm_step" << setw(20) << "norm_viol" << endl;
    }
    cout   << setw(4) <<     k  << setw(20) <<  f_k        << setw(20) <<  norm_res  << setw(20) <<  norm_step  << setw(20) <<  norm_viol  << endl;
    
    
    // Check if stopping criteria is satisfied
    if(norm_viol + norm_res  + norm_step < toldx_){
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
