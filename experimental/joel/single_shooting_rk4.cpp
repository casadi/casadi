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
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSefcn.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <casadi/casadi.hpp>
#include <nonlinear_programming/sqp_method.hpp>
#include <interfaces/qpoases/qpoases_solver.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <casadi/stl_vector_tools.hpp>

using namespace CasADi;
using namespace std;

enum LinIn{LIN_X,LIN_LAM_X,LIN_LAM_HG,LIN_NUM_IN};
enum LinOut{LIN_LHESS,LIN_LGRAD,LIN_GJAC,LIN_GLIN,LIN_OBJ,LIN_CON,LIN_NUM_OUT};

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


void liftedNewton(SXFunction &ffcn, SXFunction &gfcn, const DMatrix& x_min, const DMatrix& x_max, const DMatrix& x_init, const DMatrix& g_min, const DMatrix& g_max, int nv, QPSolverCreator qp_solver_creator, const Dictionary& qp_solver_options){
  cout << "liftedNewton" << endl;
  
  // Options
  double tol = 1e-6;     // Stopping tolerance
  int max_iter = 8;  // Maximum number of iterations

  // Extract the free variables and split into independent and dependent variables
  SXMatrix x = ffcn.inputSX();
  int nx = x.size();
  int nu = nx-nv;
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
  SXMatrixVector G_in,G_out;
  G_in.push_back(u);
  G_in.push_back(v_extended);
  G_in.push_back(lam_u);
  G_in.push_back(lam_g);
  G_in.push_back(lam_v);
  G_out.push_back(h_extended);
  G_out.push_back(lgrad_u);
  G_out.push_back(hg);
  G_out.push_back(f);
  SXFunction G(G_in,G_out);
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
  enum ZIn{Z_U,Z_D,Z_LAM_U,Z_LAM_G,Z_LAM_V,Z_NUM_IN};
  SXMatrixVector zfcn_in(Z_NUM_IN);
  zfcn_in[Z_U] = u;
  zfcn_in[Z_D] = d;
  zfcn_in[Z_LAM_U] = lam_u;
  zfcn_in[Z_LAM_G] = lam_g;
  zfcn_in[Z_LAM_V] = lam_v;
  enum ZOut{Z_D_DEF,Z_LGRAD,Z_G,Z_NUM_OUT};
  SXMatrixVector zfcn_out(Z_NUM_OUT);
  zfcn_out[Z_D_DEF] = d_def;
  zfcn_out[Z_LGRAD] = lgrad_u_z;
  zfcn_out[Z_G] = g_z;
  SXFunction zfcn(zfcn_in,zfcn_out);
  zfcn.init();

  // Matrix A and B in lifted Newton
  SXMatrix B1 = zfcn.jac(0,1);
  SXMatrix B2 = zfcn.jac(0,2);

  // Directional derivative of Z
  vector<vector<SXMatrix> > Z_fwdSeed(2,zfcn_in);
  vector<vector<SXMatrix> > Z_fwdSens(2,zfcn_out);
  vector<vector<SXMatrix> > Z_adjSeed;
  vector<vector<SXMatrix> > Z_adjSens;
  Z_fwdSeed[0][0].setZero();
  Z_fwdSeed[0][1] = -d;
  Z_fwdSeed[0][2].setZero();
  Z_fwdSeed[0][3].setZero();
  Z_fwdSeed[0][4].setZero();
  
  // Step in u
  SXMatrix du = ssym("du",nu);
  SXMatrix dlam_g = ssym("dlam_g",ng);
  
  SXMatrix b1 = lgrad_u_z;
  SXMatrix b2 = g_z;
  if(!v.empty()){
    Z_fwdSeed[1][0] = du;
    Z_fwdSeed[1][1] = -d;
    Z_fwdSeed[1][2].setZero();
    Z_fwdSeed[1][3] = dlam_g;
    Z_fwdSeed[1][4].setZero();
    zfcn.eval(zfcn_in,zfcn_out,Z_fwdSeed,Z_fwdSens,Z_adjSeed,Z_adjSens,true,false);
    b1 += Z_fwdSens[0][1];
    b2 += Z_fwdSens[0][2];
  }
  
  SXMatrixVector lfcn_out(LIN_NUM_OUT);
  lfcn_out[LIN_LHESS] = B1;
  lfcn_out[LIN_LGRAD] = b1;
  lfcn_out[LIN_GJAC] = B2;
  lfcn_out[LIN_GLIN] = b2;
  lfcn_out[LIN_OBJ] = f_z;
  lfcn_out[LIN_CON] = g_z; 
  lfcn_out.insert(lfcn_out.end(),zfcn_out.begin(),zfcn_out.end());
  SXFunction lfcn(zfcn_in,lfcn_out);
  lfcn.init();
  
  // Step expansion
  SXMatrixVector efcn_in = zfcn_in;
  efcn_in.push_back(du);
  efcn_in.push_back(dlam_g);
  SXFunction efcn(efcn_in,Z_fwdSens[1][0]);
  efcn.init();
  
  // Current guess for the primal solution
  DMatrix x_k = x_init;
  
  // Current guess for the dual solution
  DMatrix lam_x_k(lam_x.sparsity(),0);
  DMatrix lam_hg_k(lam_hg.sparsity(),0);

  // Residual
  DMatrix d_k(d.sparsity(),0);
  
  // Initialize node values manually
  copy(x_k.begin(),   x_k.begin()+nu,G.input(0).begin());
  copy(x_k.begin()+nu,x_k.end(),     G.input(1).begin());
  copy(lam_hg_k.begin(),lam_hg_k.begin()+nv,G.input(1).rbegin());
  copy(lam_x_k.begin(),lam_x_k.begin()+nu,G.input(2).begin());
  copy(lam_hg_k.begin()+nv,lam_hg_k.end(),G.input(3).begin());
  G.evaluate();
  G.getOutput(d_k,0);
      
  // Allocate a QP solver (2)
  QPSolver qp_solver = qp_solver_creator(B1.sparsity(),B2.sparsity());
  qp_solver.setOption(qp_solver_options);
  qp_solver.init();

  // Objective
  double f_k;
  
  // Primal step
  DMatrix dx_k(x_k.sparsity());

  // Dual step
  DMatrix dlam_x_k(lam_x_k.sparsity());
  DMatrix dlam_hg_k(lam_hg_k.sparsity());
  
  int k=0;
  while(true){
    
    copy(x_k.begin(),   x_k.begin()+nu,G.input(0).begin());
    copy(x_k.begin()+nu,x_k.end(),     G.input(1).begin());
    copy(lam_hg_k.begin(),lam_hg_k.begin()+nv,G.input(1).rbegin());
    copy(lam_x_k.begin(),lam_x_k.begin()+nu,G.input(2).begin());
    copy(lam_hg_k.begin()+nv,lam_hg_k.end(),G.input(3).begin());

    G.evaluate();
    G.getOutput(d_k,0);
    f_k = G.output(3).toScalar();
    const DMatrix& hg_k = G.output(2);
    
    // Get A_k and Bk
    lfcn.setInput(G.input(0),0);
    lfcn.setInput(d_k,1);
    lfcn.setInput(G.input(2),2);
    lfcn.setInput(G.input(3),3);
    lfcn.setInput(G.input(4),4);
    
    lfcn.evaluate();
    const DMatrix& B1_k = lfcn.output(LIN_LHESS);
    const DMatrix& b1_k = lfcn.output(LIN_LGRAD);
    const DMatrix& B2_k = lfcn.output(LIN_GJAC);
    const DMatrix& b2_k = lfcn.output(LIN_GLIN);

    qp_solver.setInput(B1_k,QP_H);
    qp_solver.setInput(b1_k,QP_G);
    qp_solver.setInput(B2_k,QP_A);
    std::transform(x_min.begin(),x_min.begin()+nu,x_k.begin(),qp_solver.input(QP_LBX).begin(),std::minus<double>());
    std::transform(x_max.begin(),x_max.begin()+nu,x_k.begin(),qp_solver.input(QP_UBX).begin(),std::minus<double>());
    std::transform(g_min.begin()+nv,g_min.end(), b2_k.begin(),qp_solver.input(QP_LBA).begin(),std::minus<double>());
    std::transform(g_max.begin()+nv,g_max.end(), b2_k.begin(),qp_solver.input(QP_UBA).begin(),std::minus<double>());
    
    qp_solver.evaluate();
    const DMatrix& du_k = qp_solver.output(QP_PRIMAL);
    const DMatrix& dlam_u_k = qp_solver.output(QP_LAMBDA_X);
    const DMatrix& dlam_g_k = qp_solver.output(QP_LAMBDA_A);
        
    // Expand the step
    efcn.setInput(lfcn.input(0),0);
    efcn.setInput(lfcn.input(1),1);
    efcn.setInput(lfcn.input(2),2);
    efcn.setInput(lfcn.input(3),3);
    
    efcn.setInput(lfcn.input(4),4);
    efcn.setInput(du_k,5);
    efcn.setInput(dlam_g_k,6);
    efcn.evaluate();
    const DMatrix& dv_k = efcn.output();
    
    // Expanded primal step
    copy(du_k.begin(),du_k.end(),dx_k.begin());
    copy(dv_k.begin(),dv_k.begin()+nv,dx_k.begin()+nu);

    // Expanded dual step
    copy(dlam_u_k.begin(),dlam_u_k.end(),dlam_x_k.begin());
    copy(dlam_g_k.begin(),dlam_g_k.end(),dlam_hg_k.begin()+nv);
    copy(dv_k.rbegin(),dv_k.rbegin()+nv,dlam_hg_k.begin());
    
    // Take a full step
    x_k += dx_k;
    lam_x_k += dlam_x_k;
    lam_hg_k += dlam_hg_k;

    double step_du_k = norm22(dx_k);
    double step_dmug_k = norm22(dlam_hg_k);
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
    if(norm_viol + norm_res  + norm_step < tol){
      cout << "Convergence achieved!" << endl;
      break;
    }
    
    // Increase iteration count
    k = k+1;
    
    // Check if number of iterations have been reached
    if(k >= max_iter){
      cout << "Maximum number of iterations (" << max_iter << ") reached" << endl;
      break;
    }
  }
  
  cout << "optimal cost:    " << f_k << endl;
  cout << "optimal control: " << x_k << endl;
  cout << "multipliers (x): " << lam_x_k << endl;
  cout << "multipliers (hg): " << lam_hg_k  << endl;
}

int main(){
    
  cout << "program started" << endl;
    
  // QP-solver
  QPSolverCreator qp_solver_creator = Interfaces::QPOasesSolver::creator;
  Dictionary qp_solver_options;
  qp_solver_options["printLevel"] = "none";

  // Dimensions
  int nk = 15;  // Number of control segments
  int nj = 100; // Number of integration steps per control segment

  // Control
  SXMatrix u = ssym("u",nk); // control

  // Number of states
  int nx = 3;
  
  // Intermediate variables with initial values and bounds
  SXMatrix v, v_def;
  DMatrix v_init, v_min, v_max;
  
  // Initial values and bounds for the state at the different stages
  DMatrix x_k_init =  DMatrix::zeros(nx);
  DMatrix x_k_min  = -DMatrix::inf(nx); 
  DMatrix x_k_max  =  DMatrix::inf(nx);
  
  // Initial conditions
  DMatrix x_0 = DMatrix::zeros(nx);
  x_0[0] = 0; // x
  x_0[1] = 1; // y
  x_0[2] = 0; // lterm
  
  double tf = 10;
  SX dt = tf/(nj*nk); // time step

  // For all the shooting intervals
  SXMatrix x_k = x_0;
  SXMatrix ode_rhs(x_k.sparsity(),0);
  for(int k=0; k<nk; ++k){
    // Get control
    SX u_k = u[k].at(0);
    
    // Integrate over the interval with Euler forward
    for(int j=0; j<nj; ++j){
      
      // ODE right hand side
      ode_rhs[0] = (1 - x_k[1]*x_k[1])*x_k[0] - x_k[1] + u_k;
      ode_rhs[1] = x_k[0];
      ode_rhs[2] = x_k[0]*x_k[0] + x_k[1]*x_k[1];
      
      // Take a step
      x_k += dt*ode_rhs;
    }
    
    // Lift x
    v_def.append(x_k);
    v_init.append(x_k_init);
    v_min.append(x_k_min);
    v_max.append(x_k_max);
    
    // Allocate intermediate variables
    stringstream ss;
    ss << "v_" << k;
    x_k = ssym(ss.str(),nx);
    v.append(x_k);
  }

  // Objective function
  SXMatrix f = x_k[2] + (tf/nk)*inner_prod(u,u);
  
  // Terminal constraints
  SXMatrix g;
  g.append(x_k[0]);
  g.append(x_k[1]);

  // Bounds on g
  DMatrix g_min = DMatrix::zeros(2);
  DMatrix g_max = DMatrix::zeros(2);

  // Bounds on u and initial condition
  DMatrix u_min  = -0.75*DMatrix::ones(nk);
  DMatrix u_max  =  1.00*DMatrix::ones(nk);
  DMatrix u_init =       DMatrix::zeros(nk);
  
  // Formulate the full-space NLP
  SXFunction ffcn_full(vertcat(u,v),f);
  SXFunction gfcn_full(vertcat(u,v),vertcat(v_def-v,g));
  IpoptSolver ipopt_full(ffcn_full,gfcn_full);
  
  // Set options
  ipopt_full.setOption("generate_hessian",true);
  ipopt_full.setOption("tol",1e-10);
  
  // initialize the solver
  ipopt_full.init();

  // Initial guess and bounds
  DMatrix xv_min = vertcat(u_min,v_min);
  DMatrix xv_max = vertcat(u_max,v_max);
  DMatrix xv_init = vertcat(u_init,v_init);
  DMatrix gv_min = vertcat(DMatrix::zeros(v.size()),g_min);
  DMatrix gv_max = vertcat(DMatrix::zeros(v.size()),g_max);
  
  ipopt_full.setInput(xv_min,NLP_LBX);
  ipopt_full.setInput(xv_max,NLP_UBX);
  ipopt_full.setInput(xv_init,NLP_X_INIT);
  ipopt_full.setInput(gv_min,NLP_LBG);
  ipopt_full.setInput(gv_max,NLP_UBG);

  // Solve the problem
  ipopt_full.solve();
  
  // Print the optimal solution
  cout << "F: optimal cost:    " << ipopt_full.output(NLP_COST).toScalar() << endl;
  cout << "F: optimal control: " << ipopt_full.output(NLP_X_OPT) << endl;
  cout << "F: multipliers (u): " << ipopt_full.output(NLP_LAMBDA_X) << endl;
  cout << "F: multipliers (gb): " << ipopt_full.output(NLP_LAMBDA_G) << endl;
    
  // Substitute in the lifted variables z into the expressions for z_def, f and g
  SXMatrixVector ex(2);
  ex[0] = f;
  ex[1] = g;
  substituteInPlace(v, v_def, ex, false, false);
  f = ex[0];
  g = ex[1];
  
  // Create the NLP
  SXFunction ffcn(u,f); // objective function
  SXFunction gfcn(u,g); // constraint
  SXFunction ifcn(u,v_def);
  
//   // Allocate an NLP solver
//   IpoptSolver ipopt(ffcn,gfcn);
// 
//   // Set options
//   ipopt.setOption("generate_hessian",true);
//   ipopt.setOption("tol",1e-10);
// 
//   // initialize the solver
//   ipopt.init();
// 
//   ipopt.setInput(u_min,NLP_LBX);
//   ipopt.setInput(u_max,NLP_UBX);
//   ipopt.setInput(u_init,NLP_X_INIT);
//   ipopt.setInput(g_min,NLP_LBG);
//   ipopt.setInput(g_max,NLP_UBG);

  // Solve the problem
//   ipopt.solve();
// 
//   // Print the optimal solution
//   cout << "I: optimal cost:    " << ipopt.output(NLP_COST).toScalar() << endl;
//   cout << "I: optimal control: " << ipopt.output(NLP_X_OPT) << endl;
//   cout << "I: multipliers (u): " << ipopt.output(NLP_LAMBDA_X) << endl;
//   cout << "I: multipliers (g): " << ipopt.output(NLP_LAMBDA_G) << endl;
  
  // Solve lifted problem
  liftedNewton(ffcn_full,gfcn_full,xv_min,xv_max,xv_init,gv_min,gv_max,v.size(),qp_solver_creator,qp_solver_options);

  // Solve full-space problem
  liftedNewton(ffcn_full,gfcn_full,xv_min,xv_max,xv_init,gv_min,gv_max,0,qp_solver_creator,qp_solver_options);
  
  
  return 0;
}




