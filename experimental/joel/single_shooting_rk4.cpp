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
  
void liftedNewton(SXFunction &F1, SXFunction &F2, SXFunction &ifcn, bool manual_init, bool gauss_newton, QPSolverCreator qp_solver_creator, const Dictionary& qp_solver_options,
  const DMatrix& u_min, const DMatrix& u_max, const DMatrix& u_init, const DMatrix& g_min, const DMatrix& g_max, DMatrix x_init){
  
  // Options
  double tol = 1e-6;     // Stopping tolerance
  int max_iter = 30;  // Maximum number of iterations

  // Extract the free variable and expressions for F and xdef
  SXMatrix u = F1.inputSX();
  SXMatrix f1 = F1.outputSX();
  SXMatrix f2 = F2.outputSX();
  SXMatrix xdef = ifcn.outputSX();

  // Lifted variables
  int nx = xdef.size();
  SXMatrix x = ssym("x",nx);

  // Lagrange multipliers
  SXMatrix lam = ssym("lam",f2.size());
  
  // Gradient of the lagrangian
  SXMatrix lgrad1 = gradient(f1 + inner_prod(lam,f2),u);
  
  SXMatrixVector lgradfcn_in(2);
  lgradfcn_in[0] = u;
  lgradfcn_in[1] = lam;
  SXFunction lgradfcn(lgradfcn_in,lgrad1);
  lgradfcn.init();
  
  // Substitute in the lifted variables x into the expressions for xdef, F1 and F2
  SXMatrixVector ex(2);
  ex[0] = f1;
  ex[1] = f2;
  substituteInPlace(x, xdef, ex, true, false);
  f1 = ex[0];
  f2 = ex[1];

  // Formulate the full-space problem
  SXMatrix ux = vertcat(u,x);
  SXMatrix f2_full = vertcat(f2,x-xdef);
  
  SXFunction ffcn_full(ux,f1);
  SXFunction gfcn_full(ux,f2_full);
  IpoptSolver ipopt_full(ffcn_full,gfcn_full);

//   x.printDense();
//   xdef.printDense();
//   f1.printDense();
//   return;
  
  
  // Set options
  ipopt_full.setOption("generate_hessian",true);
  ipopt_full.setOption("tol",1e-10);
  
  // initialize the solver
  ipopt_full.init();

  // Initial guess and bounds
  DMatrix ux_init = vertcat(u_init,x_init);
  DMatrix ux_min = vertcat(u_min,-DMatrix::inf(nx));
  DMatrix ux_max = vertcat(u_max, DMatrix::inf(nx));
  DMatrix g_min_full = vertcat(g_min,DMatrix::zeros(nx));
  DMatrix g_max_full = vertcat(g_max,DMatrix::zeros(nx));
  
  ipopt_full.setInput(ux_min,NLP_LBX);
  ipopt_full.setInput(ux_max,NLP_UBX);
  ipopt_full.setInput(ux_init,NLP_X_INIT);
  ipopt_full.setInput(g_min_full,NLP_LBG);
  ipopt_full.setInput(g_max_full,NLP_UBG);

  
  
//   cout << x.numel() << endl;
//   
// 
//   cout << u.numel() << endl;
//   cout << f2.numel() << endl;
//   cout << u_min << endl;
//   cout << u_max << endl;
//   cout << g_min << endl;
//   cout << g_max << endl;
// 
//   cout << ux.numel() << endl;
//   cout << f2_full.numel() << endl;
//   cout << ux_min << endl;
//   cout << ux_max << endl;
//   cout << g_min_full << endl;
//   cout << g_max_full << endl;

  //return;
  
  
  
  
  // Solve the problem
  ipopt_full.solve();
  
  // Print the optimal solution
  cout << "F: optimal cost:    " << ipopt_full.output(NLP_COST).toScalar() << endl;
  cout << "F: optimal control: " << ipopt_full.output(NLP_X_OPT) << endl;
  cout << "F: multipliers (u): " << ipopt_full.output(NLP_LAMBDA_X) << endl;
  cout << "F: multipliers (g): " << ipopt_full.output(NLP_LAMBDA_G) << endl;
  
  
  return;
  
  
  
  
  DMatrix x_init2 = x_init;
  
  SXMatrix mux, mug;
  if(!gauss_newton){
    // if Gauss-Newton no multipliers needed
    // If SQP, get the gradient of the lagrangian now
    
    // Derivatives of lifted variables
    SXMatrix xdot = ssym("xdot",x.size());
    
    // Lagrange multipliers
    mux = ssym("mux",u.size1());
    mug = ssym("mug",f2.size1());

    // Gradient of the Lagrangian
    SXMatrix xu = vertcat(u,x);
    SXMatrix lgrad = gradient(f1 + inner_prod(mug,f2) + inner_prod(xdot,xdef),xu);
    makeDense(lgrad);

    // Gradient of the Lagrangian
    f1 = lgrad(range(u.size1()),0); // + mux // NOTE: What about the mux term?

    // Definition of xdot
    SXMatrix xdotdef = lgrad(range(u.size1(),lgrad.size1()),0);
    
    // Reverse direction of x
    SXMatrix xdot_reversed = xdot;
    for(int k=0; k<xdot.size(); ++k){
      xdot_reversed[k] = xdot[-1-k];
    }
    xdot = xdot_reversed;
    
    SXMatrix xdotdef_reversed = xdotdef;
    for(int k=0; k<xdotdef.size(); ++k){
      xdotdef_reversed[k] = xdotdef[-1-k];
    }
    xdotdef = xdotdef_reversed;
    
    // Append to xdef and x
    x.append(xdot);
    xdef.append(xdotdef);
    
    // Append to initial guess
    x_init.append(DMatrix::zeros(x_init.size()));
  }

  // Residual function G
  SXMatrixVector G_in,G_out;
  G_in.push_back(u);
  G_in.push_back(x);
  G_in.push_back(mux);
  G_in.push_back(mug);
  G_out.push_back(xdef-x);
  G_out.push_back(f1);
  G_out.push_back(f2);
  SXFunction G(G_in,G_out);
  G.init();
  
  // Difference vector d
  SXMatrix d = ssym("d",xdef.size1());

  // Substitute out the x from the zdef
  SXMatrix z = xdef-d;
  ex[0] = f1;
  ex[1] = f2;
  substituteInPlace(x, z, ex, false, false);
  f1 = ex[0];
  f2 = ex[1];

  // Modified function Z
  SXMatrixVector Z_in,Z_out;
  Z_in.push_back(u);
  Z_in.push_back(d);
  Z_in.push_back(mux);
  Z_in.push_back(mug);
  Z_out.push_back(z);
  Z_out.push_back(f1);
  Z_out.push_back(f2);
  SXFunction Z(Z_in,Z_out);
  Z.init();

  // Matrix A and B in lifted Newton
  SXMatrix A = Z.jac(0,0);
  SXMatrix B1 = Z.jac(0,1);
  SXMatrix B2 = Z.jac(0,2);

  SXMatrixVector AB_out;
  AB_out.push_back(A);
  AB_out.push_back(B1);
  AB_out.push_back(B2);
  SXFunction AB(Z_in,AB_out);
  AB.init();
  
  // Variables
  DMatrix u_k = u_init;
  DMatrix x_k = x_init;
  DMatrix d_k(x.sparsity(),0);
  DMatrix mux_k(mux.sparsity(),0);
  DMatrix mug_k(mug.sparsity(),0);
  DMatrix dmux_k(mux.sparsity(),0);
  DMatrix dmug_k(mug.sparsity(),0);
  DMatrix f1_k(f1.sparsity(),0);
  DMatrix f2_k(f2.sparsity(),0);

    
  if(manual_init){
    // Initialize node values manually
    G.setInput(u_k,0);
    G.setInput(x_k,1);
    G.setInput(mux_k,2);
    G.setInput(mug_k,3);
    G.evaluate();
    G.getOutput(d_k,0);
    G.getOutput(f1_k,1); // mux is zero (initial multiplier guess)
    G.getOutput(f2_k,2);
  } else {
    // Initialize x0 by function evaluation
    Z.setInput(u_k,0);
    Z.setInput(d_k,1);
    Z.setInput(mux_k,2);
    Z.setInput(-mug_k,3);
    Z.evaluate();
    Z.getOutput(x_k,0);
    Z.getOutput(f1_k,1); // mux is zero (initial multiplier guess)
    Z.getOutput(f2_k,2);
  }
    
  // Zero seeds
  DMatrix u0seed(u.sparsity(),0);
  DMatrix d0seed(d.sparsity(),0);
  DMatrix mux0seed(mux.sparsity(),0);
  DMatrix mug0seed(mug.sparsity(),0);

  // QP solver
  QPSolver qp_solver;
  
  // Iterate
  int k = 0;
  while(true){
    
    // Get A_k and Bk
    AB.setInput(u_k,0);
    AB.setInput(d_k,1);
    AB.setInput(mux_k,2);
    AB.setInput(-mug_k,3);
    AB.evaluate();
    const DMatrix& A_k = AB.output(0);
    const DMatrix& B1_k = AB.output(1); // NOTE: # mux dissappears (constant term)
    const DMatrix& B2_k = AB.output(2);
    
    // Get a_k and b_k
    Z.setInput(u_k,0);
    Z.setInput(d_k,1);
    Z.setInput(mux_k,2);
    Z.setInput(-mug_k,3);
    Z.setFwdSeed(u0seed,0);
    Z.setFwdSeed(d_k,1);
    Z.setFwdSeed(mux0seed,2);
    Z.setFwdSeed(mug0seed,3);

    Z.evaluate(1,0);
    //Z.getOutput(x_k,0);
    Z.getOutput(f1_k,1);
    Z.getOutput(f2_k,2);
    DMatrix a_k = -Z.fwdSens(0);
    DMatrix b1_k = f1_k-Z.fwdSens(1); // mux disappears from Z (constant term)
    DMatrix b2_k = f2_k-Z.fwdSens(2);

    DMatrix H,g,A,a;
    if(gauss_newton){
      // Gauss-Newton Hessian
      H = mul(trans(B1_k),B1_k);
      g = mul(trans(B1_k),b1_k);
      A = B2_k;
      a = b2_k;
    } else {
      // Exact Hessian
      H = B1_k;
      g = b1_k; // +/- mux_k here?
      A = B2_k;
      a = b2_k;
    }

    if(k==0){
      // Allocate a QP solver
      qp_solver = qp_solver_creator(H.sparsity(),A.sparsity());
      qp_solver.setOption(qp_solver_options);
      qp_solver.init();
    }

    // Formulate the QP
    qp_solver.setInput(H,QP_H);
    qp_solver.setInput(g,QP_G);
    qp_solver.setInput(A,QP_A);
    qp_solver.setInput(u_min-u_k,QP_LBX);
    qp_solver.setInput(u_max-u_k,QP_UBX);
    qp_solver.setInput(g_min-a,QP_LBA);
    qp_solver.setInput(g_max-a,QP_UBA);

//     cout << "H   = " << qp_solver.input(QP_H) << endl;
//     cout << "g   = " << qp_solver.input(QP_G) << endl;
//     cout << "A   = " << qp_solver.input(QP_A) << endl;
//     cout << "lb  = " << qp_solver.input(QP_LBX) << endl;
//     cout << "ub  = " << qp_solver.input(QP_UBX) << endl;
//     cout << "lbA = " << qp_solver.input(QP_LBA) << endl;
//     cout << "ubA = " << qp_solver.input(QP_UBA) << endl;

    // Solve the QP
    qp_solver.evaluate();

    // Get the primal solution
    DMatrix du_k = qp_solver.output(QP_PRIMAL);
    
    // Get the dual solution
    if(!gauss_newton){
      qp_solver.getOutput(dmux_k,QP_DUAL_X);
      qp_solver.getOutput(dmug_k,QP_DUAL_A);
    }
    
    // Calculate the step in x
    Z.setFwdSeed(du_k,0);
    Z.setFwdSeed(d0seed,1); // could the a_k term be moved here?
    Z.setFwdSeed(dmux_k,2);
    Z.setFwdSeed(-dmug_k,3);
    Z.evaluate(1,0);
    DMatrix dx_k = Z.fwdSens(0);
        
    // Take a full step
    u_k =     u_k +   du_k;
    x_k =     x_k +    a_k + dx_k;
    mug_k = mug_k + dmug_k;
    mux_k = mux_k + dmux_k;

    // Call algorithm 2 to obtain new d_k and fk
    G.setInput(u_k,0);
    G.setInput(x_k,1);
    G.setInput(mux_k,2);
    G.setInput(-mug_k,3);
    G.evaluate();
    G.getOutput(d_k,0);
    G.getOutput(f1_k,1); // mux?
    G.getOutput(f2_k,2);

    // Norm of residual error
    double norm_res = norm2(d_k);
    
    // Norm of step size
    double step_du_k = norm22(du_k);
    double step_dmug_k = norm22(dmug_k);
    double norm_step = sqrt(step_du_k + step_dmug_k); // add mux

    // Norm of constraint violation
    double viol_umax = norm22(fmax(u_k-u_max,0));
    double viol_umin = norm22(fmax(u_min-u_k,0));
    double viol_gmax = norm22(fmax(f2_k-g_max,0));
    double viol_gmin = norm22(fmax(g_min-f2_k,0));
    double norm_viol = sqrt(viol_umax + viol_umin + viol_gmax + viol_gmin);

    // Print progress (including the header every 10 rows)
    if(k % 10 == 0){
      cout << setw(4) << "iter" << setw(20) << "norm_res" << setw(20) << "norm_step" << setw(20) << "norm_viol" << endl;
    }
    cout   << setw(4) <<     k  << setw(20) <<  norm_res  << setw(20) <<  norm_step  << setw(20) <<  norm_viol  << endl;
    
    // Check if stopping criteria is satisfied
    if(norm_viol + norm_res  + norm_step < tol){
      cout << "Convergens achieved!" << endl;
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
  
  F1.init();
  F1.setInput(u_k);
  F1.evaluate();
  double cost;
  F1.getOutput(cost);
  
  cout << "optimal cost:    " << cost << endl;
  cout << "optimal control: " << u_k << endl;
  cout << "multipliers (u): " << mux_k << endl;
  cout << "multipliers (g): " << mug_k << endl;
  

  lgradfcn.setInput(u_k,0);
//   lgradfcn.setInput(x_init,0);

  return;

#if 0
  {  
  
  // Dimensions
  int nu = 20;  // Number of control segments
  int nj = 100; // Number of integration steps per control segment

  // optimization variable
  vector<SX> u = ssym("u",nu).data(); // control

  SX s_0 = 0; // initial position
  SX v_0 = 0; // initial speed
  SX m_0 = 1; // initial mass
  
  SX dt = 10.0/(nj*nu); // time step
  SX alpha = 0.05; // friction
  SX beta = 0.1; // fuel consumption rate

  // Trajectory
  vector<SX> s_traj(nu), v_traj(nu), m_traj(nu);

  // Integrate over the interval with Euler forward
  SX s = s_0, v = v_0, m = m_0;
  for(int k=0; k<nu; ++k){
    for(int j=0; j<nj; ++j){
      s += dt*v;
      v += dt / m * (u[k]- alpha * v*v);
      m += -dt * beta*u[k]*u[k];
    }
    s_traj[k] = s;
    v_traj[k] = v;
    m_traj[k] = m;
  }

  // Objective function
  SX f = 0;
  for(int i=0; i<u.size(); ++i)
    f += u[i]*u[i];
    
  // Terminal constraints
  vector<SX> g(2);
  g[0] = s;
  g[1] = v;
  g.insert(g.end(),v_traj.begin(),v_traj.end());
  
  // Create the NLP
  SXFunction ffcn(u,f); // objective function
  SXFunction gfcn(u,g); // constraint
  
  // Allocate an NLP solver
  SQPMethod solver(ffcn,gfcn);

  // Set options
  solver.setOption("qp_solver",Interfaces::QPOasesSolver::creator);
  solver.setOption("generate_hessian",true);

  // initialize the solver
  solver.init();

  // Bounds on u and initial condition
  vector<double> umin(nu), umax(nu), uinit(nu);
  for(int i=0; i<nu; ++i){
    umin[i] = -10;
    umax[i] =  10;
    uinit[i] = 0.4;
  }
  solver.setInput(umin,NLP_LBX);
  solver.setInput(umax,NLP_UBX);
  solver.setInput(uinit,NLP_X_INIT);
  
  // Bounds on g
  vector<double> gmin(2), gmax(2);
  gmin[0] = gmax[0] = 10;
  gmin[1] = gmax[1] =  0;
  gmin.resize(2+nu, -numeric_limits<double>::infinity());
  gmax.resize(2+nu, 1.1);
  
  solver.setInput(gmin,NLP_LBG);
  solver.setInput(gmax,NLP_UBG);

  // Solve the problem
  solver.solve();
  
  // Print the optimal cost
  double cost;
  solver.getOutput(cost,NLP_COST);
  cout << "optimal cost: " << cost << endl;

  // Print the optimal solution
  vector<double> uopt(nu);
  solver.getOutput(uopt,NLP_X_OPT);
  cout << "optimal control: " << uopt << endl;

  // Get the state trajectory
  vector<double> sopt(nu), vopt(nu), mopt(nu);
  vector<Matrix<SX> > xfcn_out(3);
  xfcn_out[0] = s_traj;
  xfcn_out[1] = v_traj;
  xfcn_out[2] = m_traj;
  SXFunction xfcn(u,xfcn_out);
  xfcn.init();
  xfcn.setInput(uopt);
  xfcn.evaluate();
  xfcn.getOutput(sopt,0);
  xfcn.getOutput(vopt,1);
  xfcn.getOutput(mopt,2);
  cout << "position: " << sopt << endl;
  cout << "velocity: " << vopt << endl;
  cout << "mass:     " << mopt << endl;

  // Create Matlab script to plot the solution
  ofstream file;
  string filename = "rocket_ipopt_results.m";
  file.open(filename.c_str());
  file << "% Results file from " __FILE__ << endl;
  file << "% Generated " __DATE__ " at " __TIME__ << endl;
  file << endl;
  file << "cost = " << cost << ";" << endl;
  file << "u = " << uopt << ";" << endl;

  // Save results to file
  file << "t = linspace(0,10.0," << nu << ");"<< endl;
  file << "s = " << sopt << ";" << endl;
  file << "v = " << vopt << ";" << endl;
  file << "m = " << mopt << ";" << endl;
  
  // Finalize the results file
  file << endl;
  file << "% Plot the results" << endl;
  file << "figure(1);" << endl;
  file << "clf;" << endl << endl;
  
  file << "subplot(2,2,1);" << endl;
  file << "plot(t,s);" << endl;
  file << "grid on;" << endl;
  file << "xlabel('time [s]');" << endl;
  file << "ylabel('position [m]');" << endl << endl;
  
  file << "subplot(2,2,2);" << endl;
  file << "plot(t,v);" << endl;
  file << "grid on;" << endl;
  file << "xlabel('time [s]');" << endl;
  file << "ylabel('velocity [m/s]');" << endl << endl;
  
  file << "subplot(2,2,3);" << endl;
  file << "plot(t,m);" << endl;
  file << "grid on;" << endl;
  file << "xlabel('time [s]');" << endl;
  file << "ylabel('mass [kg]');" << endl << endl;
  
  file << "subplot(2,2,4);" << endl;
  file << "plot(t,u);" << endl;
  file << "grid on;" << endl;
  file << "xlabel('time [s]');" << endl;
  file << "ylabel('Thrust [kg m/s^2]');" << endl << endl;
  
  file.close();
  cout << "Results saved to \"" << filename << "\"" << endl;

  }
#endif  
  
}

int main(){
    
  cout << "program started" << endl;
    
  
  // Automatic initialization
  bool manual_init = true;

  // Use the Gauss-Newton method
  bool gauss_newton = false;
  
  // QP-solver
  QPSolverCreator qp_solver_creator = Interfaces::QPOasesSolver::creator;
  Dictionary qp_solver_options;
  qp_solver_options["printLevel"] = "none";

  // Dimensions
  int nu = 20;  // Number of control segments
  int nj = 100; // Number of integration steps per control segment

  // Control
  SXMatrix u = ssym("u",nu); // control

  // Initial conditions
  DMatrix x_0 = DMatrix::zeros(3,1);
  x_0[0] = 0; // initial position
  x_0[1] = 0; // initial speed
  x_0[2] = 1; // initial mass
  
  SX dt = 10.0/(nj*nu); // time step
  SX alpha = 0.05; // friction
  SX beta = 0.1; // fuel consumption rate

  // Trajectory
  SXMatrix s_traj, v_traj, m_traj;

  // Lifted variables
  SXMatrix L;
    
  // For all the shooting intervals
  SXMatrix x = x_0;
  SXMatrix ode_rhs(x.sparsity(),0);
  for(int k=0; k<nu; ++k){
    // Get control
    SX u_k = u[k].at(0);
    
    // Integrate over the interval with Euler forward
    for(int j=0; j<nj; ++j){
      
      // Get the state
      SX s = x.at(0);
      SX v = x.at(1);
      SX m = x.at(2);
      
      // ODE right hand side
      ode_rhs[0] = v;
      ode_rhs[1] = (u_k - alpha * v*v)/m;
      ode_rhs[2] = -beta*u_k*u_k;
      
      // Take a step
      x += dt*ode_rhs;
    }
    
    // Save trajectory
    s_traj.append(x[0]);
    v_traj.append(x[1]);
    m_traj.append(x[2]);

    // Lift x
    L.append(x);
  }
  
  // Objective function
//   SXMatrix f = inner_prod(u,u);
  SXMatrix f = -x[2];
  
  // Terminal constraints
  SXMatrix g;
  g.append(x[0]);
  g.append(x[1]);
  g.append(v_traj);
  
  // Create the NLP
  SXFunction ffcn(u,f); // objective function
  SXFunction gfcn(u,g); // constraint
  
  // Allocate an NLP solver
  IpoptSolver ipopt(ffcn,gfcn);

  // Set options
  ipopt.setOption("generate_hessian",true);
  ipopt.setOption("tol",1e-10);

  // initialize the solver
  ipopt.init();

  // Initial guess for the intermediate variables
  DMatrix x_init = repmat(x_0,nu,1);

  // Bounds on u and initial condition
  vector<double> u_min(nu), u_max(nu), u_init(nu);
  for(int i=0; i<nu; ++i){
    u_min[i] = -10;
    u_max[i] =  10;
    u_init[i] = 0.4;
  }
  ipopt.setInput(u_min,NLP_LBX);
  ipopt.setInput(u_max,NLP_UBX);
  ipopt.setInput(u_init,NLP_X_INIT);
  
  // Bounds on g
  vector<double> g_min(2), g_max(2);
  g_min[0] = g_max[0] = 10;
  g_min[1] = g_max[1] =  0;
  g_min.resize(2+nu, -numeric_limits<double>::infinity());
  g_max.resize(2+nu, 1.1);
  
  ipopt.setInput(g_min,NLP_LBG);
  ipopt.setInput(g_max,NLP_UBG);

  // Solve the problem
  ipopt.solve();

  // Print the optimal solution
  cout << "I: optimal cost:    " << ipopt.output(NLP_COST).toScalar() << endl;
  cout << "I: optimal control: " << ipopt.output(NLP_X_OPT) << endl;
  cout << "I: multipliers (u): " << ipopt.output(NLP_LAMBDA_X) << endl;
  cout << "I: multipliers (g): " << ipopt.output(NLP_LAMBDA_G) << endl;
  
  // Lifting function
  SXFunction ifcn(u,L);

  // Solve problem
  liftedNewton(ffcn,gfcn,ifcn,manual_init,gauss_newton,qp_solver_creator,qp_solver_options,u_min,u_max,u_init,g_min,g_max,x_init);
  
  return 0;
}




