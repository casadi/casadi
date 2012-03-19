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

  // Initial guess
  double x0 = 0.08;

  // Bounds on the state
  double x_min = -1;
  double x_max =  1;
  double xf_min = 0;
  double xf_max = 0;

  // Control
  double u_min_d = -1;
  double u_max_d =  1;
  
  // End time
  double T = 3;

  // Control discretization
  int nk = 30;

  // Time step
  double dT = T/nk;

  // Discretized control
  SXMatrix u = ssym("u",nk);

  // Initial guess for u
  DMatrix u_guess = DMatrix::zeros(nk);
  DMatrix u_min = u_min_d*DMatrix::ones(nk);
  DMatrix u_max = u_max_d*DMatrix::ones(nk);

  // Lifted variables
  SXMatrix L;

  // Objective terms
  SXMatrix F;

  // NLP functions
  SXFunction F1, F2;
  
  // Get an expression for the state that the final time
  SXMatrix x = x0;
  
  for(int k=0; k<nk; ++k){
    // Get new value for X
    x = x + dT*(x*(x+1)+u[k]);
    
    // Append terms to objective function
    F.append(u[k]);
    F.append(x);
    
    // Lift x
    L.append(x);
  }

  if(gauss_newton){
    // Objective function (GN)
    F1 = SXFunction(u,F);
    
  } else {
    // Objective function (SQP)
    F1 = SXFunction(u,inner_prod(F,F));
  }

  // Constraint function
  F2 = SXFunction(u,x);
  DMatrix g_min = xf_min;
  DMatrix g_max = xf_max;

  // Solve with ipopt
  IpoptSolver nlp_solver(F1,F2);
  nlp_solver.init();
  nlp_solver.setInput(u_guess,NLP_X_INIT);
  nlp_solver.setInput(u_min,NLP_LBX);
  nlp_solver.setInput(u_max,NLP_UBX);
  nlp_solver.setInput(xf_min,NLP_LBG);
  nlp_solver.setInput(xf_max,NLP_UBG);
  nlp_solver.solve();

  // Lifting function
  SXFunction ifcn(u,L);

  // Problem formulation ends
  // Everything below should go into a lifted newton SQP solver class

  // Options
  double tol = 1e-6;     // Stopping tolerance
  int max_iter = 30;  // Maximum number of iterations

  // Extract the free variable and expressions for F and xdef
  u = F1.inputSX();
  SXMatrix f1 = F1.outputSX();
  SXMatrix f2 = F2.outputSX();
  SXMatrix xdef = ifcn.outputSX();

  // Lifted variables
  x = ssym("x",xdef.size());

  // Substitute in the lifted variables x into the expressions for xdef, F1 and F2
  SXMatrixVector ex(2);
  ex[0] = f1;
  ex[1] = f2;
  substituteInPlace(x, xdef, ex, true, false);
  f1 = ex[0];
  f2 = ex[1];

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
    SXMatrix lgrad = gradient(f1 - inner_prod(mug,f2) + inner_prod(xdot,xdef),xu);

    // Gradient of the Lagrangian
    f1 = lgrad(range(u.size1()),0); // + mux // NOTE: What about the mux term?

    // Definition of xdot
    SXMatrix xdotdef = lgrad(range(u.size1(),lgrad.size1()),0);
    
    // Reverse direction of x
    SXMatrix xdot_reversed = xdot;
    for(int k=0; k<xdot.size(); ++k){
      xdot_reversed[k] = xdot[xdot.size()-1-k];
    }
    xdot = xdot_reversed;
    
    SXMatrix xdotdef_reversed = xdotdef;
    for(int k=0; k<xdotdef.size(); ++k){
      xdotdef_reversed[k] = xdotdef[xdotdef.size()-1-k];
    }
    xdotdef = xdotdef_reversed;
    
    // Append to xdef and x
    x.append(xdot);
    xdef.append(xdotdef);
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
  
  // Dimensions
  int nx = x.size1();

  // Variables
  DMatrix u_k = u_guess;
  DMatrix x_k(x.sparsity(),0);
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
    Z.setInput(mug_k,3);
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
    AB.setInput(mug_k,3);
    AB.evaluate();
    DMatrix A_k = AB.output(0);
    DMatrix B1_k = AB.output(1); // NOTE: # mux dissappears (constant term)
    DMatrix B2_k = AB.output(2);
    
    // Get a_k and b_k
    Z.setInput(u_k,0);
    Z.setInput(d_k,1);
    Z.setInput(mux_k,2);
    Z.setInput(mug_k,3);
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
    Z.setFwdSeed(dmug_k,3);
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
    G.setInput(mug_k,3);
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
    cout   << setw(4) <<     k << setw(20) <<   norm_res  << setw(20) <<  norm_step  << setw(20) <<  norm_viol  << endl;
    
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

  
  
  
  return 0;
  
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

  return 0;
  
  }
}





