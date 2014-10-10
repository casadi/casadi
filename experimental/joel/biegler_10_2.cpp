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


#include <iostream>
#include "core/std_vector_tools.hpp"
#include "interfaces/ipopt/ipopt_solver.hpp"
#include "interfaces/ipopt/ipopt_internal.hpp"
#include "core/sx/sx_tools.hpp"
#include "core/function/sx_function.hpp"

/** Excercise 2, chapter 10 from Larry Biegler's book */

using namespace casadi;
using namespace std;

// Legendre collocation points
  static double legendre_points1[] = {0,0.500000};
  static double legendre_points2[] = {0,0.211325,0.788675};
  static double legendre_points3[] = {0,0.112702,0.500000,0.887298};
  static double legendre_points4[] = {0,0.069432,0.330009,0.669991,0.930568};
  static double legendre_points5[] = {0,0.046910,0.230765,0.500000,0.769235,0.953090};
  static double* legendre_points[] = {0,legendre_points1,legendre_points2,legendre_points3,legendre_points4,legendre_points5};

  // Radau collocation points
  static double radau_points1[] = {0,1.000000};
  static double radau_points2[] = {0,0.333333,1.000000};
  static double radau_points3[] = {0,0.155051,0.644949,1.000000};
  static double radau_points4[] = {0,0.088588,0.409467,0.787659,1.000000};
  static double radau_points5[] = {0,0.057104,0.276843,0.583590,0.860240,1.000000};
  static double* radau_points[] = {0,radau_points1,radau_points2,radau_points3,radau_points4,radau_points5};

  // Type of collocation points
  enum CollocationPoints{LEGENDRE,RADAU};
  static double** collocation_points[] = {legendre_points,radau_points};

void get_coeff(int K, vector<vector<double> >& C, vector<double>& D, CollocationPoints cp){
  D.resize(K+1);
  C.resize(K+1);
  
  // Collocation point
  SX tau("tau");
  
  // Roots
  double *tau_root = collocation_points[cp][K];

  for(int j=0; j<=K; ++j){
    // Lagrange polynomials
    SX L = 1;
    for(int k=0; k<=K; ++k)
      if(k != j)
        L *= (tau-tau_root[k])/(tau_root[j]-tau_root[k]);
  
    // Lagrange polynomial function
    SXFunction lfcn(tau,L);
    lfcn.init();
    
    // Get the coefficients of the continuity equation
    lfcn.setInput(1.0);
    lfcn.evaluate();
    lfcn.getOutput(D[j]);

    // Get the coefficients of the collocation equation
    C[j].resize(K+1);
    for(int k=0; k<=K; ++k){
      lfcn.setInput(tau_root[k]);
      lfcn.setFwdSeed(1.0);
      lfcn.evaluate(1,0);
      lfcn.getFwdSens(C[j][k]);
    }
  }
}

int main(){
  cout << "program started" << endl;

    // Time 
  SX t("t");
  
  // Differential states
  SX Ls("Ls");    // mean crystal size
  SX Nc("Nc");    // number of nuclei per liter of solvent
  SX L("L");      // total length of crystals per liter of solvent
  SX Ac("Ac");    // total surface area of the crystals per liter of solvent
  SX Vc("Vc");    // total volume of the crysals per liter of solvent
  SX Mc("Mc");    // total mass of the crystals
  SX Cc("Cc");    // solute concentration
  SX Tc("Tc");    // cystillizer temperature

  // State vector
  vector<SX> x;
  x.push_back(Ls);
  x.push_back(Nc);
  x.push_back(L);
  x.push_back(Ac);
  x.push_back(Vc);
  x.push_back(Mc);
  x.push_back(Cc);
  x.push_back(Tc);
  
  // Bounds on the states
  vector<double> x_lb, x_ub;
  x_lb = vector<double>(x.size(),-numeric_limits<double>::infinity());
  x_ub = vector<double>(x.size(),numeric_limits<double>::infinity());

  // Initial values
  double Ls_init = 0.0005;
  double Nc_init = 0;
  double L_init = 0;
  double Ac_init = 0;
  double Vc_init = 0;
  double Mc_init = 2.0;
  double Cc_init = 5.4;
  double Tc_init = 75;
  double x_init_array[] = {Ls_init,Nc_init,L_init,Ac_init,Vc_init,Mc_init,Cc_init,Tc_init};

  // Initial values of the states
  vector<double> x_init(x_init_array, x_init_array+8);

  // Control
  SX Tj("Tj"); // jacket temperature
  double Tj_lb = 10, Tj_ub = 60, Tj_init = 30;

  // Control vector
  vector<SX> u;  u.push_back(Tj);

  // Bounds on the control
  vector<double> u_lb, u_ub;
  u_lb.push_back(Tj_lb);
  u_ub.push_back(Tj_ub);

  // Initial values for the controls
  vector<double> u_init;
  u_init.push_back(Tj_init);

  // Constants
  double Vs = 300; // volume of the solvent
  double W = 2025; // the total mass in the crysallizer
  double a[] = {-66.4309, 2.8604, -0.022579, 6.7117e-5};
  double b[] = {16.08852, -2.708263, 0.0670694, -3.5685e-4};
  double Kg = 0.00418;
  double Bn = 385;
  double Cp = 0.4;
  double Kc = 35;
  double Ke = 377;
  double eta1 = 1.1;
  double eta2 = 5.72;
  double Ls0 = 5e-4; // initial crystal size
  double L0 = 5e-5; // nucleate crystal size
  double Ws0 = 2; // weight of the seed crystals
  double rho = 1.58; // specific gravity of crystals
  double alpha = 0.2; // shape factor for area of crystals
  double beta = 1.2; // shape factor for volume of crystals

  // Time horizon
  double tf = 16;

  // Dependent variables
  SX C_bar = 100*Cc/(1.35+Cc);
  SX Tequ = a[0] + a[1]*C_bar + a[2]*C_bar*C_bar + a[3]*C_bar*C_bar*C_bar; // equilibrium temperature
  SX Ta = b[0] + b[1]*C_bar + b[2]*C_bar*C_bar + b[3]*C_bar*C_bar*C_bar; // lower bound of Tj
  SX DeltaT; // degree of supercooling
  // DeltaT = fmax(0,Tequ-Tc);   // Original formulation
//  DeltaT = fmax(1e-8,Tequ-Tc); // "epsilon" to avoid divide by zero
  DeltaT = log(exp(1e-8)+exp(Tequ-Tc));   // Log sum exp

  // Differential equations
  SX Ls_dot = Kg*sqrt(Ls)*pow(DeltaT,eta1);
  SX Nc_dot = Bn*pow(DeltaT,eta2);
  SX L_dot = Nc*Ls_dot + L0*Nc_dot;
  SX Ac_dot = 2*alpha*Nc*Ls_dot + L0*L0*Nc_dot;
  SX Vc_dot = 3*beta*Ac*Ls_dot + L0*L0*L0*Nc_dot;
  SX Mc_dot = 3*(Ws0/(Ls0*Ls0*Ls0))*Ls*Ls*Ls_dot + rho*Vs*Vc_dot;
  SX Cc_dot = -1/Vs*Mc_dot;
  SX Tc_dot = (Kc*Mc_dot - Ke*(Tc-Tj))/(W*Cp);
  
  SX xdot_array[] = {Ls_dot, Nc_dot, L_dot, Ac_dot, Vc_dot, Mc_dot, Cc_dot, Tc_dot};
  vector<SX> xdot(xdot_array,xdot_array+8);

  // Right hand side of the ODE
  vector<vector<SX> > y(3);
  y[0] = vector<SX>(1,t);
  y[1] = x;
  y[2] = u;
  SXFunction ffcn(y,xdot);
  ffcn.init();
  
  // Objective function (meyer term)
  SXFunction mfcn(y,-Ls);
  mfcn.init();

  // Nonlinear constraint function
  SXFunction cfcn(y,Tj-Ta);
  cfcn.init();
  
  // Degree of interpolating polynomial
  int K = 3;
  
  // Number of finite elements
  int N = 15;

  // Radau collocation points
  CollocationPoints cp = RADAU;

  // Size of the finite elements
  double h = tf/N;

  // Coefficients of the collocation equation
  vector<vector<double> > C;

  // Coefficients of the continuity equation
  vector<double> D;

  // Get the coefficients
  get_coeff(K,C,D,cp);
                      
  // Collocated times
  vector<vector<double> > T;
  T.resize(N);
  for(int i=0; i<N; ++i){
    T[i].resize(K+1);
    for(int j=0; j<=K; ++j){
      T[i][j] = h*(i + collocation_points[cp][K][j]);
    }
  }
  
  // Collocated states
  vector< vector< Matrix<SX> > > X = ssym("X",x.size(),1,K+1,N);
  
  // Collocated control (piecewice constant)
  vector< Matrix<SX> > U = ssym("U", u.size(),1,N);
  
  // State at end time
  Matrix<SX> XF = ssym("XF",x.size());
    
  // All variables with bounds and initial guess
  Matrix<SX> vars;
  vector<double> vars_lb;
  vector<double> vars_ub;
  vector<double> vars_sol;

  // Loop over the finite elements
  for(int i=0; i<N; ++i){
    // collocated controls 
    vars.append( U[i] );
    for(int r=0; r<u.size(); ++r){
      // Add to list of NLP variables
      vars_lb.push_back(u_lb[r]);
      vars_ub.push_back(u_ub[r]);
      vars_sol.push_back(u_init[r]);
    }

    // Collocated states
    for(int j=0; j<=K; ++j){
      vars.append(X[i][j]);
      for(int r=0; r<x.size(); ++r){
        // Add to list of NLP variables
        vars_sol.push_back(x_init[r]);
        if(i==0 && j==0){
          // Initial constraints
          vars_lb.push_back(x_init[r]);
          vars_ub.push_back(x_init[r]);
        } else {
          // Variable bounds
          vars_lb.push_back(x_lb[r]);
          vars_ub.push_back(x_ub[r]);
        }
      }

    }
  }
  
  // Add states at end time
  vars.append(XF);
  for(int r=0; r<x.size(); ++r){
    vars_sol.push_back(x_init[r]);
    vars_lb.push_back(x_lb[r]);
    vars_ub.push_back(x_ub[r]);
  }
  
  // Constraint function for the NLP
  Matrix<SX> g;
  
  vector<double> lbg,ubg;
  for(int i=0; i<N; ++i){
    for(int k=1; k<=K; ++k){
      // augmented state vector
      vector<Matrix<SX> > y_ik(3);
      y_ik[0] = T[i][k];
      y_ik[1] = X[i][k];
      y_ik[2] = U[i];

      // Add collocation equations to NLP
      Matrix<SX> temp = ffcn.eval(y_ik)[0];
      temp *= h;
      for(int j=0; j<=K; ++j)
        temp -= X[i][j]*C[j][k];
      
      g.append(temp);
      lbg.insert(lbg.end(),x.size(),0); // equality constraints
      ubg.insert(ubg.end(),x.size(),0); // equality constraints
      
      // Add nonlinear constraints
      temp = cfcn.eval(y_ik)[0];
      g.append(temp);
      lbg.insert(lbg.end(),1,0);
      ubg.insert(ubg.end(),1,numeric_limits<double>::infinity());
    }

   // Add continuity equation to NLP
   Matrix<SX> temp = i<N-1 ? X[i+1][0] : XF;
   for(int j=0; j<=K; ++j)
     temp -= D[j]*X[i][j];

   g.append(temp);
   lbg.insert(lbg.end(),x.size(),0); // equality constraints
   ubg.insert(ubg.end(),x.size(),0); // equality constraints
  }
  
  // Nonlinear constraint function
  SXFunction gfcn_nlp(vars,g);
  
  // Objective function of the NLP
  vector<Matrix<SX> > y_f(3);
  y_f[0] = T.back().back();
  y_f[1] = XF;
  y_f[2] = U.back();
  Matrix<SX> f = mfcn.eval(y_f)[0];
  SXFunction ffcn_nlp(vars, f);
  
  // Hessian of the Lagrangian:
  // Lagrange multipliers
  Matrix<SX> lambda = ssym("lambda",g.size());

  // Objective function scaling
  Matrix<SX> sigma = ssym("sigma");

  // Lagrangian function
  vector< Matrix<SX> > lfcn_input(3);
  lfcn_input[0] = vars;
  lfcn_input[1] = lambda;
  lfcn_input[2] = sigma;
  SXFunction lfcn(lfcn_input, sigma*f + inner_prod(lambda,g));
  lfcn.init();
  
  // Hessian of the Lagrangian
  SXFunction HL = lfcn.hessian();

    // ----
  // SOLVE THE NLP
  // ----
  
  // Allocate an NLP solver
  IpoptSolver solver(ffcn_nlp,gfcn_nlp,HL);

  // Set options
  solver.setOption("tol",1e-6);
/*  solver.setOption("hessian_approximation","limited-memory");
  solver.setOption("pass_nonlinear_variables",false);*/
//  solver.setOption("derivative_test","first-order");

  // initialize the solver
  solver.init();
  
  // Initial condition
  solver.setInput(vars_sol,"x0");

  // Bounds on x
  solver.setInput(vars_lb,"lbx");
  solver.setInput(vars_ub,"ubx");

  // Bounds on g
  solver.setInput(lbg,"lbg");
  solver.setInput(ubg,"ubg");

  // Solve the problem
  solver.solve();

  // Print the optimal cost
  double cost;
  solver.getOutput(cost,"f");
  cout << "optimal cost: " << cost << endl;

  // Get the solution
  solver.getOutput(vars_sol,"x");

  // ----
  // SAVE SOLUTION TO DISK
  // ----

  // Get the optimal solution
  vector<double> Tj_opt(N);
  vector<double> Ls_opt(N*(K+1));
  vector<double> Nc_opt(N*(K+1));
  vector<double> L_opt(N*(K+1));
  vector<double> Ac_opt(N*(K+1));
  vector<double> Vc_opt(N*(K+1));
  vector<double> Mc_opt(N*(K+1));
  vector<double> Cc_opt(N*(K+1));
  vector<double> Tc_opt(N*(K+1));
  vector<double> t_opt(N*(K+1));
  int ind = 0; // index of nlp->x
  for(int i=0; i<N; ++i){
    Tj_opt[i] = vars_sol[ind++];
    for(int j=0; j<=K; ++j){
      int ij = (K+1)*i+j;
      Ls_opt[ij] = vars_sol[ind++];
      Nc_opt[ij] = vars_sol[ind++];
      L_opt[ij] = vars_sol[ind++];
      Ac_opt[ij] = vars_sol[ind++];
      Vc_opt[ij] = vars_sol[ind++];
      Mc_opt[ij] = vars_sol[ind++];
      Cc_opt[ij] = vars_sol[ind++];
      Tc_opt[ij] = vars_sol[ind++];
      t_opt[ij] = T[i][j];
    }
  }
  
  std::ofstream resfile;
  resfile.open ("results_biegler_10_2.txt");

  // Save optimal solution to disk
  resfile << "Ls_opt " << Ls_opt << endl;
  resfile << "Nc_opt " << Nc_opt << endl;
  resfile << "L_opt " << L_opt << endl;
  resfile << "Ac_opt " << Ac_opt << endl;
  resfile << "Vc_opt " << Vc_opt << endl;
  resfile << "Mc_opt " << Mc_opt << endl;
  resfile << "Cc_opt " << Cc_opt << endl;
  resfile << "Tc_opt " << Tc_opt << endl;
  resfile << "Tj_opt " << Tj_opt << endl;
  resfile << "t_opt " << t_opt << endl;
  
  resfile.close();

  return 0;
}
