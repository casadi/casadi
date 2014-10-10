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
#include <ctime>
#include "core/std_vector_tools.hpp"
#include "interfaces/ipopt/ipopt_solver.hpp"
#include "core/sx/sx_tools.hpp"
#include "core/function/sx_function.hpp"

/** Excercise 1, chapter 10 from Larry Biegler's book */

using namespace std;
using namespace casadi;

int main(){
  cout << "program started" << endl;

  std::ofstream resfile;
  resfile.open ("results_biegler_10_1.txt");

  // Test with different number of elements
  for(int N=1; N<=10; ++N){
  
  // Degree of interpolating polynomial
  int K = 2;
  
  // Legrandre roots
  vector<double> tau_root(K+1);
  tau_root[0] = 0.;
  tau_root[1] = 0.211325;
  tau_root[2] = 0.788675;

  // Radau roots (K=3)
  /*  tau_root[0] = 0;
  tau_root[1] = 0.155051;
  tau_root[2] = 0.644949;
  tau_root[3] = 1;*/

  
  // Time
  SX t("t");
  
  // Differential equation
  SX z("z");
  SXFunction F(z,z*z - 2*z + 1);
  F.setOption("name","dz/dt");
  F.init();
  cout << F << endl;
  
  double z0 = -3;
  
  // Analytic solution
  SXFunction z_analytic(t, (4*t-3)/(3*t+1));
  z_analytic.setOption("name","analytic solution");
  z_analytic.init();
  cout << z_analytic << endl;
  
  // Collocation point
  SX tau("tau");

  // Step size
  double h = 1.0/N;
  
  // Lagrange polynomials
  vector<SXFunction> l(K+1);
  for(int j=0; j<=K; ++j){
    SX L = 1;
    for(int k=0; k<=K; ++k)
      if(k != j)
        L *= (tau-tau_root[k])/(tau_root[j]-tau_root[k]);
  
    l[j] = SXFunction(tau,L);
    stringstream ss;
    ss << "l(" << j << ")";
    l[j].setOption("name",ss.str());    
    l[j].init();
    cout << l[j] << endl;
  }
  
  // Get the coefficients of the continuity equation
  vector<double> D(K+1);
  for(int j=0; j<=K; ++j){
    l[j].setInput(1.0);
    l[j].evaluate();
    l[j].getOutput(D[j]);
  }
  cout << "D = " << D << endl;

  // Get the coefficients of the collocation equation
  vector<vector<double> > C(K+1);
  for(int j=0; j<=K; ++j){
    C[j].resize(K+1);
    for(int k=0; k<=K; ++k){
      l[j].setInput(tau_root[k]);
      l[j].setFwdSeed(1.0);
      l[j].evaluate(1,0);
      l[j].getFwdSens(C[j][k]);
    }
  }
  cout << "C = " << C << endl;
  
  // Collocated states
  SX Z = ssym("Z",N,K+1);
  
  // State at final time
// SX ZF("ZF");
  
  // All variables
  SX x;
  x << vec(trans(Z));
  // x << vec(ZF);  
  cout << "x = " << x << endl;
  
  // Construct the "NLP"
  SX g;
  for(int i=0; i<N; ++i){
    for(int k=1; k<=K; ++k){
      
      // Add collocation equations to NLP
      SX rhs = 0;
      for(int j=0; j<=K; ++j)
        rhs += Z(i,j)*C[j][k];
      g << (h*F.eval(SX(Z(i,k))) - rhs);
    }
    
   // Add continuity equation to NLP
   SX rhs = 0;
   for(int j=0; j<=K; ++j)
     rhs += D[j]*Z(i,j);

   if(i<N-1)
     g << (SX(Z(i+1,0)) - rhs);
/*   else
    g << (ZF - rhs);*/
         
  }
  cout << "g = " << g << endl;
    
  SXFunction gfcn(x,g);

  // Dummy objective function
  SXFunction obj(x, Z(0,0)*Z(0,0));
  
  // ----
  // SOLVE THE NLP
  // ----
  
  // Allocate an NLP solver
  IpoptSolver solver(obj,gfcn);

  // Set options
  solver.setOption("tol",1e-10);
  solver.setOption("hessian_approximation","limited-memory");
//   pass_nonlinear_variables

  // initialize the solver
  solver.init();

  // Initial condition
  vector<double> xinit(x.numel(),0);
  solver.setInput(xinit,"x0");

  // Bounds on x
  vector<double> lbx(x.numel(),-100);
  vector<double> ubx(x.numel(), 100);
  lbx[0] = ubx[0] = z0;
  solver.setInput(lbx,"lbx");
  solver.setInput(ubx,"ubx");
  
  // Bounds on the constraints
  vector<double> lubg(g.numel(),0);
  solver.setInput(lubg,"lbg");
  solver.setInput(lubg,"ubg");
  
  // Solve the problem
  solver.solve();
  
  // Print the time points
  vector<double> t_opt(N*(K+1)+1);
  for(int i=0; i<N; ++i)
    for(int j=0; j<=K; ++j)
      t_opt[j + (K+1)*i] = h*(i + tau_root[j]);
  t_opt.back() = 1;
  
  cout << "time points: " << t_opt << endl;
  resfile << t_opt << endl;
  
  // Print the optimal cost
  cout << "optimal cost: " << solver.output(NLP_SOLVER_F) << endl;

  // Print the optimal solution
  vector<double> xopt(x.numel());
  solver.getOutput(xopt,"x");
  cout << "optimal solution: " << xopt << endl;
  resfile << xopt << endl;
  
  }
 
 resfile.close();
  
  return 0;
}
