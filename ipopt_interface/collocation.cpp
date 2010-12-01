/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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

#include "collocation.hpp"
#include "casadi/stl_vector_tools.hpp"
#include "casadi/expression_tools.hpp"

using namespace std;
namespace IpoptInterface{

Collocation::Collocation(){
  K = 3;
  N = 15;
  cp = RADAU;
}


void Collocation::optimize(){
  // The right hand side of the ODE
  vector<SXMatrix> y(3);
  y[0] = t;
  y[1] = x;
  y[2] = u;

  ffcn = SXFunction(y,xdot);

  // Objective function (meyer term)
  mfcn = SXFunction(y,m);
  
  // Nonlinear constraint function
  cfcn = SXFunction(y,c);

  // Step size
  double h = tf/N;

  // Coefficients of the collocation equation
  vector<vector<double> > C;

  // Coefficients of the continuity equation
  vector<double> D;

  // Get the coefficients
  get_collocation_coeff(K,C,D,cp);
  
  // All variables with bounds and initial guess
  SXMatrix vars;
  vector< double > vars_lb;
  vector< double > vars_ub;
  vector< double > vars_sol;
  
  // Collocated times
  T.resize(N);
  for(int i=0; i<N; ++i){
    T[i].resize(K+1);
    for(int j=0; j<=K; ++j){
      T[i][j] = h*(i + collocation_points[cp][K][j]);
    }
  }
  
  // Collocated states
  vector< vector< SXMatrix > > X;
  collocate(x,X,N,K);

  // State at end time
  SXMatrix XF;
  collocate_final(x,XF);
  
  // Collocated control (piecewice constant)
  vector< SXMatrix > U;  
  collocate(u,U,N);
    
  // Loop over the finite elements
  for(int i=0; i<N; ++i){
    // collocated controls 
    for(int r=0; r<u.size(); ++r){
      // Add to list of NLP variables
      vars << U[i][r];
      vars_lb.push_back(u_lb[r]);
      vars_ub.push_back(u_ub[r]);
      vars_sol.push_back(u_init[r]);
    }

    // Collocated states
    for(int j=0; j<=K; ++j){
      for(int r=0; r<x.size(); ++r){
        // Add to list of NLP variables
        vars << X[i][j][r];
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
  vars << XF;
  for(int r=0; r<x.size(); ++r){
    vars_sol.push_back(x_init[r]);
    vars_lb.push_back(x_lb[r]);
    vars_ub.push_back(x_ub[r]);
  }
  
  // Constraint function for the NLP
  SXMatrix g;
  vector<double> lbg,ubg;
  for(int i=0; i<N; ++i){
    for(int k=1; k<=K; ++k){
      // augmented state vector
      vector<SXMatrix> y_ik(3);
      y_ik[0] = T[i][k];
      y_ik[1] = X[i][k];
      y_ik[2] = U[i];

      // Add collocation equations to NLP
      SXMatrix rhs(x.size());
      for(int j=0; j<=K; ++j)
        rhs += X[i][j]*C[j][k];
      g << h*ffcn(y_ik) - rhs;
      lbg.insert(lbg.end(),x.size(),0); // equality constraints
      ubg.insert(ubg.end(),x.size(),0); // equality constraints
      
      // Add nonlinear constraints
      g << cfcn(y_ik);
      lbg.insert(lbg.end(),1,0);
      ubg.insert(ubg.end(),1,numeric_limits<double>::infinity());
      
    }

   // Add continuity equation to NLP
   SXMatrix rhs(x.size());
   for(int j=0; j<=K; ++j)
     rhs += D[j]*X[i][j];

   if(i<N-1)
     g << X[i+1][0] - rhs;
   else
     g << XF - rhs;
   lbg.insert(lbg.end(),x.size(),0); // equality constraints
   ubg.insert(ubg.end(),x.size(),0); // equality constraints
  }
  SXFunction gfcn_nlp(vars,g);
  
  // Objective function of the NLP
  vector<SXMatrix> y_f(3);
  y_f[0] = T.back().back();
  y_f[1] = XF;
  y_f[2] = U.back();
  SXMatrix f = mfcn(y_f);
  SXFunction ffcn_nlp(vars, f);

  // Hessian of the Lagrangian (null by default, i.e. BFGS)
  SXFunction hfcn_nlp;
  if(0){
    // Lagrange multipliers
    SXMatrix lambda("lambda",g.size());

    // Objective function scaling
    SXMatrix sigma("sigma");
  
    // Hessian of the Lagrangian
    SXMatrix H = hessian(sigma*f + trans(lambda)*g, vars);
    vector<SXMatrix> hfcn_input(3);
    hfcn_input[0] = vars;
    hfcn_input[1] = lambda;
    hfcn_input[2] = sigma;
    hfcn_nlp = SXFunction(hfcn_input,H);
  }
  
  // the NLP
  NLP nlp(ffcn_nlp,gfcn_nlp,hfcn_nlp);

  // Bounds on x and initial guess
  copy(vars_lb.begin(),vars_lb.end(),nlp->lbx.begin());
  copy(vars_ub.begin(),vars_ub.end(),nlp->ubx.begin());
  copy(vars_sol.begin(),vars_sol.end(),nlp->x.begin());
  
  // Bounds on the constraints
  copy(lbg.begin(),lbg.end(),nlp->lbg.begin());
  copy(ubg.begin(),ubg.end(),nlp->ubg.begin());
  
  // Allocate an NLP solver
  IpoptSolver solver(nlp);

  // Set options
  solver.setOption("abstol",1e-6);

  // initialize the solver
  solver.init();
  
  // Solve the problem
  solver.solve();

  // Print the optimal cost
  cout << "optimal cost: " << nlp->cost << endl;

  // Save optimal solution to disk
  x_opt = nlp->x;
}



} // namespace IpoptInterface

