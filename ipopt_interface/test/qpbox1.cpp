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

#include <iostream>
#include <ctime>
#include "casadi/stl_vector_tools.hpp"
#include "ipopt_interface/ipopt_solver.hpp"
#include "casadi/expression_tools.hpp"

using namespace CasADi;
using namespace std;

int main(){
  try{

    // QP example from Matlab

    // Variable
    SXMatrix x("x",400);
  
    // Hessian
    SXMatrix H(400,400);
    for(int i=0; i<400; ++i){
      if(i>0) H(i,i-1) = -2;
      H(i,i) = 5;
      if(i<400-1) H(i,i+1) = -2;
    }
    
    SXMatrix f(400,1);
    f[0] = -2;
    f[399] = -2;

    // Objective function
    SXMatrix F = trans(x)*H*x + trans(f)*x;
    SXFunction ffcn(x,F);
    ffcn.setOption("ad_order",1);
    ffcn.init();

    // Nonlinear constraint 
    SXMatrix G = x[0];
    SXFunction gfcn(x,G);
    gfcn.setOption("ad_order",1);
    gfcn.init();
    
    // Hessian
    SXFunction hfcn(x,H);
    hfcn.setOption("ad_order",1);
    hfcn.init();    
    
    // Create an NLP
    NLP nlp(ffcn,gfcn,hfcn);
    
    // Upper variable bounds
    fill(nlp->lbx.begin(),nlp->lbx.end(),0);
    nlp->lbx[399] = -1000;

    // Lower variable bounds
    fill(nlp->ubx.begin(),nlp->ubx.end(),0.9);
    nlp->ubx[399] = 1000;

    // Upper nonlinear bound
    nlp->lbg[0] = -1; // never active

    // Lower nonlinear bound
    nlp->ubg[0] = 2; // never active
    
    // Set starting value
    fill(nlp->x.begin(),nlp->x.end(),0.5);

    // Allocate an NLP solver
    IpoptSolver solver(nlp);
    
    // Set options
    solver.setOption("exact_hessian",true);
    solver.setOption("abstol",1e-8);

    // initialize the solver
    solver.init();

    // Solve the problem
    solver.solve();
  
    
    
//    cout << F << endl;
    
    
    
    
    
    
    
//    cout << "H = " << H << endl;
    

return 0;

    
    
    
#if 0    
    cout << "program started" << endl;
    
  // Dimensions
  int nu = 1000;  // Number of control segments
  int nj = 1000; // 10000;  // // Number of integration steps per control segment

  // optimization variable
  SXMatrix u("u",nu); // control

  SX s_0 = 0; // initial position
  SX v_0 = 0; // initial speed
  SX m_0 = 1; // initial mass
  
  SX dt = 10.0/(nj*nu); // time step
  SX alpha = 0.05; // friction
  SX beta = 0.1; // fuel consumption rate

  // Integrate over the interval with Euler forward
  SX s = s_0, v = v_0, m = m_0;
  for(int k=0; k<nu; ++k){
    for(int j=0; j<nj; ++j){
      s += dt*v;
      v += dt / m * (u[k] - alpha * v*v);
      m += -dt * beta*u[k]*u[k];
    }
  }

// Objective function
  SXMatrix f = trans(u)*u;

  // Terminal constraints
  SXMatrix g(2,1);
  g[0] = s;
  g[1] = v;

  // Create the NLP
  SXFunction ffcn(u,f); // objective function
  ffcn->setOption("ad_order",1);
  ffcn->init();

  SXFunction gfcn(u,g); // constraint
  gfcn->setOption("ad_order",1);
  gfcn->init();

#if 0
  // Create an exact Hessian of the lagrangian function

  // objective scaling function
  SX sigma("sigma"); 
  
  // Lagrange multipliers
  SXMatrix lambda("lambda",g.size());

  // Lagrangian function
  SXMatrix L = sigma*f + trans(lambda)*g;
  
  // Calculate the Hessian of the Lagrangian symbolically
  SXMatrix HL = L.hessian(u);
  
  // Create a Hessian of the Lagrangian function
  vector<SXMatrix> hfcn_in(3);
  hfcn_in[0] = u;
  hfcn_in[1] = lambda;
  hfcn_in[2] = sigma;  
  SXFunction hfcn(hfcn_in,HL);
  hfcn.init();
  
  // Create the nonlinear program
  NLP nlp(ffcn,gfcn,hfcn);
#else
  NLP nlp(ffcn,gfcn);  
#endif
  
  
  // Bounds on u and initial condition
  vector<double> umin(nu), umax(nu), usol(nu);
  for(int i=0; i<nu; ++i){
    umin[i] = -10;
    umax[i] =  10;
    usol[i] = 0.4;
  }
  copy(umin.begin(),umin.end(),nlp->lbx.begin());
  copy(umax.begin(),umax.end(),nlp->ubx.begin());
  copy(usol.begin(),usol.end(),nlp->x.begin());

  // Upper and lower bounds
  vector<double> gmin(2), gmax(2);
  gmin[0] = gmax[0] = 10;
  gmin[1] = gmax[1] =  0;
  copy(gmin.begin(),gmin.end(),nlp->lbg.begin());
  copy(gmax.begin(),gmax.end(),nlp->ubg.begin());

  // Allocate an NLP solver
  IpoptSolver solver(nlp);

  // Set options
  solver.setOption("exact_hessian",false);
  solver.setOption("abstol",1e-10);

  // initialize the solver
  solver.init();

  // Solve the problem
  solver.solve();

  // Get the optimal cost
//  double fopt = nlp->eval_f();
//  cout << "optimal cost: " << fopt << endl;

  // Get the optimal solution
  copy(nlp->x.begin(),nlp->x.end(),usol.begin());
  cout << "optimal solution: " << usol << endl;
#endif
  
  return 0;
  }
  catch (exception& e){
    cout << "rocket_ipopt failed: " << e.what() << endl;
    return 1;
  }
  catch (const char * str){
    cerr << "rocket_ipopt failed (OLD EXCEPTION TYPE!): " << str << endl;
    return 1;
  }
}