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
#include <casadi/stl_vector_tools.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <casadi/mx/mx_tools.hpp>
#include <casadi/fx/mx_function.hpp>
#include <casadi/expression_tools.hpp>
#include "interfaces/sundials/cvodes_integrator.hpp"
#include "interfaces/sundials/idas_integrator.hpp"

using namespace CasADi;
using namespace CasADi::Sundials;
using namespace std;

// CONSTRUCT THE INTEGRATOR
FX create_integrator_euler(){
  SX u("u"); // control for one segment

  // Initial position
  SX s0("s0"); // initial position
  SX v0("v0"); // initial speed
  SX m0("m0"); // initial mass

  SX t0("t0"); // initial time
  SX tf("tf"); // final time

  int nj = 1000; // Number of integration steps per control segment
  SX dt = (tf-t0)/nj; // time step
  SX alpha = 0.05; // friction
  SX beta = 0.1; // fuel consumption rate

  // Integrate over the interval with Euler forward
  SX s = s0, v = v0, m = m0;

  SX dm = -dt*beta*u*u;
  for(int j=0; j<nj; ++j){
    s += dt*v;
    v += dt * (u-alpha*v*v)/m;
    m += dm;
  }

  // State vector
  SXMatrix x, x0;
  x << s << v << m;
  x0 << s0 << v0 << m0;

  // Integrator
  vector<SXMatrix> input(INTEGRATOR_NUM_IN);
  input[INTEGRATOR_T0] = t0;
  input[INTEGRATOR_TF] = tf;
  input[INTEGRATOR_X0] = x0;
  input[INTEGRATOR_P] = u;

  vector<SX> xp0(x0.size()); make_symbolic(xp0.begin(),xp0.end(),"xp0");
  input[INTEGRATOR_XP0] = vertcat(xp0);

  SXMatrix output = x;
  SXFunction integrator(input,output);
  integrator.setOption("ad_order",1);
  integrator.init();

//  integrator->generateCode("rocket.c");

  return integrator;
}

FX create_integrator_sundials(bool explicit_integrator){
  // Time 
  SX t("t");

  // Differential states
  SX s("s"), v("v"), m("m");
  vector<SX> y(3); 
  y[0] = s;
  y[1] = v;
  y[2] = m;

  // Control
  SX u("u");

  SX alpha = 0.05; // friction
  SX beta = 0.1; // fuel consumption rate
  
  // Differential equation
  vector<SX> rhs(3);
  rhs[0] = v;              // sdot
  rhs[1] = (u-alpha*v*v)/m; // vdot
  rhs[2] = -beta*u*u;      // mdot

  // Initial conditions
  vector<double> y0(3);
  y0[0] = 0;
  y0[1] = 0;
  y0[2] = 1;
  
  // Integrator
  Integrator integrator;

  // Create an integrator
  if(explicit_integrator){
    
    // Input of the ode rhs
    vector<vector<SX> > ffcn_in(ODE_NUM_IN);
    ffcn_in[ODE_T].push_back(t);
    ffcn_in[ODE_Y] = y;
    ffcn_in[ODE_P].push_back(u);
  
    // ODE right hand side
    SXFunction ffcn(ffcn_in,rhs);
    ffcn.setOption("name","ODE right hand side");
    ffcn.setOption("ad_order",1);

    // Explicit integrator (CVODES)
    integrator = Sundials::CVodesIntegrator(ffcn);
    // integrator.setOption("exact_jacobian",true);
    // integrator.setOption("linear_multistep_method","bdf"); // adams or bdf
    // integrator.setOption("nonlinear_solver_iteration","newton"); // newton or functional
  } else {
       // Implicit integrator (IDAS)

      // State derivative
      SX sdot("sdot"), vdot("vdot"), mdot("mdot");
      vector<SX> ydot(3); 
      ydot[0] = sdot;
      ydot[1] = vdot;
      ydot[2] = mdot;

      // Input of the dae residual
      vector<vector<SX> > ffcn_in(DAE_NUM_IN);
      ffcn_in[DAE_T].push_back(t);
      ffcn_in[DAE_Y] = y;
      ffcn_in[DAE_YDOT] = ydot;
      ffcn_in[DAE_P].push_back(u);

      // DAE residual function
      vector<SX> res = ydot;
      for(int i=0; i<res.size(); ++i)
        res[i] -= rhs[i];
      SXFunction ffcn(ffcn_in,res);
      ffcn.setOption("name","ODE right hand side");
      ffcn.setOption("ad_order",1);

      // Create an integrator
      integrator = Sundials::IdasIntegrator(ffcn);
      integrator.setOption("calc_ic",false);
  }

  integrator.setOption("ad_order",1);
  integrator.setOption("fsens_err_con",true);
  integrator.setOption("quad_err_con",true);
  integrator.setOption("abstol",1e-6);
  integrator.setOption("reltol",1e-6);
  integrator.setOption("stop_at_end",false);
//  integrator.setOption("fsens_all_at_once",false);
  integrator.setOption("steps_per_checkpoint",100); // BUG: Too low number causes segfaults

  integrator.init();

  return integrator;
}


int main(){
  
  // Time length
  double T = 10.0;

  // Shooting length
  int nu = 20; // Number of control segments
  double DT = T/nu;

  // Initial position
  vector<double> X0(3);
  X0[0] = 0; // initial position
  X0[1] = 0; // initial speed
  X0[2] = 1; // initial mass


  // Create an integrator
//  FX integrator = create_integrator_euler();
//  FX integrator = create_integrator_sundials(true);
   FX integrator = create_integrator_sundials(false);

  // control for all segments
  MX U("U",nu); 

  // Integrate over all intervals
  MX X=X0;
  MX T0 = 0;
  MX TF = DT;
  for(int k=0; k<nu; ++k){
    // Assemble the input
    vector<MX> input(INTEGRATOR_NUM_IN);
    input[INTEGRATOR_T0] = T0; // k*DT
    input[INTEGRATOR_TF] = TF; // (k+1)*DT
    input[INTEGRATOR_X0] = X;
    input[INTEGRATOR_P] = U[k];

    vector<double> xp(X.numel());
    input[INTEGRATOR_XP0] = xp;
    
    // Integrate
    X = integrator(input);
  }

  // Objective function
  MX F = inner_prod(U,U);

  // Terminal constraints
  MX G = vertcat(X[0],X[1]);
  
  // Create the NLP
  MXFunction ffcn(U,F); // objective function
  MXFunction gfcn(U,G); // constraint function

  // Allocate an NLP solver
//  LiftedNewtonSolver solver(ffcn,gfcn);
  IpoptSolver solver(ffcn,gfcn);
  
  // Set options
  solver.setOption("tol",1e-10);
  solver.setOption("hessian_approximation","limited-memory");
  
  // initialize the solver
  solver.init();

  // Bounds on u and initial condition
  vector<double> Umin(nu), Umax(nu), Usol(nu);
  for(int i=0; i<nu; ++i){
    Umin[i] = -10;
    Umax[i] =  10;
    Usol[i] = 0.4;
  }
  solver.setInput(Umin,NLP_LBX);
  solver.setInput(Umax,NLP_UBX);
  solver.setInput(Usol,NLP_X_INIT);

  // Bounds on g
  vector<double> Gmin(2), Gmax(2);
  Gmin[0] = Gmax[0] = 10;
  Gmin[1] = Gmax[1] =  0;
  solver.setInput(Gmin,NLP_LBG);
  solver.setInput(Gmax,NLP_UBG);

  // Solve the problem
  solver.solve();

  // Get the solution
  solver.getOutput(Usol,NLP_X_OPT);
  cout << "optimal solution: " << Usol << endl;

  return 0;

}

