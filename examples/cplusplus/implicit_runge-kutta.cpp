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

/**
   Demonstration on how to construct a fixed-step implicit Runge-Kutta integrator
   @author: Joel Andersson, K.U. Leuven 2013
*/

#include <symbolic/casadi.hpp>
#include <integration/integration_tools.hpp>
#include <nonlinear_programming/newton_implicit_solver.hpp>
#include <interfaces/csparse/csparse.hpp>
#include <interfaces/sundials/cvodes_integrator.hpp>
#include <iomanip>

using namespace CasADi;
using namespace std;

int main(){

  // End time
  double tf = 10.0;

  // Dimensions
  int nx = 3;
  int np = 1;

  // Declare variables
  SXMatrix x  = ssym("x",nx);  // state
  SXMatrix p  = ssym("u",np);  // control

  // ODE right hand side function
  SXMatrix ode;
  ode.append( (1 - x[1]*x[1])*x[0] - x[1] + p );
  ode.append( x[0]                            );
  ode.append( x[0]*x[0] + x[1]*x[1] + p*p     );
  SXFunction f(daeIn("x",x,"p",p),daeOut("ode",ode));
  f.init();

  // Number of finite elements
  int n = 100;

  // Size of the finite elements
  double h = tf/n;

  // Choose collocation points
  vector<double> tau_root = collocationPoints(4,"legendre");

  // Degree of interpolating polynomial
  int d = tau_root.size()-1;

  // Coefficients of the collocation equation
  vector<vector<double> > C(d+1,vector<double>(d+1,0));

  // Coefficients of the continuity equation
  vector<double> D(d+1,0);

  // Dimensionless time inside one control interval
  SXMatrix tau = ssym("tau");
  
  // For all collocation points
  for(int j=0; j<d+1; ++j){
    // Construct Lagrange polynomials to get the polynomial basis at the collocation point
    SXMatrix L = 1;
    for(int r=0; r<d+1; ++r){
      if(r!=j){
        L *= (tau-tau_root[r])/(tau_root[j]-tau_root[r]);
      }
    }
    SXFunction lfcn(tau,L);
    lfcn.init();
  
    // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
    lfcn.setInput(1.0);
    lfcn.evaluate();
    lfcn.getOutput(D[j]);

    // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
    for(int r=0; r<d+1; ++r){
      lfcn.setInput(tau_root[r]);
      lfcn.setFwdSeed(1.0);
      lfcn.evaluate(1,0);
      lfcn.getFwdSens(C[j][r]);
    }
  }

  // Total number of variables for one finite element
  MX X0 = msym("X0",nx);
  MX P  = msym("P",np);
  MX V = msym("V",d*nx);
  
  // Get the state at each collocation point
  vector<MX> X(1,X0);
  for(int r=0; r<d; ++r){
    X.push_back(V[Slice(r*nx,(r+1)*nx)]);
  }
  
  // Get the collocation quations (that define V)
  MX V_eq;
  for(int j=1; j<d+1; ++j){
    // Expression for the state derivative at the collocation point
    MX xp_j = 0;
    for(int r=0; r<d+1; ++r){
      xp_j += C[r][j]*X[r];
    }
    
    // Append collocation equations
    MX f_j = daeOut(f.call(daeIn("x",X[j],"p",P)),"ode").front();
    V_eq.append(h*f_j - xp_j);
  }

  // Root-finding function, implicitly defines V as a function of X0 and P
  vector<MX> vfcn_in;
  vfcn_in.push_back(V);
  vfcn_in.push_back(X0);
  vfcn_in.push_back(P);
  MXFunction vfcn(vfcn_in,V_eq);
  vfcn.init();
  
  // Convert to SXFunction to decrease overhead
  SXFunction vfcn_sx(vfcn);

  // Create a implicit function instance to solve the system of equations
  NewtonImplicitSolver ifcn(vfcn_sx);
  ifcn.setOption("linear_solver",CSparse::creator);
  ifcn.init();
  vector<MX> ifcn_arg;
  ifcn_arg.push_back(X0);
  ifcn_arg.push_back(P);
  V = ifcn.call(ifcn_arg).front();
  X.resize(1);
  for(int r=0; r<d; ++r){
    X.push_back(V[Slice(r*nx,(r+1)*nx)]);
  }

  // Get an expression for the state at the end of the finie element
  MX XF = 0;
  for(int r=0; r<d+1; ++r){
    XF += D[r]*X[r];
  }
  
  // Get the discrete time dynamics
  vector<MX> F_arg;
  F_arg.push_back(X0);
  F_arg.push_back(P);
  MXFunction F(F_arg,XF);
  F.init();

  // Do this iteratively for all finite elements
  for(int i=0; i<n; ++i){
    F_arg[0] = F.call(F_arg).front();
  }

  // Fixed-step integrator
  MXFunction irk_integrator(integratorIn("x0",X0,"p",P),integratorOut("xf",F_arg[0]));
  irk_integrator.setOption("name","irk_integrator");
  irk_integrator.setOption("number_of_fwd_dir",2);
  irk_integrator.init();

  // Create a convensional integrator for reference
  CVodesIntegrator ref_integrator(f);
  ref_integrator.setOption("name","ref_integrator");
  ref_integrator.setOption("number_of_fwd_dir",2);
  ref_integrator.setOption("tf",tf);
  ref_integrator.init();

  // Test values
  double x0_val[] = {0,1,0};
  double p_val = 0.2;

  // Make sure that both integrators give consistent results
  for(int integ=0; integ<2; ++integ){
    FX integrator = integ==0 ? FX(irk_integrator) : FX(ref_integrator);
    cout << "-------" << endl;
    cout << "Testing " << integrator.getOption("name") << endl;
    cout << "-------" << endl;

    // Pass arguments
    integrator.setInput(x0_val,"x0");
    integrator.setInput(p_val,"p");
  
    // Forward sensitivity analysis, first direction: seed p
    integrator.setFwdSeed(0.0,"x0",0);
    integrator.setFwdSeed(1.0,"p",0);
  
    // Forward sensitivity analysis, second direction: seed x0[0]
    double x0_seed[] = {1,0,0};
    integrator.setFwdSeed(x0_seed,"x0",1);
    integrator.setFwdSeed(0.0,"p",1);
  
    // Adjoint sensitivity analysis, seed xf[2]
    double xf_seed[] = {0,0,1};
    integrator.setAdjSeed(xf_seed,"xf");

    // Integrate
    integrator.evaluate(2,1);

    // Get the nondifferentiated results
    cout << setw(15) << "xf = " << integrator.output("xf") << endl;

    // Get the forward sensitivities
    cout << setw(15) << "d(xf)/d(p) = " <<  integrator.fwdSens("xf",0) << endl;
    cout << setw(15) << "d(xf)/d(x0[0]) = " << integrator.fwdSens("xf",1) << endl;

    // Get the adjoint sensitivities
    cout << setw(15) << "d(xf[2])/d(x0) = " << integrator.adjSens("x0") << endl;
    cout << setw(15) << "d(xf[2])/d(p) = " << integrator.adjSens("p") << endl;
  }
  return 0;
}

