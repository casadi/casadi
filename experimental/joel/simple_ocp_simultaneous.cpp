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

/*
 *
 *		Attempt to solve the OCP: (Simultaneous approach; control and state as variables)
 *
 *			min sum (x_i ^2 + u_i ^2)
 *			 x in R^31
 *			 u in R^30
 *
 *			subject to:
 *
 *			x_0 = 0.2
 *
 *			-1 <= x_i <= 1							for all i in {1,...,30}
 *			-1 <= u_i <= 1							for all i in {1,...,30}
 *
 *			x_{i+1} = x_i+0.1*(x_i+x_i ^2 + u_i) 	for all i in {1,...,30}
 *
 */

#include <iostream>
#include <ctime>
#include <core/casadi.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <core/std_vector_tools.hpp>

using namespace casadi;
using namespace std;

int main(){
    
  cout << "program started" << endl;
      
  // Dimensions
  int nu = 30;  	// Number of control segments
  int nx = nu+1;	// State

  // Optimization variable
  vector<SX> v = ssym("v",nu+nx).data(); // control and state

  // Objective function
  SX f = 0; //Could take x_0 into account, but since it is a constant...
  for(int i=0; i<nu; ++i)
    f += v[i]*v[i]+v[i+nu+1]*v[i+nu+1]; 
  
  // Constraints (xEq <= g <= xEq (equality constraints, due to initial condition and state evolution))
  vector<SX> g(nx);
  g[0]=v[nu]-0.2;
  for (int i=1; i<nx; ++i)
	g[i]=v[i+nu]-(v[i+nu-1]+0.1*(v[i+nu-1]*(1+v[i+nu-1])+v[i-1]));

  SXFunction ffcn(v,f); // objective function
  SXFunction gfcn(v,g); // constraints
 
  
  // Allocate an NLP solver
  IpoptSolver solver(ffcn,gfcn);

  // Set options
  solver.setOption("tol",1e-6);
  solver.setOption("generate_hessian",true);

  // Initialize the solver
  solver.init();

  // State evolution constraints (difference equation)
  vector<double> xEq(nx,0.);
  	
  	
  solver.setInput(xEq,"lbg");
  solver.setInput(xEq,"ubg");
  
  // Bounds on u and u, and initial guess
  vector<double> vmin(nu+nx), vmax(nu+nx), vinit(nu+nx);

  for(int i=0; i<nu+nx; ++i){ // Control and state bounds
    vmin[i] = -1;
    vmax[i] =  1;
  }
  
  for(int i=0; i<nu; ++i) { // Initial guess on control
  	vinit[i]=-0.2;
  }
  
  vinit[0]=0.2; // Initial guess on state
  for(int i=1; i<nx; ++i){  // Forward simulation for initial guess on state
  	vinit[i+nu]=vinit[i+nu-1]+0.1*(vinit[i+nu-1]*(1+vinit[i+nu-1])+vinit[i]);
  }


  solver.setInput(vmin,"lbx");
  solver.setInput(vmax,"ubx");
  solver.setInput(vinit,"x0");

  // Solve the problem
  solver.solve();

  // Print the optimal cost
  double cost;
  solver.getOutput(cost,"f");
  cout << "optimal cost: " << cost << endl;

  // Print the optimal solution
  vector<double> varopt(nu+nx), uopt(nu), xopt(nx);
  solver.getOutput(varopt,"x");
  
  for(int i=0; i<nu; ++i)
  	uopt[i]=varopt[i];
  
  for(int i=0; i<nx; ++i)
  	xopt[i]=varopt[i+nu];
  	
  cout << "Optimal control: " << uopt << endl;
  cout << "Optimal state: " << xopt << endl;
  
  return 0;
}
