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
 * 		
 *		Attempt to solve the OCP: (Sequential approach; only control as variable)
 *
 *			min sum (x_i ^2 + u_i ^2)
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
  int nu = 30;  // Number of control segments

  // Optimization variable
  vector<SX> u = ssym("u",nu),data(); // control

  double x_0 = 0.2; // Initial state
  
  // Trajectory
  vector<SX> x_traj(nu);

  SX x = x_0;
  for(int k=0; k<nu; ++k){
	x+=0.1*(x*(1+x)+u[k]);
    x_traj[k] = x;
  }

  // Objective function
  SX f = 0; // Could take x_0 into account, but since it is a constant...
  for(int i=0; i<u.size(); ++i)
    f += u[i]*u[i] + x_traj[i]*x_traj[i]; 
    
  // Constraints (xmin <= g <= xmax, see later)
  vector<SX> g;
  g.insert(g.end(),x_traj.begin(),x_traj.end());
  
  
  
  // Create the NLP
  SXFunction ffcn(u,f); // objective function
  SXFunction gfcn(u,g); // constraint
  
  // Allocate an NLP solver
  IpoptSolver solver(ffcn,gfcn);

  // Set options
  solver.setOption("tol",1e-6);
  solver.setOption("generate_hessian",true);

  // Initialize the solver
  solver.init();

  // Bounds on u and initial guess
  vector<double> umin(nu), umax(nu), uinit(nu);
  for(int i=0; i<nu; ++i){
    umin[i] = -1;
    umax[i] =  1;
    uinit[i] = -0.2;
  }
  solver.setInput(umin,"lbx");
  solver.setInput(umax,"ubx");
  solver.setInput(uinit,"x0");
  
  // Bounds on the state (xmin < g < xmax)
  vector<double> xmin(nu,-1.), xmax(nu,1.);
  
  solver.setInput(xmin,"lbg");
  solver.setInput(xmax,"ubg");

  // Solve the problem
  solver.solve();

  // Print the optimal cost
  double cost;
  solver.getOutput(cost,"f");
  cout << "optimal cost: " << cost << endl;

  // Print the optimal solution
  vector<double> uopt(nu);
  solver.getOutput(uopt,"x");
  cout << "optimal control: " << uopt << endl;

  // Get the state trajectory
  vector<double> xopt_evol(nu), xopt(1);
  xopt[0]=x_0;
  vector<Matrix<SX> > xfcn_in(1,u);
  vector<Matrix<SX> > xfcn_out(1);
  xfcn_out[0] = x_traj;
  SXFunction xfcn(xfcn_in,xfcn_out);
  xfcn.init();
  xfcn.setInput(uopt);
  xfcn.evaluate();
  xfcn.getOutput(xopt_evol,0);
  xopt.insert(xopt.end(),xopt_evol.begin(),xopt_evol.end());  
  cout << "Optimal state: " << xopt << endl;
  

  return 0;
}

