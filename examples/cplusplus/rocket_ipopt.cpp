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
#include <symbolic/casadi.hpp>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <symbolic/stl_vector_tools.hpp>

using namespace CasADi;
using namespace std;

int main(){
    
  cout << "program started" << endl;
      
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
  IpoptSolver solver(ffcn,gfcn);

  // initialize the solver
  solver.init();

  // Bounds on u and initial condition
  vector<double> umin(nu), umax(nu), uinit(nu);
  for(int i=0; i<nu; ++i){
    umin[i] = -10;
    umax[i] =  10;
    uinit[i] = 0.4;
  }
  solver.setInput(umin,"lbx");
  solver.setInput(umax,"ubx");
  solver.setInput(uinit,"x0");
  
  // Bounds on g
  vector<double> gmin(2), gmax(2);
  gmin[0] = gmax[0] = 10;
  gmin[1] = gmax[1] =  0;
  gmin.resize(2+nu, -numeric_limits<double>::infinity());
  gmax.resize(2+nu, 1.1);
  
  solver.setInput(gmin,"lbg");
  solver.setInput(gmax,"ubg");

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

