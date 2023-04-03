/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
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
#include <fstream>
#include <ctime>
#include <casadi/casadi.hpp>

using namespace casadi;

int main(){

  std::cout << "program started" << std::endl;

  // Dimensions
  int nu = 20;  // Number of control segments
  int nj = 100; // Number of integration steps per control segment

  // optimization variable
  SX u = SX::sym("u",nu); // control

  SX s_0 = 0; // initial position
  SX v_0 = 0; // initial speed
  SX m_0 = 1; // initial mass

  SX dt = 10.0/(nj*nu); // time step
  SX alpha = 0.05; // friction
  SX beta = 0.1; // fuel consumption rate

  // Trajectory
  SX s_traj = SX::zeros(nu);
  SX v_traj = SX::zeros(nu);
  SX m_traj = SX::zeros(nu);

  // Integrate over the interval with Euler forward
  SX s = s_0, v = v_0, m = m_0;
  for(int k=0; k<nu; ++k){
    for(int j=0; j<nj; ++j){
      s += dt*v;
      v += dt / m * (u(k)- alpha * v*v);
      m += -dt * beta*u(k)*u(k);
    }
    s_traj(k) = s;
    v_traj(k) = v;
    m_traj(k) = m;
  }

  // Objective function
  SX f = dot(u,u);

  // Terminal constraints
  SX g = vertcat(s, v, v_traj);

  // Create an NLP
  SXDict nlp = {{"x", u}, {"f", f}, {"g", g}};

  // Create an NLP solver and buffers
  Function solver = nlpsol("solver", "snopt", nlp);
  std::map<std::string, DM> arg, res;

  // Bounds on u and initial condition
  arg["lbx"] = -10;
  arg["ubx"] = 10;
  arg["x0"] = 0.4;

  // Bounds on g
  std::vector<double> gmin(2), gmax(2);
  gmin[0] = gmax[0] = 10;
  gmin[1] = gmax[1] =  0;
  gmin.resize(2+nu, -std::numeric_limits<double>::infinity());
  gmax.resize(2+nu, 1.1);
  arg["lbg"] = gmin;
  arg["ubg"] = gmax;

  // Solve the problem
  res = solver(arg);

  // Print the optimal cost
  double cost(res.at("f"));
  std::cout << "optimal cost: " << cost << std::endl;

  // Print the optimal solution
  std::vector<double> uopt(res.at("x"));
  std::cout << "optimal control: " << uopt << std::endl;

  // Get the state trajectory
  Function xfcn("xfcn", {u}, {s_traj, v_traj, m_traj});
  std::vector<double> sopt, vopt, mopt;
  xfcn({uopt}, {&sopt, &vopt, &mopt});
  std::cout << "position: " << sopt << std::endl;
  std::cout << "velocity: " << vopt << std::endl;
  std::cout << "mass:     " << mopt << std::endl;

  // Create Matlab script to plot the solution
  std::ofstream file;
  std::string filename = "rocket_snopt_results.m";
  file.open(filename.c_str());
  file << "% Results file from " __FILE__ << std::endl;
  file << "% Generated " __DATE__ " at " __TIME__ << std::endl;
  file << std::endl;
  file << "cost = " << cost << ";" << std::endl;
  file << "u = " << uopt << ";" << std::endl;

  // Save results to file
  file << "t = linspace(0,10.0," << nu << ");"<< std::endl;
  file << "s = " << sopt << ";" << std::endl;
  file << "v = " << vopt << ";" << std::endl;
  file << "m = " << mopt << ";" << std::endl;

  // Finalize the results file
  file << std::endl;
  file << "% Plot the results" << std::endl;
  file << "figure(1);" << std::endl;
  file << "clf;" << std::endl << std::endl;

  file << "subplot(2,2,1);" << std::endl;
  file << "plot(t,s);" << std::endl;
  file << "grid on;" << std::endl;
  file << "xlabel('time [s]');" << std::endl;
  file << "ylabel('position [m]');" << std::endl << std::endl;

  file << "subplot(2,2,2);" << std::endl;
  file << "plot(t,v);" << std::endl;
  file << "grid on;" << std::endl;
  file << "xlabel('time [s]');" << std::endl;
  file << "ylabel('velocity [m/s]');" << std::endl << std::endl;

  file << "subplot(2,2,3);" << std::endl;
  file << "plot(t,m);" << std::endl;
  file << "grid on;" << std::endl;
  file << "xlabel('time [s]');" << std::endl;
  file << "ylabel('mass [kg]');" << std::endl << std::endl;

  file << "subplot(2,2,4);" << std::endl;
  file << "plot(t,u);" << std::endl;
  file << "grid on;" << std::endl;
  file << "xlabel('time [s]');" << std::endl;
  file << "ylabel('Thrust [kg m/s^2]');" << std::endl << std::endl;

  file.close();
  std::cout << "Results saved to \"" << filename << "\"" << std::endl;

  return 0;
}
