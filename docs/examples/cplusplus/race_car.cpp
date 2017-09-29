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
#include <fstream>
#include <ctime>
#include <casadi/casadi.hpp>

using namespace casadi;
using namespace std;

// dx/dt = f(x,u)
MX f(const MX& x, const MX& u) {
  return vertcat(x(1), u-x(1));
}

int main(){


  // Car race along a track
  // ----------------------
  // An optimal control problem (OCP),
  // solved with direct multiple-shooting.
  //
  // For more information see: http://labs.casadi.org/OCP

  int N = 100; // number of control intervals

  auto opti = casadi::Opti(); // Optimization problem

  Slice all;
  // ---- decision variables ---------
  auto X = opti.variable(2,N+1); // state trajectory
  auto pos   = X(0,all);
  auto speed = X(1,all);
  auto U = opti.variable(1,N);   // control trajectory (throttle)
  auto T = opti.variable();      // final time

  // ---- objective          ---------
  opti.minimize(T); // race in minimal time

  // ---- dynamic constraints --------
  auto dt = T/N;
  for (int k=0;k<N;++k) {
    auto k1 = f(X(all,k),         U(all,k));
    auto k2 = f(X(all,k)+dt/2*k1, U(all,k));
    auto k3 = f(X(all,k)+dt/2*k2, U(all,k));
    auto k4 = f(X(all,k)+dt*k3,   U(all,k));
    auto x_next = X(all,k) + dt/6*(k1+2*k2+2*k3+k4);
    opti.subject_to(X(all,k+1)==x_next); // close the gaps 
  }

  // ---- path constraints -----------
  opti.subject_to(speed<=1-sin(2*pi*pos)/2); // track speed limit
  opti.subject_to(0<=U<=1);           // control is limited

  // ---- boundary conditions --------
  opti.subject_to(pos(1)==0);   // start at position 0 ...
  opti.subject_to(speed(1)==0); // ... from stand-still 
  opti.subject_to(pos(N)==1); // finish line at position 1

  // ---- misc. constraints  ----------
  opti.subject_to(T>=0); // Time must be positive

  // ---- initial values for solver ---
  opti.set_initial(speed, 1);
  opti.set_initial(T, 1);

  // ---- solve NLP              ------
  opti.solver("ipopt"); // set numerical backend
  auto sol = opti.solve();   // actual solve

  // Create Matlab script to plot the solution
  ofstream file;
  string filename = "race_car_results.m";
  file.open(filename.c_str());
  file << "% Results file from " __FILE__ << endl;
  file << "% Generated " __DATE__ " at " __TIME__ << endl;
  file << endl;
  
  // Save results to file
  file << "t = linspace(0," << sol.value(T) << "," << N << "+1);"<< endl;
  file << "speed = " << std::vector<double>(sol.value(speed)) << ";" << endl;
  file << "pos = " << std::vector<double>(sol.value(pos)) << ";" << endl;
  file << "U = " << std::vector<double>(sol.value(U)) << ";" << endl;

  file << "figure;" << endl;
  file << "hold on;" << endl;
  file << "plot(t,speed);" << endl;
  file << "plot(t,pos);" << endl;
  file << "plot(t,1-sin(2*pi*pos)/2,'r--');" << endl;
  file << "stairs(t(1:end-1),U,'k');" << endl;
  file << "xlabel('Time [s]');" << endl;
  file << "legend('speed','pos','speed limit','throttle','Location','northwest');" << endl;

  // Have a look at the constraint Jacobian
  jacobian(opti.g(), opti.x()).sparsity().spy_matlab("race_car_jac_g.m");
  
  file << "figure" << std::endl;
  file << "race_car_jac_g;" << std::endl;
  file << "xlabel('decision variables');" << std::endl;
  file << "ylabel('constraints');" << std::endl;
  file << "print('jac_sp','-dpng');" << std::endl;

  return 0;
}
