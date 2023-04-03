/*
 *    MIT No Attribution
 *
 *    Copyright 2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
 *
 *    Permission is hereby granted, free of charge, to any person obtaining a copy of this
 *    software and associated documentation files (the "Software"), to deal in the Software
 *    without restriction, including without limitation the rights to use, copy, modify,
 *    merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
 *    permit persons to whom the Software is furnished to do so.
 *
 *    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 *    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 *    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 *    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */


#include <iostream>
#include <fstream>
#include <ctime>
#include <casadi/casadi.hpp>

using namespace casadi;

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
  opti.subject_to(pos(0)==0);   // start at position 0 ...
  opti.subject_to(speed(0)==0); // ... from stand-still 
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
  std::ofstream file;
  std::string filename = "race_car_results.m";
  file.open(filename.c_str());
  file << "% Results file from " __FILE__ << std::endl;
  file << "% Generated " __DATE__ " at " __TIME__ << std::endl;
  file << std::endl;
  
  // Save results to file
  file << "t = linspace(0," << sol.value(T) << "," << N << "+1);"<< std::endl;
  file << "speed = " << std::vector<double>(sol.value(speed)) << ";" << std::endl;
  file << "pos = " << std::vector<double>(sol.value(pos)) << ";" << std::endl;
  file << "U = " << std::vector<double>(sol.value(U)) << ";" << std::endl;

  file << "figure;" << std::endl;
  file << "hold on;" << std::endl;
  file << "plot(t,speed);" << std::endl;
  file << "plot(t,pos);" << std::endl;
  file << "plot(t,1-sin(2*pi*pos)/2,'r--');" << std::endl;
  file << "stairs(t(1:end-1),U,'k');" << std::endl;
  file << "xlabel('Time [s]');" << std::endl;
  file << "legend('speed','pos','speed limit','throttle','Location','northwest');" << std::endl;

  // Have a look at the constraint Jacobian
  jacobian(opti.g(), opti.x()).sparsity().spy_matlab("race_car_jac_g.m");
  
  file << "figure" << std::endl;
  file << "race_car_jac_g;" << std::endl;
  file << "xlabel('decision variables');" << std::endl;
  file << "ylabel('constraints');" << std::endl;
  file << "print('jac_sp','-dpng');" << std::endl;

  return 0;
}
