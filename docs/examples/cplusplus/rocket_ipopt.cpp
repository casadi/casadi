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

int main(){
  std::cout << "program started" << std::endl;

  // Dimensions
  int nu = 20;  // Number of control segments
  int nj = 100; // Number of integration steps per control segment

  // optimization variable
  SX u = SX::sym("u", nu); // control

  SX s_0 = 0; // initial position
  SX v_0 = 0; // initial speed
  SX m_0 = 1; // initial mass

  SX dt = 10.0/(nj*nu); // time step
  SX alpha = 0.05; // friction
  SX beta = 0.1; // fuel consumption rate

  // Trajectory
  std::vector<SX> s_k, v_k, m_k;

  // Integrate over the interval with Euler forward
  SX s = s_0, v = v_0, m = m_0;
  for(int k=0; k<nu; ++k){
    for(int j=0; j<nj; ++j){
      s += dt*v;
      v += dt / m * (u(k)- alpha * v*v);
      m += -dt * beta*u(k)*u(k);
    }
    s_k.push_back(s);
    v_k.push_back(v);
    m_k.push_back(m);
  }
  SX s_all=vertcat(s_k), v_all=vertcat(v_k), m_all=vertcat(m_k);

  // Objective function
  SX f = dot(u, u);

  // Terminal constraints
  SX g = vertcat(s, v, v_all);

  // Create the NLP
  SXDict nlp = {{"x", u}, {"f", f}, {"g", g}};

  // Allocate an NLP solver and buffers
  Function solver = nlpsol("solver", "ipopt", nlp);

  // Bounds on g
  std::vector<double> gmin = {10, 0};
  std::vector<double> gmax = {10, 0};
  gmin.resize(2+nu, -std::numeric_limits<double>::infinity());
  gmax.resize(2+nu, 1.1);

  // Solve the problem
  DMDict arg = {{"lbx", -10},
                     {"ubx", 10},
                     {"x0", 0.4},
                     {"lbg", gmin},
                     {"ubg", gmax}};
  DMDict res = solver(arg);

  // Print the optimal cost
  double cost(res.at("f"));
  std::cout << "optimal cost: " << cost << std::endl;

  // Print the optimal solution
  std::vector<double> uopt(res.at("x"));
  std::cout << "optimal control: " << uopt << std::endl;

  // Get the state trajectory
  Function xfcn("xfcn", {u}, {s_all, v_all, m_all});
  std::vector<double> sopt, vopt, mopt;
  xfcn({uopt}, {&sopt, &vopt, &mopt});
  std::cout << "position: " << sopt << std::endl;
  std::cout << "velocity: " << vopt << std::endl;
  std::cout << "mass:     " << mopt << std::endl;

  // Create Matlab script to plot the solution
  std::ofstream file;
  std::string filename = "rocket_ipopt_results.m";
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
