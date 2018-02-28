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


#include <casadi/casadi.hpp>
#include <string.h>
/** Solve a simple QP
    We want to model a chain attached to two supports and hanging in between. Let us discretise
    it with N mass points connected by N-1 springs. Each mass i has position (yi,zi), i=1,...,N.
    The equilibrium point of the system minimises the potential energy.

    The potential energy of each spring is
    Vi=D_i/2 * ((y_i-y_{i+1})^2 + (z_i-z_{i+1})^2)

    The gravitational potential energy of each mass is
    Vg_i = m_i*g0*z_i

    The total potential energy is thus given by:

    Vchain(y,z) = 1/2*sum{i=1,...,N-1} D_i ((y_i-y_{i+1})^2+(z_i-z_{i+1})^2) + g0 * sum{i=1,...,N} m_i * z_i

    where y=[y_1,...,y_N] and z=[z_1,...,z_N]

    We wish to solve
    minimize{y,z} Vchain(y, z)

    Subject to the piecewise linear ground constraints:
    z_i >= zin
    z_i - 0.1*y_i >= 0.5

    Joel Andersson, 2015
*/

using namespace casadi;
using namespace std;

int main(int argc, char* argv[]) {
  // Default options
  string plugin_name = "qpoases";
  bool schur = false;
  bool nlp = false;
  bool large = false;
  // Read options, first is plugin name
  if (argc>1) plugin_name = argv[1];
  for (int i=2; i<argc; ++i) {
    if (strcmp(argv[i], "schur")==0) schur = true;
    if (strcmp(argv[i], "nlp")==0) nlp = true;
    if (strcmp(argv[i], "large")==0) large = true;
  }

  // Constants
  int N = large ? 200 : 40;
  double m_i = 40.0/N;
  double D_i = 70.0*N;
  double g0 = 9.81;
  //double zmin = -inf; // unbounded
  double zmin = 0.5; // ground

  // Objective function
  SX Vchain = 0;

  // Variables
  vector<SX> x;

  // Variable bounds
  vector<double> lbx, ubx;

  // Constraints
  vector<SX> g;

  // Constraint bounds
  vector<double> lbg, ubg;

  // Loop over all chain elements
  SX y_prev, z_prev, y_i, z_i;
  for (int i=1; i<N+1; ++i) {
    // Previous point
    if (i>1) {
      y_prev = y_i;
      z_prev = z_i;
    }

    // Create variables for the (y_i, z_i) coordinates
    y_i = SX::sym("y_" + str(i));
    z_i = SX::sym("z_" + str(i));

    // Add to the list of variables
    x.push_back(y_i);
    x.push_back(z_i);
    if (i==1) {
      lbx.push_back(-2);
      ubx.push_back(-2);
      lbx.push_back( 1);
      ubx.push_back( 1);
    } else if (i==N) {
      lbx.push_back( 2);
      ubx.push_back( 2);
      lbx.push_back( 1);
      ubx.push_back( 1);
    } else {
      lbx.push_back(-inf);
      ubx.push_back( inf);
      lbx.push_back(zmin);
      ubx.push_back( inf);
    }

    // Spring potential
    if (i>1) {
      Vchain += D_i/2*(sq(y_prev-y_i) + sq(z_prev-z_i));
    }

    // Graviational potential
    Vchain += g0 * m_i * z_i;

    // Slanted ground constraints
    g.push_back(z_i - 0.1*y_i);
    lbg.push_back( 0.5);
    ubg.push_back( inf);
  }

  // Formulate QP
  SXDict qp = {{"x", vertcat(x)}, {"f", Vchain}, {"g", vertcat(g)}};

  // Solver specific options
  Dict solver_options;
  if (plugin_name == "qpoases") {
    solver_options["sparse"] = true;
    solver_options["schur"] = schur;
    solver_options["print_time"] = true;
  }

  // Create solver instance
  Function solver = qpsol("solver", plugin_name, qp, solver_options);

  // Get the optimal solution
  DMDict arg = {{"lbx", lbx},
                {"ubx", ubx},
                {"lbg", lbg},
                {"ubg", ubg}};
  DMDict res = solver(arg);
  DM x_opt = res["x"];
  double f_opt(res["f"]);
  cout << "f_opt = " << f_opt << endl;

  // Retrieve the result
  DM y_opt = x_opt(Slice(0, x_opt.size1(), 2));
  DM z_opt = x_opt(Slice(1, x_opt.size1(), 2));

  // Create Matlab script to plot the solution
  ofstream file;
  string filename = "chain_qp_results.m";
  file.open(filename.c_str());
  file << "% Results file from " __FILE__ << endl;
  file << "% Generated " __DATE__ " at " __TIME__ << endl;
  file << endl;
  file << "t = linspace(-2,2," << N << ");"<< endl;
  file << "f_opt = " << f_opt << ";" << endl;
  file << "x_opt = " << x_opt << ";" << endl;
  file << "y_opt = " << y_opt << ";" << endl;
  file << "z_opt = " << z_opt << ";" << endl;

  // Finalize the results file
  file << endl;
  file << "% Plot the results" << endl;
  file << "figure(1);" << endl;
  file << "clf;" << endl;
  file << "plot(y_opt, z_opt);" << endl;
  file << "grid on;" << endl;
  file << "xlabel('y [m]');" << endl;
  file << "ylabel('z [m]');" << endl;

  file.close();
  cout << "Results saved to \"" << filename << "\"" << endl;

  return 0;
}
