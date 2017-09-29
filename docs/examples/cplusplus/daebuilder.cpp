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

using namespace std;
using namespace casadi;

int main(){

  // Example on how to use the DaeBuilder class
  // Joel Andersson, UW Madison 2017

  // Start with an empty DaeBuilder instance
  DaeBuilder dae;

  // Add input expressions
  auto a = dae.add_p("a");
  auto b = dae.add_p("b");
  auto u = dae.add_u("u");
  auto h = dae.add_x("h");
  auto v = dae.add_x("v");
  auto m = dae.add_x("m");

  // Constants
  double g = 9.81; // gravity

  // Add output expressions
  auto hdot = v;
  auto vdot = (u-a*pow(v,2))/m-g;
  auto mdot = -b*pow(u,2);
  dae.add_ode("hdot", hdot);
  dae.add_ode("vdot", vdot);
  dae.add_ode("mdot", mdot);

  // Specify initial conditions
  dae.set_start("h", 0);
  dae.set_start("v", 0);
  dae.set_start("m", 1);

  // Add meta information
  dae.set_unit("h","m");
  dae.set_unit("v","m/s");
  dae.set_unit("m","kg");

  // Print DAE
  dae.disp(cout, true);

  return 0;
}
