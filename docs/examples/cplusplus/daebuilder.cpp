/*
 *    MIT No Attribution
 *
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
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


#include <casadi/casadi.hpp>

using namespace casadi;

int main(){

  // Example on how to use the DaeBuilder class
  // Joel Andersson, UW Madison 2017

  // Start with an empty DaeBuilder instance
  DaeBuilder dae("rocket");

  // Add input expressions
  auto a = dae.add_p("a");
  auto b = dae.add_p("b");
  auto u = dae.add_u("u");
  auto h = dae.add_x("h");
  auto v = dae.add_x("v");
  auto m = dae.add_x("m");

  // Constants
  double g = 9.81; // gravity

  // Set ODE right-hand-side
  dae.set_ode("h", v);
  dae.set_ode("v", (u-a*pow(v,2))/m-g);
  dae.set_ode("m", -b*pow(u,2));

  // Specify initial conditions
  dae.set_start("h", 0);
  dae.set_start("v", 0);
  dae.set_start("m", 1);

  // Add meta information
  dae.set_unit("h","m");
  dae.set_unit("v","m/s");
  dae.set_unit("m","kg");

  // Print DAE
  dae.disp(std::cout, true);

  // Generate FMU
  std::vector<std::string> files = dae.export_fmu();
  std::cout << "generated files: " << files << "\n";
  return 0;
}
