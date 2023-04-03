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
#include <iomanip>

/**
* This example demonstrates how NL-files, which can be generated
* by AMPl or Pyomo, can be imported in CasADi and solved using
* e.g. the interface to AMPL

 \author Joel Andersson, Vyacheslav Kungurtsev
 \date 2013
*/


using namespace casadi;

int main(int argc, char **argv){

  // Get the problem
  std::string problem = (argc==2) ? argv[1] : "../docs/examples/nl_files/hs107.nl";

  // Parse an NL-file
  NlpBuilder nl;
  nl.import_nl(problem);

  // Set options
  Dict opts;
  opts["expand"] = true;
  //  opts["verbose"] = true;

  // Allocate NLP solver and buffers
  Function solver = nlpsol("nlpsol", "worhp", nl, opts);
  std::map<std::string, DM> arg, res;

  // Solve NLP
  arg["lbx"] = nl.x_lb;
  arg["ubx"] = nl.x_ub;
  arg["lbg"] = nl.g_lb;
  arg["ubg"] = nl.g_ub;
  arg["x0"] = nl.x_init;
  res = solver(arg);
  for (auto&& s : res) {
    std::cout << std::setw(10) << s.first << ": " << std::vector<double>(s.second) << std::endl;
  }

  return 0;
}
