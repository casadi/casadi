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
#include <fstream>

using namespace casadi;
using namespace std;
/**
 *  Solve an NLP using codegen  
 *  Part 1: generation
 *  Joel Andersson, K.U. Leuven 2013
 */

int main(){
    
  /** Test problem 
   * 
   *    min x0^2 + x1^2
   *    s.t.    x0 + x1 - 10 = 0
   */

  // Optimization variables
  MX x = MX::sym("x",2);

  // Objective
  MX f = x[0]*x[0] + x[1]*x[1];

  // Constraints
  MX g = x[0]+x[1]-10;
    
  // Infinity
  double inf = numeric_limits<double>::infinity();

  // Create IPOPT instance
  Function solver = nlpsol("solver", "ipopt", {{"x", x}, {"f", f}, {"g", g}});

  // Generate C code for the NLP functions
  solver.generate_dependencies("nlp.c");

  // Generate Makefile
  ofstream makefile;
  makefile.open("./CMakeLists.txt");
  makefile << "cmake_minimum_required(VERSION 2.8.6)" << endl;
  makefile << "project(nlp-codegen-autogen C)" << endl;
  
  // Generate compilation instructions
  makefile << "add_library(nlp SHARED nlp.c)" << endl;
  makefile << "set_target_properties(nlp PROPERTIES PREFIX \"\")" << endl;
  makefile << "set_target_properties(nlp PROPERTIES SUFFIX \".casadi\")" << endl;
  makefile << endl;

  // Finalize makefile
  makefile.close();

  return 0;
}

