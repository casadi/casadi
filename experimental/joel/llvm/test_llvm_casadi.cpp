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

#include <core/casadi.hpp>

using namespace casadi;
using namespace std;

int main(){
  
  // Construct a simple function
  SX x1 = ssym("x1");
  SX x2 = ssym("x2");
  SX r1 = sin(x2);
  SX r2 = x1+5;
  
  // Input arguments
  vector<SX> F_in(2);
  F_in[0] = x1;
  F_in[1] = x2;
  
  // Output arguments
  vector<SX> F_out(2);
  F_out[0] = r1;
  F_out[1] = r2;
  
  // Create function
  SXFunction F(F_in,F_out);
  F.setOption("just_in_time",true);
  F.init();

  // Generate C code
  F.generateCode("test.c");
  
  // Pass inputs
  F.setInput(10,0);
  F.setInput(20,1);
  
  // Evaluate
  F.evaluate();
  
  // Print results
  cout << F.output(0) << endl;
  cout << F.output(1) << endl;
  
  // Print the LLVM IR
  F.print();
  
  return 0;
} 
