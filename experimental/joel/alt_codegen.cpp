/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSefcn.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include <iostream>
using namespace std;

// Define the source code to be both compiled and code generated
#define MY_SRC(_n) _n \
void hello_world(){ _n	      \
  cout << "Hello, world!" << endl; _n \
  /*this is a comment which will disappear in the generated code */ \
  cout << "Hello second line, world!" << endl; _n \
  cout << "Hello again, world!" << endl; _n \
} _n \

// Declare the source for compilation
MY_SRC()

// String trick
#define STR1(x) #x
#define STR(x) STR1(x)

// Source code as a string
const char* my_src_str = STR(MY_SRC(\n));

int main(){
  // Print the code
  cout << "This is the source: " << endl;
  cout << my_src_str << endl;

  // Execute the code
  cout << "This is the execution: " << endl;
  hello_world();
    
  return 0;
}

