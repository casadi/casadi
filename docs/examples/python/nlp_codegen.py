#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             K.U. Leuven. All rights reserved.
#     Copyright (C) 2011-2014 Greg Horn
#
#     CasADi is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 3 of the License, or (at your option) any later version.
#
#     CasADi is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public
#     License along with CasADi; if not, write to the Free Software
#     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
#
from casadi import *

# Test problem
#
#    min x^2 + y^2
#    s.t.    x + y - 10 = 0
#

# Optimization variables
x = MX.sym("x")
y = MX.sym("y")

# Objective
f = x*x + y*y

# Constraints
g = x+y-10

# Create an NLP problem structure
nlp = {"x": vertcat(x,y), "f": f, "g": g}

mode = "jit"

# Pick a compiler
compiler = "gcc"    # Linux
# compiler = "clang"  # OSX
# compiler = "cl.exe" # Windows

# Run this script in an environment that recognised the compiler as command.
# On Windows, the suggested way is to run this script from a "x64 Native Tools Command Promt for VS" (Requires Visual C++ components or Build Tools for Visual Studio, available from Visual Studio installer. You also need SDK libraries in order to access stdio and math.)

flags = ["-O3"] # Linux/OSX

for mode in ["jit","external"]:

  if mode=="jit":
    # By default, the compiler will be gcc or cl.exe
    jit_options = {"flags": flags, "verbose": True, "compiler": compiler}
    options = {"jit": True, "compiler": "shell", "jit_options": jit_options}
    
    # Create an NLP solver instance
    solver = nlpsol("solver", "ipopt", nlp, options)

  elif mode=="external":
    # Create an NLP solver instance
    solver = nlpsol("solver", "ipopt", nlp)

    # Generate C code for the NLP functions
    solver.generate_dependencies("nlp.c")

    import subprocess
    # On Windows, use other flags
    subprocess.Popen([compiler,"-fPIC","-shared"]+flags+["nlp.c","-o","nlp.so"]).wait()

    # Create a new NLP solver instance from the compiled code
    solver = nlpsol("solver", "ipopt", "nlp.so")

  arg = {}

  arg["lbx"] = -DM.inf()
  arg["ubx"] =  DM.inf()
  arg["lbg"] =  0
  arg["ubg"] =  0
  arg["x0"] = 0

  # Solve the NLP
  res = solver(**arg)

  # Print solution
  print("-----")
  print("objective at solution =", res["f"])
  print("primal solution =", res["x"])
  print("dual solution (x) =", res["lam_x"])
  print("dual solution (g) =", res["lam_g"])
