#
#     MIT No Attribution
#
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
#
#     Permission is hereby granted, free of charge, to any person obtaining a copy of this
#     software and associated documentation files (the "Software"), to deal in the Software
#     without restriction, including without limitation the rights to use, copy, modify,
#     merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
#     permit persons to whom the Software is furnished to do so.
#
#     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#     INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#     PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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
    cmd_args = [compiler,"-fPIC","-shared"]+flags+["nlp.c","-o","nlp.so"]
    subprocess.run(cmd_args)

    # Create a new NLP solver instance from the compiled code
    solver = nlpsol("solver", "ipopt", "./nlp.so")

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
