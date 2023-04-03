#
#     MIT No Attribution
#
#     Copyright 2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
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

"""
Solve the Rosenbrock problem, formulated as the NLP:

minimize     x^2 + 100*z^2
subject to   z+(1-x)^2-y == 0

Joel Andersson, 2015
"""

# Declare variables
x = SX.sym("x")
y = SX.sym("y")
z = SX.sym("z")

# Formulate the NLP
f = x**2 + 100*z**2
g = z + (1-x)**2 - y
nlp = {'x':vertcat(x,y,z), 'f':f, 'g':g}

# Create an NLP solver
solver = nlpsol("solver", "ipopt", nlp)

# Solve the Rosenbrock problem
res = solver(x0  = [2.5,3.0,0.75],
             ubg = 0,
             lbg = 0)

# Print solution
print()
print("%50s " % "Optimal cost:", res["f"])
print("%50s " % "Primal solution:", res["x"])
print("%50s " % "Dual solution (simple bounds):", res["lam_x"])
print("%50s " % "Dual solution (nonlinear bounds):", res["lam_g"])
