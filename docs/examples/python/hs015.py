#
#     MIT No Attribution
#
#     Copyright (C) 2010-2026 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
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
import casadi as cs

"""
Solve the hs015 problem, formulated as the NLP:

minimize     100*(x[1] - x[0]**2)**2 + (1 - x[0])**2
subject to   x[0]*x[1] >= 1
             x[0] + x[1]**2 >= 0
             x[0] <= 1/2

David Kiessling, 2026
"""

# Declare variable x
x = cs.MX.sym('x', 2)

# Formulate the NLP
f = 100*(x[1] - x[0]**2)**2 + (1 - x[0])**2
g = cs.MX.zeros(2, 1)
g[0] = x[0]*x[1]
g[1] = x[0] + x[1]**2

ubx = cs.vertcat(1/2, cs.inf)
lbx = cs.vertcat(-cs.inf, -cs.inf)
lbg = cs.vertcat(1, 0)
ubg = cs.vertcat(cs.inf, cs.inf)
nlp = {'x':x, 'f':f, 'g':g}

# Create an NLP solver
solver = cs.nlpsol("solver", "uno", nlp)

# Solve the problem
res = solver(x0  = [-2.0,1.0],
             ubg = ubg,
             lbg = lbg,
             ubx = ubx)

# Print solution
print()
print("%50s " % "Optimal cost:", res["f"])
print("%50s " % "Primal solution:", res["x"])
print("%50s " % "Dual solution (simple bounds):", res["lam_x"])
print("%50s " % "Dual solution (nonlinear bounds):", res["lam_g"])
