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
This example demonstrates how NL-files, which can be generated
by AMPl or Pyomo, can be imported in CasADi and solved using
e.g. the interface to AMPL

Joel Andersson
2012

"""

# Create an NLP instance
nl = NlpBuilder()

# Parse an NL-file
nl.import_nl("../nl_files/hs107.nl",{"verbose":False})

# NLP solver options
opts = {}
opts["expand"] = True
# opts["verbose"] = True
# opts["ipopt"] = dict(max_iter=10, linear_solver="ma57", hessian_approximation="limited-memory")

# Create an NLP solver
nlpsol = nlpsol("nlpsol", "ipopt", nl, opts)

# Solve NLP
res = nlpsol(lbx=nl.x_lb,
             ubx=nl.x_ub,
             lbg=nl.g_lb,
             ubg=nl.g_ub,
             x0=nl.x_init)
