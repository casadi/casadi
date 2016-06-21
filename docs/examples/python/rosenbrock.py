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
