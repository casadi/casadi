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
This example mainly intended for CasADi presentations. It contains a compact
implementation of a direct single shooting method for DAEs using a minimal
number of CasADi concepts.

It solves the following optimal control problem (OCP) in differential-algebraic
equations (DAE):

minimize     integral_{t=0}^{10} x0^2 + x1^2 + u^2  dt
x0,x1,z,u

subject to   dot(x0) == z*x0-x1+u     \
             dot(x1) == x0             }  for 0 <= t <= 10
                   0 == x1^2 + z - 1  /
             x0(t=0) == 0
             x1(t=0) == 1
             x0(t=10) == 0
             x1(t=10) == 0
             -0.75 <= u <= 1  for 0 <= t <= 10

Note that other methods such as direct collocation or direct multiple shooting
are usually preferably to the direct single shooting method in practise.

Joel Andersson, 2012-2015
"""

# Declare variables
x = SX.sym("x",2) # Differential states
z = SX.sym("z")   # Algebraic variable
u = SX.sym("u")   # Control

# Differential equation
f_x = vertcat(z*x[0]-x[1]+u, x[0])

# Algebraic equation
f_z = x[1]**2 + z - 1

# Lagrange cost term (quadrature)
f_q = x[0]**2 + x[1]**2 + u**2

# Create an integrator
dae = {'x':x, 'z':z, 'p':u, 'ode':f_x, 'alg':f_z, 'quad':f_q}
opts = {"tf":0.5} # interval length
I = integrator('I', "idas", dae, opts)

# All controls
U = MX.sym("U", 20)

# Construct graph of integrator calls
X  = [0,1]
J = 0
for k in range(20):
  Ik = I(x0=X, p=U[k])
  X = Ik['xf']
  J += Ik['qf']   # Sum up quadratures

# Allocate an NLP solver
nlp = {'x':U, 'f':J, 'g':X}
opts = {"ipopt.linear_solver":"ma27"}
solver = nlpsol("solver", "ipopt", nlp, opts)

# Pass bounds, initial guess and solve NLP
sol = solver(lbx = -0.75, # Lower variable bound
             ubx =  1.0,  # Upper variable bound
             lbg =  0.0,  # Lower constraint bound
             ubg =  0.0,  # Upper constraint bound
             x0  =  0.0) # Initial guess
print(sol)
