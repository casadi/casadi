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
# interval length 0.5s
I = integrator('I', "idas", dae, 0, 0.5)

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
