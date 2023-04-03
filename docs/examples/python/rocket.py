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
from pylab import *

# Control
u = MX.sym("u")

# State
x = MX.sym("x",3)
s = x[0] # position
v = x[1] # speed
m = x[2] # mass

# ODE right hand side
sdot = v
vdot = (u - 0.05 * v*v)/m
mdot = -0.1*u*u
xdot = vertcat(sdot,vdot,mdot)

# ODE right hand side function
f = Function('f', [x,u],[xdot])

# Integrate with Explicit Euler over 0.2 seconds
dt = 0.01  # Time step
xj = x
for j in range(20):
  fj = f(xj,u)
  xj += dt*fj

# Discrete time dynamics function
F = Function('F', [x,u],[xj])

# Number of control segments
nu = 50 

# Control for all segments
U = MX.sym("U",nu) 
 
# Initial conditions
X0 = MX([0,0,1])

# Integrate over all intervals
X=X0
for k in range(nu):
  X = F(X,U[k])

# Objective function and constraints
J = mtimes(U.T,U) # u'*u in Matlab
G = X[0:2]     # x(1:2) in Matlab

# NLP
nlp = {'x':U, 'f':J, 'g':G}
 
# Allocate an NLP solver
opts = {"ipopt.tol":1e-10, "expand":True}
solver = nlpsol("solver", "ipopt", nlp, opts)
arg = {}

# Bounds on u and initial condition
arg["lbx"] = -0.5
arg["ubx"] =  0.5
arg["x0"] =   0.4

# Bounds on g
arg["lbg"] = [10,0]
arg["ubg"] = [10,0]

# Solve the problem
res = solver(**arg)

# Get the solution
plot(res["x"])
plot(res["lam_x"])
grid()
show()
