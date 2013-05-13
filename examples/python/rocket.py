#
#     This file is part of CasADi.
# 
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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
from pylab import *

# Control
u = msym("u")

# State
x = msym("x",3)
s = x[0] # position
v = x[1] # speed
m = x[2] # mass

# ODE right hand side
sdot = v
vdot = (u - 0.05 * v*v)/m
mdot = -0.1*u*u
xdot = vertcat([sdot,vdot,mdot])

# ODE right hand side function
f = MXFunction([x,u],[xdot])
f.setOption("name","f")
f.init()

# Integrate with Explicit Euler over 0.2 seconds
dt = 0.01  # Time step
xj = x
for j in range(20):
  [fj] = f.call([xj,u])
  xj += dt*fj

# Discrete time dynamics function
F = MXFunction([x,u],[xj])
F.init()

# Number of control segments
nu = 50 

# Control for all segments
U = msym("U",nu) 
 
# Initial conditions
X0 = MX([0,0,1])

# Integrate over all intervals
X=X0
for k in range(nu):
  [X] = F.call([X,U[k]])

# Objective function and constraints
J = mul(U.T,U) # u'*u in Matlab
G = X[0:2]     # x(1:2) in Matlab

# NLP
nlp = MXFunction(nlpIn(x=U),nlpOut(f=J,g=G))
  
# Allocate an NLP solver
solver = IpoptSolver(nlp)
solver.setOption("tol",1e-10)
solver.setOption("expand",True)
solver.init()

# Bounds on u and initial condition
solver.setInput(-0.5, "lbx")
solver.setInput( 0.5, "ubx")
solver.setInput( 0.4, "x0")

# Bounds on g
solver.setInput([10,0],"lbg")
solver.setInput([10,0],"ubg")

# Solve the problem
solver.solve()

# Get the solution
plot(solver.output("x"))
plot(solver.output("lam_x"))
grid()
show()
