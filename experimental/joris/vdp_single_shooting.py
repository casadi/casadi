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
import numpy as NP
import matplotlib.pyplot as plt

from casadi.tools import *

nk = 4    # Control discretization
tf = 10.0  # End time

# Declare variables (use scalar graph)
t  = ssym("t")    # time
u  = ssym("u")    # control
x  = ssym("x",3)  # state
xd = ssym("xd",3) # state derivative

# ODE right hand side
rhs = vertcat( [(1 - x[1]*x[1])*x[0] - x[1] + u, \
                x[0], \
                x[0]*x[0] + x[1]*x[1] + u*u] )

# DAE residual function
f = SXFunction([t,x,u,xd],[rhs-xd])

# Create an integrator
f_d = CVodesIntegrator(f)
f_d.setOption("abstol",1e-8) # tolerance
f_d.setOption("reltol",1e-8) # tolerance
f_d.setOption("steps_per_checkpoint",1000)
f_d.setOption("tf",tf/nk) # final time
f_d.init()

# All controls (use matrix graph)
U = msym("U",nk) # nk-by-1 symbolic variable

# The initial state (x_0=0, x_1=1, x_2=0)
X  = msym([0,1,0])

# State derivative (only relevant for DAEs)
Xp = msym([0,0,0])

# Build a graph of integrator calls
for k in range(nk):
  [X,Xp] = f_d.call([X,U[k],Xp])
  
dotsave(X,filename='single.pdf')
  
# Objective function: x_2(T)
F = MXFunction([U],[X[2]])

# Terminal constraints: x_0(T)=x_1(T)=0
X_01 = X[0:2] # first two components of X
G = MXFunction([U],[X_01])

# Allocate an NLP solver
solver = IpoptSolver(F,G)
solver.setOption("monitor",["eval_f"])
solver.init()

# Set bounds and initial guess
solver.setInput(-0.75*NP.ones(nk), "lbx")
solver.setInput(1.0*NP.ones(nk), "ubx")
solver.setInput(NP.zeros(nk),"x0")
solver.setInput(NP.zeros(2),"lbg")
solver.setInput(NP.zeros(2),"ubg")

# Solve the problem
solver.solve()

# Retrieve the solution
u_opt = NP.array(solver.getOutput("x"))

# Time grid
tgrid_x = NP.linspace(0,10,nk+1)
tgrid_u = NP.linspace(0,10,nk)

# Plot the results
plt.figure(1)
plt.clf()
plt.plot(tgrid_u,u_opt,'-.')
plt.title("Van der Pol optimization - single shooting")
plt.xlabel('time')
plt.legend(['u trajectory'])
plt.grid()
plt.show()
