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

# Time
t = ssym("t")

# States
x = ssym("x",2)

# control
u = ssym("u")

# algebraic variable
z = ssym("z")

# DAE
f_x = vertcat(( x[1], z - x[0] + u ))
f_z = z + (x[0]**2 - 1)*x[1]
f_q = x[0]**2 + x[1]**2 + u**2

# DAE callback function
f = SXFunction([x,z,u,t],[f_x,f_z,f_q])

# Create an integrator
DT = 0.5   # Length of a control interval
I = IdasIntegrator(f)
I.setOption("tf",DT)
I.init()

# All controls
N =  20    # Number of control intervals
U = msym("U",N)

# Eliminate all intermediate X:es and sum up cost contributions
X  = msym([1,0])
J = 0
for k in range(N):
  # Call the integrator
  (X,Q,_,_) = I.call( (X,U[k]) )
  
  # Add quadrature to objective
  J += Q
  
# Objective function: x_2(T)
jfcn = MXFunction([U],[J])

# Terminal constraints: x_0(T)=x_1(T)=0
gfcn = MXFunction([U],[X])

# Allocate an NLP solver
solver = IpoptSolver(jfcn,gfcn)
solver.setOption("generate_hessian",True)
solver.init()

# Variable bounds
solver.setInput(-0.8, NLP_LBX)
solver.setInput( 1.0, NLP_UBX)

# Constraint bounds
solver.setInput( 0.0, NLP_LBG)
solver.setInput( 0.0, NLP_UBG)

# Initial guess
solver.setInput( 0.0, NLP_X_INIT)

# Solve the problem
solver.solve()



import numpy as NP
import matplotlib.pyplot as plt

# Retrieve the solution
u_opt = NP.array(solver.output(NLP_X_OPT))

# End time
T =  N*DT  

# Time grid
tgrid_x = NP.linspace(0,T,N+1)
tgrid_u = NP.linspace(0,T,N)

# Plot the results
plt.figure(1)
plt.clf()
plt.plot(tgrid_u,u_opt,'-.')
plt.title("Van der Pol optimization - single shooting")
plt.xlabel('time')
plt.legend(['u trajectory'])
plt.grid()
plt.show()
