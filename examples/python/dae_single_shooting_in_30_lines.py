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

"""
This example mainly intended for CasADi presentations. 
It implements direct single shooting for DAEs in only about 30 lines of code, 
using a minimal number of concepts.

We do not recommend using this example as a template for solving optimal control problems in practise.

Joel Andersson, 2012
"""

# Declare variables
x = ssym("x",2) # Differential states
z = ssym("z")   # Algebraic variable
u = ssym("u")   # Control
t = ssym("t")   # Time

# Differential equation
f_x = vertcat((z*x[0]-x[1]+u, x[0]))

# Algebraic equation
f_z = x[1]**2 + z - 1

# Lagrange cost term (quadrature)
f_q = x[1]**2 + x[1]**2 + u**2

# DAE callback function
f = SXFunction([x,z,u,t],[f_x,f_z,f_q])

# Create an integrator
I = IdasIntegrator(f)
I.setOption("tf",0.5) # interval length
I.init()

# All controls
U = msym("U",20)

# Construct graph of integrator calls
X  = msym([0,1])
J = 0
for k in range(20):
  (X,Q,_,_) = I.call( (X,U[k]) )   # Call the integrator
  J += Q                           # Sum up quadratures
  
# NLP callback functions
nlp = MXFunction(nlIn(x=U),nlOut(f=J,g=X))

# Allocate an NLP solver
solver = IpoptSolver(nlp)
solver.init()

# Pass bounds, initial guess and solve NLP
solver.setInput(-0.75, "lbx")    # Lower variable bound
solver.setInput( 1.0,  "ubx")    # Upper variable bound
solver.setInput( 0.0,  "lbg")    # Lower constraint bound
solver.setInput( 0.0,  "ubg")    # Upper constraint bound
solver.setInput( 0.0,  "x0")     # Initial guess
solver.solve()
