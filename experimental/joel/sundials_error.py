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

# Declare variables
x = ssym("x",2) # Differential states
z = ssym("z")   # Algebraic variable
u = ssym("u")   # Control

# Differential equation
f_x = vertcat((x[1], z*x[1]-x[0]+u ))

# Algebraic equation
f_z = x[0]**2 + z - 1

# DAE callback function
f = SXFunction(daeIn(x=x,z=z,p=u),daeOut(ode=f_x,alg=f_z))

# Create an integrator
I = IdasIntegrator(f)
I.setOption("tf",2.0) # interval length
I.init()

# All controls
U = msym("U",5)

# Construct graph of integrator calls
X  = msym([1,0])
for k in range(5):
  (X,_,_,_) = I.call( (X,U[k]) )   # Call the integrator
  
# NLP callback functions
jfcn = MXFunction([U],[inner_prod(U,U)]) # Objective 
gfcn = MXFunction([U],[X]) # Constraint

# Allocate an NLP solver
solver = IpoptSolver(jfcn,gfcn)
solver.setOption("generate_hessian",True)
solver.init()

# Pass bounds, initial guess and solve NLP
solver.setInput(-0.8, NLP_LBX)    # Lower variable bound
solver.setInput( 1.0, NLP_UBX)    # Upper variable bound
solver.setInput( 0.0, NLP_LBG)    # Lower constraint bound
solver.setInput( 0.0, NLP_UBG)    # Upper constraint bound
solver.setInput( 0.0, NLP_X_INIT) # Initial guess
solver.solve()
