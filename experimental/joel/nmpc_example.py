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

# Formulate the NLP
u = msym("u",30)     # Control
p = msym("p")        # Parameter
J = inner_prod(u,u)  # NLP objective
G = MX.zeros(0,1)    # NLP constraints
x = p                # State
for k in range(30):
    x = x +  0.1*(x*(x+1) + u[k])
    x.lift(0.0)      # Treat as NLP variable
    J = J + x*x      # Add to objective
    G.append(x)      # Add to constraints

# Setup the NLP solver
f = MXFunction([u,p],[J]) # Objective function
g = MXFunction([u,p],[G]) # Constraint function
S = SCPgen(f,g)      # NLP solver instance
S.setOption("parametric",True)
S.setOption("qp_solver",QPOasesSolver)
S.setOption("qp_solver_options",\
                    {"printLevel":"none"})
S.init()

# Pass bounds and solve the NLP
S.setInput(0.30, NLP_P)    # p
S.setInput(-1.0, NLP_LBX)  # u_min
S.setInput( 1.0, NLP_UBX)  # u_max
S.setInput(-1.0, NLP_LBG)  # x_min
S.setInput( 1.0, NLP_UBG)  # x_max
S.solve()

# Visualize the trajectory
from matplotlib.pylab import *
plot(S.output(NLP_X_OPT))
show()

