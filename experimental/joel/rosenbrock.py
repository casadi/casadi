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
import numpy as NP
import matplotlib.pyplot as plt

x = ssym("x")
y = ssym("y")
z = ssym("z")
v = vertcat([x,y,z])

f = SXFunction([v],[x**2 + 100*z**2])
g = SXFunction([v],[z + (1-x)**2 - y])

# Choose NLP solver
#nlp_solver = IpoptSolver
#nlp_solver = SQPMethod
nlp_solver = LiftedSQP

# Choose a qp solver
#qp_solver = QPOasesSolver
qp_solver = IpoptQPSolver
#qp_solver = OOQPSolver

# Set solver specific options
if qp_solver == QPOasesSolver:
  qp_solver_options = {"printLevel" : "none"}
elif qp_solver == IpoptQPSolver:
  qp_solver_options = {"print_level" : 0}
else:
  qp_solver_options = {}

# Create solver
solv = nlp_solver(f,g)

# Pass options
if nlp_solver == IpoptSolver:
  solv.setOption("generate_hessian",True)
else:
  solv.setOption("qp_solver",qp_solver)
  solv.setOption("qp_solver_options",qp_solver_options)

# Init solver  
solv.init()

# Solve the rosenbrock problem
solv.setInput([2.5,3.0,0.75],NLP_X_INIT)
solv.setInput(0,NLP_UBG)
solv.setInput(0,NLP_LBG)
solv.evaluate()

# Print solution
print
print 
print "%50s " % "Optimal cost:", solv.output(NLP_COST)
print "%50s " % "Primal solution:", solv.output(NLP_X_OPT)
print "%50s " % "Dual solution (simple bounds):", solv.output(NLP_LAMBDA_X)
print "%50s " % "Dual solution (nonlinear bounds):", solv.output(NLP_LAMBDA_G)


