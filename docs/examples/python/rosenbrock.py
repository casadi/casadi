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

# Declare variables
x = SX.sym("x")
y = SX.sym("y")
z = SX.sym("z")
v = vertcat([x,y,z])

# Form NLP functions
nlp = SXFunction(nlpIn(x=v),nlpOut(f=x**2 + 100*z**2, g=z + (1-x)**2 - y))

# Choose NLP solver
nlp_solver = "ipopt"
#nlp_solver = "worhp"
#nlp_solver = "sqpmethod"
#nlp_solver = "scpgen"

# Choose a qp solver (for CasADi NLP methods)
#qp_solver = "qpoases"
#qp_solver_options = {"printLevel" : "none"}

#qp_solver = "nlp"
#qp_solver_options = {"nlp_solver":"ipopt", "nlp_solver_options": {"print_level" : 0}}

#qp_solver = ooqp"
#qp_solver_options = {}

# Create solver
solv = NlpSolver(nlp_solver, nlp)

# NLP solver options
if nlp_solver in ("sqpmethod", "scpgen"):
  solv.setOption("qp_solver",qp_solver)
  solv.setOption("qp_solver_options",qp_solver_options)
  solv.setOption("max_iter",5)
  
# Init solver  
solv.init()

# Solve the rosenbrock problem
solv.setInput([2.5,3.0,0.75],"x0")
solv.setInput(0,"ubg")
solv.setInput(0,"lbg")
solv.evaluate()

# Print solution
print
print 
print "%50s " % "Optimal cost:", solv.getOutput("f")
print "%50s " % "Primal solution:", solv.getOutput("x")
print "%50s " % "Dual solution (simple bounds):", solv.getOutput("lam_x")
print "%50s " % "Dual solution (nonlinear bounds):", solv.getOutput("lam_g")

