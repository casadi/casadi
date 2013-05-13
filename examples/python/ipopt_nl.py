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
This example demonstrates how NL-files, which can be generated
by AMPl or Pyomo, can be imported in CasADi and solved using
e.g. the interface to AMPL

Joel Andersson
2012

"""

# Create an NLP instance
nl = SymbolicNLP()

# Parse an NL-file
nl.parseNL("../nl_files/hs107.nl",{"verbose":False})

# NLP function
nlp = SXFunction(nlpIn(x=nl.x),nlpOut(f=nl.f,g=nl.g))
  
# NLP solver
nlp_solver = IpoptSolver(nlp)
  
# Set options
# nlp_solver.setOption("max_iter",10)
#nlp_solver.setOption("verbose",True)
# nlp_solver.setOption("linear_solver","ma57")
# nlp_solver.setOption("hessian_approximation","limited-memory")
  
# Initialize NLP solver
nlp_solver.init()
  
# Pass the bounds and initial guess
nlp_solver.setInput(nl.x_lb,"lbx")
nlp_solver.setInput(nl.x_ub,"ubx")
nlp_solver.setInput(nl.g_lb,"lbg")
nlp_solver.setInput(nl.g_ub,"ubg")
nlp_solver.setInput(nl.x_init,"x0")
  
# Solve NLP
nlp_solver.solve()

