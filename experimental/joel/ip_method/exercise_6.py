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
import numpy as N

# Formulate the LP
n = 10
m = 8

# Initial guess
x0 = N.random.rand(n,1)+0.1
A = N.random.randn(m,n)
b = mul(A,x0)
c = N.random.randn(n)

# Formulate the NLP
x = ssym("x",n)
ffcn = SXFunction([x],[inner_prod(c,x)])
gfcn = SXFunction([x],[mul(A,x)])

for solver in ['ipopt','ip_method']:
  if solver=='ipopt':
    # Solve with IPOPT
    nlp_solver = IpoptSolver(ffcn,gfcn)
  elif solver =='ip_method':
    # Solve with IP-method
    nlp_solver = IPMethod(ffcn,gfcn)
    nlp_solver.setOption("linear_solver",CSparse)
    #nlp_solver.setOption("verbose",True)
  
  # Initialize the NLP solver
  nlp_solver.init()
  
  # Solve the NLP
  nlp_solver.setInput(x0,"x0")
  nlp_solver.setInput(N.inf,"ubx")
  nlp_solver.input("lbx").setZero()
  nlp_solver.setInput(b,"lbg")
  nlp_solver.setInput(b,"ubg")  
  nlp_solver.solve()
  
  # Print the solution
  print "Solution for ", solver, ": ", nlp_solver.getOutput("x_opt")
  
  # Check if feasible
  print "Residual = ", mul(A,nlp_solver.getOutput("x_opt"))-b








