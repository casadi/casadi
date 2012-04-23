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

#solv = IpoptSolver(f,g)

solv = SQPMethod(f,g)
solv.setOption("qp_solver",IpoptQPSolver)

solv.init()

solv.setInput([2.5,3.0,0.75],NLP_X_INIT)
solv.setInput(0,NLP_UBG)
solv.setInput(0,NLP_LBG)
solv.solve()

print solv.output(NLP_X_OPT)


