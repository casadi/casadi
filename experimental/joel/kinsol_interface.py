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

# Implicitly defined variables
x = SX("a") 
y = SX("b")

# Parameters
p = SX("p")

# Residual function 
impfun = SXFunction([[x,y],[p]],[[x*x - y,x + 2*y + p]])

# Implicitly defined function
fun = KinsolSolver(impfun) 
fun.setOption("linear_solver",CSparse)
fun.setOption("abstol",1e-10)
fun.init() 

# Give the parameter a value
fun.setInput(0.1)

# Evalute
fun.setOutput([1E-7,1E-7]) # initial guess to the implicitly defined variable ([x,y])
fun.evaluate()
print fun.getOutput()

# Change output initial guess and evaluate again
fun.setOutput([1,1],0) # initial guess to the implicitly defined variable ([x,y])
fun.evaluate()
print fun.getOutput()

# Give forward seeds
fun.setFwdSeed(1.0)

# Give adjoint seeds
fun.setAdjSeed([0.0,1.0])

# Evaluate with sensitivities
fun.evaluate(1,1)

# Print sensitivities
print "forward sensitivities = ", fun.getFwdSens()
print "adjoint sensitivities = ", fun.getAdjSens()



