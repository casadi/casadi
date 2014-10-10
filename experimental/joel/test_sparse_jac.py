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
from numpy import *
from casadi import *
 
# create a simple function with sparse jacobian
x = ssym("x",10)
f = x+5
f[2] = 44
f[3] = 55
f[4:7] += inner_prod(x[1:4],x[1:4])
fcn = SXFunction([x],[f])

# Use forward or adjoint mode ad to calculate directional derivatives (uncomment one to force a mode)
#fcn.setOption("ad_mode","forward")
#fcn.setOption("ad_mode","reverse")

fcn.init()

# some point to evaluate the jacobian
x0 = linspace(-3.,12.,10)

# input index (i.e. with respect to which (matrix-valued) variable are we differentiating
iind = 0

# output index (i.e. with (matrix-valued) function output are we differentiating
oind = 0

# create the jacobian in the convensional way (i.e. using source code transformation)
J1 = fcn.jacobian(iind,oind)
J1.init()
J1.setInput(x0)
J1.evaluate()
print "J1.output().size() = ", J1.output().size()
print "J1.output().numel() = ", J1.output().numel()
print "J1(x0)", array(J1.getOutput())

# create the jacobian using compression techniques
J2 = Jacobian(fcn,iind,oind) # create a "Jacobian" function instance explicitly
J2.setOption("verbose",True) # so that it prints the number of directions
J2.init()
J2.setInput(x0)
J2.evaluate()
print "J2.output().size() = ", J2.output().size()
print "J2.output().numel() = ", J2.output().numel()
print "J2(x0)", array(J2.getOutput())

# Print difference
print "Difference: ", J2.getOutput()-J1.getOutput()
assert isEqual(J1.getOutput(),J2.getOutput())

