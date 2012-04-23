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
from numpy import *

x = ssym("x",5)
f = sin(x) + x
print "f = ", f

fcn = SXFunction([x],[f])
fcn.init()

[f_recon,J] = fcn.jac([(0,-1),(0,0)])
print "f_recon = ", f_recon
print "J = ", J

x = SX("x")
y = SX("y")
z = SX("z")
w = SX("w")
f = [x, (x+(2*(y*y))), ((x+(2*(y*(y*y))))+(3*((z*z)*(z*z)))), w]
print "f = ", f

fcn = SXFunction([[x,y,z,w]],[f])
fcn.init()

[f_recon,J] = fcn.jac([(0,-1),(0,0)])
print "f_recon = ", f_recon
print "J = ", J



##q = SXMatrix(1,1,SX(5)) # segfault
#q = SXVector(1,SX(5)) # segfault
##q = DVector(1,5) # works
#d = q[:][0]
#print d
