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
import casadi
ocp = casadi.SymbolicOCP()
#ocp.parseFMI('modelDescription.xml')
ocp.parseFMI('modelDescription.xml',{'sort_equations':False,'eliminate_dependent':False})
ocp.sortType(True) # temporary solution: enables the new sorting
print ocp

x = ocp.variable('x')
x_start = ocp.variable('x_start')
u = ocp.variable('u')
u_cost = ocp.variable('u_cost')

print(x_start.getStart()), " == 1.0"
casadi.updateDependent(ocp)
print(u_cost.getStart()), " == 2.0?"
print(u.getMax()), " == 0.2"
x_start.setStart(2)
print(x_start.getStart()), " == 2.0"
print(x.getStart()), " == 2.0"
print(u.getMax()), " == 0.4"
print(u_cost.getStart()), " == 4.0"
u_cost.setStart(3) # Error since u_cost is a dependent parameter?
x.setStart(4) # Not sure what we want here! What should happen if x_start changes value afterwards?

# JModelica
#from pyjmi import CasadiModel
#from pymodelica import compile_fmux
#jn = compile_fmux("Parameters", "parameters.mop")
#model = CasadiModel(jn)
#import matplotlib.pyplot as plt
#res = model.optimize()
#x = res['x']
#u = res['u']
#time = res['time']
#plt.figure(1)
#plt.clf()
#plt.subplot(2, 1, 1)
#plt.plot(time, x)
#plt.xlabel('$t$')
#plt.ylabel('$x$')
#plt.grid(True)

#plt.subplot(2, 1, 2)
#plt.plot(time, u)
#plt.grid(True)
#plt.xlabel('$t$')
#plt.ylabel('$u$')
#plt.show()
