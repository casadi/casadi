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

x=SX.sym("x")
rx = SX.sym("rx")
tstart = SX.sym("tstart")
tend = SX.sym("tend")
f = SXFunction(daeIn(x=x),daeOut(ode=x))
f.init()

g = SXFunction(rdaeIn(x=x,rx=rx),rdaeOut(ode=1))
g.init()

integrator = CollocationIntegrator(f,g)
integrator.setOption("implicit_solver",KinsolSolver)
integrator.init()

integrator.setInput(1.1,"x0")

integrator.evaluate()

print integrator.getOutput("rxf")
