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
#! Integrator tolerances
#! =====================
from casadi import *
from numpy import *
from pylab import *

x=SX.sym("x") 
dx=SX.sym("dx")
states = vertcat([x,dx])

f=SXFunction(daeIn(x=states),daeOut(ode=vertcat([dx,-x])))
f.init()

tend = 2*pi*3
ts = linspace(0,tend,1000)

tolerances = [-10,-5,-4,-3,-2,-1]

figure()

for tol in tolerances:
  integrator = Integrator("cvodes", f)
  integrator.setOption("reltol",10.0**tol)
  integrator.setOption("abstol",10.0**tol)
  integrator.init()

  sim=Simulator(integrator,ts)
  sim.init()
  sim.setInput([1,0],"x0")
  sim.evaluate()

  plot(ts,array(sim.getOutput())[0,:].T,label="tol = 1e%d" % tol)

legend( loc='upper left')
xlabel("Time [s]")
ylabel("State x [-]")
show()


tolerances = logspace(-15,1,500)
endresult=[]

for tol in tolerances:
  integrator = Integrator("cvodes", f)
  integrator.setOption("reltol",tol)
  integrator.setOption("abstol",tol)
  integrator.setOption("tf",tend)
  integrator.init()
  integrator.setInput([1,0],"x0")
  integrator.evaluate()
  endresult.append(integrator.getOutput()[0])
  
figure()
loglog(tolerances,(array(endresult)-1),'b',label="Positive error")
loglog(tolerances,-(array(endresult)-1),'r',label="Negative error")
xlabel("Integrator relative tolerance")
ylabel("Error at the end of integration time")
legend(loc='upper left')
show()
