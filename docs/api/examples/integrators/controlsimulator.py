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
#! ControlSimulator
#! =====================
from casadi import *
from numpy import *
from pylab import *

#! The ControlSimulator is typically used for dynamic systems with piecewise constant control.
#! We consider a very system, an excitated linear spring 
#! m x'' + c x' + k x = u(t) 
	 
t = SX.sym("t")

#! states 
x  = SX.sym("x") 
v  = SX.sym("v") 
 
 
x0 = SX.sym("x0")
v0 = SX.sym("v0")

#! parameters 
k = SX.sym("k") 
c = SX.sym("c") 
m = SX.sym("m") 
 
#! controls 
u = SX.sym("u") 

#! Create the dae
rhs = vertcat((v, (u -  c*v - k*x)/m ))
fin = controldaeIn(
    t = t,
    x = vertcat((x,v)),
    p = vertcat((k,c,m)),
    u = u,
    x_major = vertcat((x0,v0))
  )
f=SXFunction(fin,daeOut(ode=rhs))
f.init()

#! Choose a time grid, we will have 10-1 = 9 control intervals
ts = linspace(0,50,10)

sim=ControlSimulator(f,ts)
sim.setOption("integrator", "cvodes")

#! Each control interval will be subdived in 8
sim.setOption("nf",8) 
sim.init()
sim.setInput([0,0],"x0")
sim.setInput([1,0.1,1],"p")
#! Our 9 control intervals have the following prescribed values for u:
sim.setInput([0,-0.2,0,0.5,0,0,0,0.2,-0.8],"u") 
sim.evaluate()

#! Obtain the fine time grid
tsf = sim.getMinorT()

figure(1)
#! Plot the default output, i.e. the states
plot(tsf,sim.getOutput()[0,:].T)
xlabel("t")
ylabel("x")

#! Plot the controls
plot(ts[:-1],sim.getInput("u")[0,:].T,'o') # Sampled on the coarse grid
plot(tsf[:-1],sim.getMinorU()[0,:].T,'.')           # Sampled on the fine grid 
legend(('x','u (coarse)','u (fine)'))

show()


#! Custom output function
#! =======================

fin = controldaeIn(
    t = t,
    x = vertcat([x,v]),
    p = vertcat([k,c,m]),
    u = u,
    x_major = vertcat([x0,v0])
  )
h=SXFunction(fin,[x0,u])
h.init()

sim=ControlSimulator(f,h,ts)
sim.setOption("integrator", "cvodes")

#! Each control interval will be subdived in 8
sim.setOption("nf",8) 
sim.init()
sim.setInput([0,0],"x0")
sim.setInput([1,0.1,1],"p")
#! Our 9 control intervals have the following prescribed values for u:
sim.setInput([0,-0.2,0,0.5,0,0,0,0.2,-0.8],"u") 
sim.evaluate()

figure(1)

plot(tsf,sim.getOutput(1).T,'x') 
plot(tsf,sim.getOutput(0).T,'*') 
legend(('x','u (coarse)','u (fine)','u (output)','x0 (output)'),loc='lower left')
show()

#! Working with interpolation
#! ===========================

fin = controldaeIn(
    t = t,
    x = vertcat([x,v]),
    p = vertcat([k,c,m]),
    u = u,
    x_major = vertcat([x0,v0])
  )
  
f=SXFunction(fin,[rhs])
f.init()

ui = SX.sym("ui")

fin = controldaeIn(
    t = t,
    x = vertcat([x,v]),
    p = vertcat([k,c,m]),
    u = ui,
    x_major = vertcat([x0,v0])
  )
  
h=SXFunction(fin,[x,ui])
h.init()

sim=ControlSimulator(f,h,ts)
sim.setOption("integrator", "cvodes")

#! Each control interval will be subdived in 8
sim.setOption("control_interpolation","linear")
sim.setOption("control_endpoint",True)
sim.setOption("nf",8) 
sim.init()
sim.setInput([0,0],"x0")
sim.setInput([1,0.1,1],"p")
#! CONTROLSIMULATOR_U is larger, it has a value at the end of the last control interval, such that interpolation can happen
sim.setInput([0,-0.2,0,0.5,0,0,0,0.2,-0.8,0],"u") 
sim.evaluate()

#! Obtain the fine time grid
tsf = sim.getMinorT()

figure(2)
#! Plot the default output, i.e. the states
plot(tsf,sim.getOutput()[0,:].T)
xlabel("t")
ylabel("x")

#! Plot the controls
plot(ts,sim.getInput("u")[0,:].T,'o') # Sampled on the coarse grid
plot(tsf,sim.getOutput(1).T,'-')           # Sampled on the fine grid 
legend(('x','u (coarse)','u (fine)'))

show()

