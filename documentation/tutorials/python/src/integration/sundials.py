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
#! CasADi tutorial
#! ==================
from numpy import *
import numpy
from casadi import *
from pylab import *
#! ODE integration
#! -----------------
#! Let's construct a simple Van der Pol oscillator.
u = ssym("u")
x = ssym("x")
y = ssym("y")
f  = SXFunction([vertcat([x,y]),u], [vertcat([(1-y*y)*x-y+u,x])])
#! Manipulate the function to adhere to the integrator's
#! input/output signature
#! f(time;states;parameters)
fmod=SXFunction(daeIn(x=f.inputExpr(0),p=f.inputExpr(1)),daeOut(ode=f.outputExpr(0)))
fmod.setOption("name","ODE right hand side")
#! Create the CVodesIntegrator
integrator = CVodesIntegrator(fmod)
#! The whole series of sundials options are available for the user
integrator.setOption("fsens_err_con",True)
integrator.setOption("quad_err_con",True)
integrator.setOption("abstol",1e-6)
integrator.setOption("reltol",1e-6)
tend=10
integrator.setOption("t0",0)
integrator.setOption("tf",tend)
integrator.init()
#$ The integrator is really just a special kind of FX. Assume that we have an ODE/DAE in either explicit form:
#$ \begin{verbatim}
#$   der(x) = fx(x,z,p,t)
#$ \end{verbatim}
#$  or in implicit form
#$ \begin{verbatim}
#$   0 = fx(x,z,p,t,der(x))
#$   0 = fz(x,z,p,t)
#$ \end{verbatim}
#$ 
#$ Also assume a fixed time interval [t0,tf] and that the
#$ state x is given at the t0.
#$ 
#$ An integrator can then be thought of as a mapping from the state at
#$ the initial time and the parameters to the state at the final time:
#$ \begin{verbatim}
#$   x(tf) = F(x(t0),p)
#$ \end{verbatim}
#$ 
#$ Integrators in CasADi are a little more general than this, allowing
#$ the user to augment the ODE/DAE with a set of quadrature equations:
#$ \begin{verbatim}
#$     der(q) = fq(x,z,p,t), q(t0) = 0
#$ \end{verbatim}
#$ 
#$ As well as another DAE with terminal value constraints depending on the
#$ solution of the forward DAE:
#$ \begin{verbatim}
#$           0 = gx(rx,rz,rp,x,z,p,t,der(rx))
#$           0 = gz(rx,rz,rp,x,z,p,t)
#$     der(rq) = gq(rx,rz,rp,x,z,p,t)
#$ \end{verbatim}
#$ 
#$ This gives us in total four inputs and four outputs:
#$ \begin{verbatim}
#$   [x(tf),q(tf),rx(t0),rq(t0)]  = F(x(t0),p,rx(0),rp)
#$ \end{verbatim}
#$ 
print isinstance(integrator,FX)
print "%d -> %d" % (integrator.getNumInputs(),integrator.getNumOutputs())
#! Setup the Integrator to integrate from 0 to t=tend, starting at [x0,y0]
#! The output of Integrator is the state at the end of integration.
#! There are two basic mechanisms two retrieve the whole trajectory of states:
#!  - method A, using reset + integrate
#!  - method B, using Simulator
#!
#! We demonstrate first method A:
ts=numpy.linspace(0,tend,100)
x0 = 0;y0 = 1
integrator.input("x0").set([x0,y0])
integrator.input("p").set(0)
integrator.evaluate()
integrator.reset(0)
	
#! Define a convenience function to acces x(t)
def out(t):
	integrator.integrate(t)
	return integrator.output().toArray()

sol = array([out(t) for t in ts]).squeeze()
	
#! Plot the trajectory
figure()
plot(sol[:,0],sol[:,1])
title('Van der Pol phase space')
xlabel('x')
ylabel('y')
show()

#! We demonstrate method B:
sim=Simulator(integrator,ts)
sim.init()
sim.input("x0").set([x0,y0])
sim.input("p").set(0)
sim.evaluate()

sol2 = sim.output().toArray()
#! sol and sol2 are exactly the same
print linalg.norm(sol-sol2)

#! Sensitivity for initial conditions
#! ------------------------------------
#$ Let's see how a perturbation $\delta x(0)$ on the initial state $x(0)=x'(0) + \delta x(0)$
#$ affects the solution at tend $x(tend)=x'(tend)+\delta x(tend)$.
#$ We plot the map $\delta x_0 \mapsto \delta x(tend) $

def out(dx0):
	integrator.input("x0").set([x0+dx0,y0])
	integrator.evaluate()
	return integrator.output().toArray()
dx0=numpy.linspace(-2,2,100)

out = array([out(dx) for dx in dx0]).squeeze()
	
dxtend=out[:,0]-sol[-1,0]

figure()
plot(dx0,dxtend)
grid()
title('Initial perturbation map')
xlabel('dx(0)')
ylabel('dx(tend)')
show()
#$ By definition, this mapping goes through the origin. In the limit of $dx0 \to 0$, this map is purely linear. The slope at the origin is exactly what we call 'sensitivity'
#

integrator.input("x0").set([x0,y0])
integrator.fwdSeed("x0").set([1,0])
integrator.evaluate(1,0)
A = integrator.fwdSens()[0]
plot(dx0,A*dx0)
legend(('True sensitivity','Linearised sensitivity'))
plot(0,0,'o')
show()
#$ A common approach to visualise sensitivity for initial conditions is to overlay the phase space solution with ellipses defined by the local jacobian $\frac{\partial x(t)}{\partial x(0)} = [\frac{dx_1(t)}{dx_1(0)}\quad\frac{dx_1(t)}{dx_2(0)};\frac{dx_2(t)}{dx_1(0)}\quad\frac{dx_2(t)}{dx_2(0)}]$
#! The interpetation is that a small initial circular patch of phase space evolves into ellipsoid patches at later stages.

def out(t):
	integrator.fwdSeed("x0").set([1,0])
        integrator.evaluate(1,0)
	A=integrator.fwdSens().toArray()
	integrator.fwdSeed("x0").set([0,1])
	integrator.evaluate(1,0)
	B=integrator.fwdSens().toArray()
	return array([A,B]).squeeze().T

circle = array([[sin(x),cos(x)] for x in numpy.linspace(0,2*pi,100)]).T

figure()
plot(sol[:,0],sol[:,1])
grid()
for i in range(10):
	J=out(ts[10*i])
	e=0.1*numpy.dot(J,circle).T+sol[10*i,:]
	plot(e[:,0],e[:,1],color='red')
	
show()


#J=integrator.jacobian(INTEGRATOR_X0,0)
#print J
#print type(J)
#J.init()
#J.input("t0").set(0)
#J.input("tf").set(tend)
#J.input("x0").set([x0,y0])
#J.input("p").set(0)
#J.evaluate()
#print J.output().toArray()

#! The figure reveals that perturbations perpendicular to the phase space trajectory shrink.

#! Symbolic intergator results
#! ---------------------------
#! Since CVodesIntegrator is just another FX, 
#! the usual CasADi rules for symbolic evaluation are active.
#!
#! We create an MX 'w' that contains the result of a time integration with:
#! - a fixed integration start time , t=0s
#! - a fixed integration end time, t=10s
#! - a fixed initial condition (1,0)
#! - a free symbolic input, held constant during integration interval
u=MX("u")
w, = integratorOut(integrator.call(integratorIn(x0=MX([1,0]),p=u)),"xf")

#! We construct an MXfunction and a python help function 'out'
f=MXFunction([u],[w])
f.init()

def out(u):
	f.setInput(u)
	f.evaluate()
	return f.output().toArray()

print out(0)
print out(1)

#! Let's plot the results
uv=numpy.linspace(-1,1,100)

out = array([out(i) for i in uv]).squeeze()
figure()
plot(uv,out)
grid()
title('Dependence of final state on input u')
xlabel('u')
ylabel('state')
legend(('x','y'))
show()



