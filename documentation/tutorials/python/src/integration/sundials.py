#! CasADi tutorial
#! ==================
from numpy import *
from casadi import *
#from matplotlib.pylab import *
from pylab import *
#! ODE integration
#! -----------------
#! Let's construct a simple Van der Pol oscillator.
u = SX("u")
x = SX("x")
y = SX("y")
f  = SXFunction([[x,y],[u]], [[(1-y*y)*x-y+u,x]])
#! Manipulate the function to adhere to the integrator's
#! input/output signature
#! f(time;states;parameters)
t = SX("t")
fmod=SXFunction({'NUM': DAE_NUM_IN, DAE_T: t, DAE_Y: f.inputSX(0), DAE_P: f.inputSX(1)},[f.outputSX(0)])
fmod.setOption("name","ODE right hand side")
#! Create the CVodesIntegrator
integrator = CVodesIntegrator(fmod)
#! The whole series of sundials options are available for the user
integrator.setOption("fsens_err_con",True)
integrator.setOption("quad_err_con",True)
integrator.setOption("abstol",1e-6)
integrator.setOption("reltol",1e-6)
tend=10
integrator.setOption("tf",tend)
integrator.init()
#! The integrator is really just a special kind of FX.
#$ Quoting the integrator.hpp header documentation:
#$ \begin{verbatim}
#$ An "integrator" is a function that solves
#$   an initial value problem (IVP) of the generic form:
#$   
#$   F(x,der(x),p,t) == 0
#$   x(t0) = x0
#$ 
#$   It has 4 inputs, initial time, final time,
#$   initial state (vector-valued) and parameter (vector-valued)
#$   and one output, the state at the final time. 
#$   In addition to this, the integrator provides some additional functionality,
#$   such as getting the value of the state and/or sensitivities
#$   at certain time points.
#$   Controls are assumed to be parametrized at this point. 
#$     
#$   \end{verbatim}
print isinstance(integrator,FX)
print "%d -> %d" % (integrator.getNumInputs(),integrator.getNumOutputs())
#! Setup the Integrator to integrate from 0 to t=tend, starting at [x0,y0]
#! The output of Integrator is the state at the end of integration.
#! There are two basic mechanisms two retrieve the whole trajectory of states:
#!  - method A, using reset + integrate
#!  - method B, using Simulator
#!
#! We demonstrate first method A:
ts=linspace(0,tend,100)
x0 = 0;y0 = 1
integrator.input(INTEGRATOR_X0).set([x0,y0])
integrator.input(INTEGRATOR_P).set(0)
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
sim.input(SIMULATOR_X0).set([x0,y0])
sim.input(SIMULATOR_P).set(0)
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
	integrator.input(INTEGRATOR_X0).set([x0+dx0,y0])
	integrator.evaluate()
	return integrator.output().toArray()
dx0=linspace(-2,2,100)

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

integrator.input(INTEGRATOR_X0).set([x0,y0])
integrator.fwdSeed(INTEGRATOR_X0).set([1,0])
integrator.evaluate(1,0)
A = integrator.fwdSens()[0]
plot(dx0,A*dx0)
legend(('True sensitivity','Linearised sensitivity'))
plot(0,0,'o')
show()
#$ A common approach to visualise sensitivity for initial conditions is to overlay the phase space solution with ellipses defined by the local jacobian $\frac{\partial x(t)}{\partial x(0)} = [\frac{dx_1(t)}{dx_1(0)}\quad\frac{dx_1(t)}{dx_2(0)};\frac{dx_2(t)}{dx_1(0)}\quad\frac{dx_2(t)}{dx_2(0)}]$
#! The interpetation is that a small initial circular patch of phase space evolves into ellipsoid patches at later stages.

integrator.reset() # start integration from time zero again

def out(t):
	integrator.setFinalTime(t)
	integrator.fwdSeed(INTEGRATOR_X0).set([1,0])
	integrator.integrate(t)
	A=integrator.fwdSens().toArray()
	integrator.fwdSeed(INTEGRATOR_X0).set([0,1])
	integrator.evaluate(1,0)
	B=integrator.fwdSens().toArray()
	return array([A,B]).squeeze().T

circle = array([[sin(x),cos(x)] for x in linspace(0,2*pi,100)]).T

figure()
plot(sol[:,0],sol[:,1])
grid()
for i in range(10):
	J=out(ts[10*i])
	e=0.1*dot(J,circle).T+sol[10*i,:]
	plot(e[:,0],e[:,1],color='red')
	
show()


#J=integrator.jacobian(INTEGRATOR_X0,0)
#print J
#print type(J)
#J.init()
#J.input(INTEGRATOR_T0).set(0)
#J.input(INTEGRATOR_TF).set(tend)
#J.input(INTEGRATOR_X0).set([x0,y0])
#J.input(INTEGRATOR_P).set(0)
#J.evaluate()
#print J.output().toArray()

#! The figure reveals that perturbations perpendicular to the phase space trajectory shrink.

#! Symbolic intergator results
#! ---------------------------
#! Since CVodesIntegrator is just another FX, 
#! the usual CasAdi rules for symbolic evaluation are active.
#!
#! We create an MX 'w' that contains the result of a time integration with:
#! - a fixed integration start time , t=0s
#! - a fixed integration end time, t=10s
#! - a fixed initial condition (1,0)
#! - a free symbolic input, held constant during integration interval
u=MX("u")
integrator.setFinalTime(tend)
w=integrator({'NUM': INTEGRATOR_NUM_IN, INTEGRATOR_X0: MX([1,0]), INTEGRATOR_P: u})

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
uv=linspace(-1,1,100)

out = array([out(i) for i in uv]).squeeze()
figure()
plot(uv,out)
grid()
title('Dependence of final state on input u')
xlabel('u')
ylabel('state')
legend(('x','y'))
show()



