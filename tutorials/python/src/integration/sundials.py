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
fmod = ODE_NUM_IN * [[]]
fmod[ODE_T] = SXMatrix(t)
fmod[ODE_Y] = f.getArgumentIn(0)
fmod[ODE_P] = f.getArgumentIn(1)
print fmod
fmod=SXFunction(fmod,[f.getArgumentOut(0)])
fmod.setOption("name","ODE right hand side")
fmod.setOption("ad_order",1)
#! Create the CVodesIntegrator
integrator = CVodesIntegrator(fmod)
integrator.setOption("ad_order",1)
integrator.setOption("fsens_err_con",True)
integrator.setOption("quad_err_con",True)
integrator.setOption("abstol",1e-6)
integrator.setOption("reltol",1e-6)
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
tend=10
ts=linspace(0,tend,100)
x0 = 0;y0 = 1
integrator.setInput(0,INTEGRATOR_T0);   
integrator.setInput([x0,y0],INTEGRATOR_X0);
integrator.setInput(0,INTEGRATOR_P); 

#! Define a convenience function to acces x(t)
def out(t):
	integrator.setInput(t,INTEGRATOR_TF);
	integrator.evaluate()
	return integrator.getOutputData()

sol = array(map(out,ts))
figure()
plot(sol[:,0],sol[:,1])
title('Van der Pol phase space')
xlabel('x')
ylabel('y')
show()
#! Sensitivity for initial conditions
#! ------------------------------------
#$ Let's see how a perturbation $\delta x(0)$ on the initial state $x(0)=x'(0) + \delta x(0)$
#$ affects the solution at tend $x(tend)=x'(tend)+\delta x(tend)$.
#$ We plot the map $\delta x_0 \mapsto \delta x(tend) $
def out(dx0):
	integrator.setInput([x0+dx0,y0],INTEGRATOR_X0)
	integrator.evaluate()
	return integrator.getOutputData()
dx0=linspace(-2,2,100)

out = array(map(out,dx0))
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

integrator.setInput([x0,y0],INTEGRATOR_X0)
integrator.setFwdSeed([1,0],INTEGRATOR_X0)
integrator.evaluate(1,0)
A = integrator.getFwdSens()[0]
print A
plot(dx0,A*dx0)
legend(('True sensitivity','Linearised sensitivity'))
plot(0,0,'o')
show()
#$ A common approach to visualise sensitivity for initial conditions is to overlay the phase space solution with ellipses defined by the local jacobian $\frac{\partial x(t)}{\partial x(0)} = [\frac{dx_1(t)}{dx_1(0)}\quad\frac{dx_1(t)}{dx_2(0)};\frac{dx_2(t)}{dx_1(0)}\quad\frac{dx_2(t)}{dx_2(0)}]$
#! The interpetation is that a small initial circular patch of phase space evolves into ellipsoid patches at later stages.

def out(t):
	integrator.setInput(t,INTEGRATOR_TF);
	integrator.setFwdSeed([1,0],INTEGRATOR_X0)
	integrator.evaluate(1,0)
	A=integrator.getFwdSens()
	integrator.setFwdSeed([0,1],INTEGRATOR_X0)
	integrator.evaluate(1,0)
	B=integrator.getFwdSens()
	return array([A,B]).T

circle = array([[sin(x),cos(x)] for x in linspace(0,2*pi,100)]).T

figure()
plot(sol[:,0],sol[:,1])
grid()
for i in range(10):
	J=out(ts[10*i])
	e=0.1*dot(J,circle).T+sol[10*i,:]
	plot(e[:,0],e[:,1],color='red')
	
show()


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
w=integrator([MX(0),MX(10),MX([1,0]),u])

#! We construct an MXfunction and a python help function 'out'
f=MXFunction([u],[w])
f.init()

def out(u):
	f.setInput(u)
	f.evaluate()
	return f.getOutputData()

print out(0)
print out(1)

#! Let's plot the results
uv=linspace(-1,1,100)

out = array(map(out,uv))
figure()
plot(uv,out)
grid()
title('Dependence of final state on input u')
xlabel('u')
ylabel('state')
legend(('x','y'))
show()



