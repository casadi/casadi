#! Simulator
#! =====================
from casadi import *
from numpy import *
from pylab import *

#! We will investigate the working of Simulator with the help of the parametrically exited Duffing equation:
#!
#$ $\ddot{u}+\dot{u}-\epsilon (2 \mu \dot{u}+\alpha u^3+2 k u \cos(\Omega t))$ with $\Omega = 2 + \epsilon \sigma$.

t = SX("t")

u = SX("u") 
v = SX("v") 

eps   = SX("eps")
mu    = SX("mu")
alpha = SX("alpha")
k     = SX("k")
sigma = SX("sigma")
Omega = 2 + eps*sigma

params = [eps,mu,alpha,k,sigma]
rhs    = [v,-u-eps*(2*mu*v+alpha*u**3+2*k*u*cos(Omega*t))]

f=SXFunction({'NUM': DAE_NUM_IN, DAE_T: t, DAE_Y: [u,v], DAE_P: params},[rhs])
f.init()

integrator = CVodesIntegrator(f)

#! We will simulate over 50 seconds, 1000 timesteps.
ts = linspace(0,50,1000)

sim=Simulator(integrator,ts)
sim.init()
sim.input(INTEGRATOR_X0).set([1,0])
sim.input(INTEGRATOR_P).set([0.1,0.1,0.1,0.3,0.1])
sim.evaluate()

#! Plot the solution
plot(array(sim.output())[:,0],array(sim.output())[:,1])
xlabel("u")
ylabel("u_dot")
show()

