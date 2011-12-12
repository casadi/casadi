#! Simulator
#! =====================
from casadi import *
from numpy import *
from pylab import *

#! Simulation without controls
#! ===========================
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
sim.input(SIMULATOR_X0).set([1,0])
sim.input(SIMULATOR_P).set([0.1,0.1,0.1,0.3,0.1])
print sim.input(SIMULATOR_P)
print sim.input(SIMULATOR_V)
sim.evaluate()

#! Plot the solution
#plot(array(sim.output())[:,0],array(sim.output())[:,1])
xlabel("u")
ylabel("u_dot")
show()

#! Simulation with controls
#! ===========================
#! We consider a very system, an excitated linear spring
#! m x'' + c x' + k x = u(t)

#! states
x  = SX("x")
dx = SX("dx")
v  = SX("v")
dv = SX("dv")

#! parameters
k = SX("k")
c = SX("c")
m = SX("m")

#! controls
u = SX("u")

rhs = vertcat([v - dx, (u -  c*v - k*x)/m - dv ])
f=SXFunction({'NUM': DAE_NUM_IN, DAE_T: t, DAE_YDOT: [dx,dv], DAE_Y: [x,v], DAE_P: [k,c,m,u]},[rhs])
f.init()

ts = linspace(0,50,10)

integrator = CVodesIntegrator(f)
integrator.setOption("monitor",["res"])
sim=Simulator(integrator,ts)
sim.setOption("np",3)
sim.setOption("nf",10)
sim.init()
sim.input(SIMULATOR_X0).set([0,0])
sim.input(SIMULATOR_P).set([1,0.1,1])
sim.input(SIMULATOR_V).set([0,-0.2,0,0.5,0,0,0,0.2,-0.8])
tsf = sim.getGrid()
sim.evaluate()

plot(tsf,array(sim.output())[:,0])
xlabel("t")
ylabel("x")

#! This example is plaged by a bug: http://sundials.2283335.n4.nabble.com/ReInit-functions-in-sundialsTB-td3239946.html

plot(ts[:-1],array(sim.input(SIMULATOR_V))[:,0])
show()

