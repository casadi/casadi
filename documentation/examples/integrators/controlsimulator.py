#! ControlSimulator
#! =====================
from casadi import *
from numpy import *
from pylab import *

#! The ControlSimulator is typically used for dynamic systems with piecewise constant control.
#! We consider a very system, an excitated linear spring 
#! m x'' + c x' + k x = u(t) 
	 
t = SX("t")

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

#! Create the dae
rhs = vertcat([v - dx, (u -  c*v - k*x)/m - dv ])
f=SXFunction({'NUM': DAE_NUM_IN, DAE_T: t, DAE_YDOT: [dx,dv], DAE_Y: [x,v], DAE_P: [k,c,m,u]},[rhs])
f.init()

#! Choose a time grid, we will have 10-1 = 9 control intervals
ts = linspace(0,50,10)

sim=ControlSimulator(f,ts)
sim.setOption("integrator",CVodesIntegrator)

#! Of the list [k,c,m,u] the first 3 are static parameters
sim.setOption("np",3) 
#! Each control interval will be subdived in 8
sim.setOption("nf",8) 
sim.init()
sim.input(PW_SIMULATOR_X0).set([0,0])
sim.input(PW_SIMULATOR_P).set([1,0.1,1])
#! Our 9 control intervals have the following prescribed values for u:
sim.input(PW_SIMULATOR_V).set([0,-0.2,0,0.5,0,0,0,0.2,-0.8]) 
sim.evaluate()

#! Obtain the fine time grid
tsf = sim.getGrid()

#! Plot the default output, i.e. the states
plot(tsf,array(sim.output())[:,0])
xlabel("t")
ylabel("x")

#! Plot the controls
plot(ts[:-1],array(sim.input(PW_SIMULATOR_V))[:,0],'o') # Sampled on the coarse grid
plot(tsf[:-1],array(sim.getVFine())[:,0],'.')           # Sampled on the fine grid 
legend(('x','u (coarse)','u (fine)'))
show()
