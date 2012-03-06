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
 
 
x0 = SX("x0")
v0 = SX("v0")

#! parameters 
k = SX("k") 
c = SX("c") 
m = SX("m") 
 
#! controls 
u = SX("u") 

#! Create the dae
rhs = vertcat([v - dx, (u -  c*v - k*x)/m - dv ])
f=SXFunction({'NUM': CONTROL_DAE_NUM_IN, CONTROL_DAE_T: t, CONTROL_DAE_YDOT: [dx,dv], CONTROL_DAE_Y: [x,v], CONTROL_DAE_P: [k,c,m], CONTROL_DAE_U: [u] ,CONTROL_DAE_Y_MAJOR: [x0,v0]},[rhs])
f.init()

#! Choose a time grid, we will have 10-1 = 9 control intervals
ts = linspace(0,50,10)

sim=ControlSimulator(f,ts)
sim.setOption("integrator",CVodesIntegrator)

#! Each control interval will be subdived in 8
sim.setOption("nf",8) 
sim.init()
sim.input(CONTROLSIMULATOR_X0).set([0,0])
sim.input(CONTROLSIMULATOR_P).set([1,0.1,1])
#! Our 9 control intervals have the following prescribed values for u:
sim.input(CONTROLSIMULATOR_U).set([0,-0.2,0,0.5,0,0,0,0.2,-0.8]) 
sim.evaluate()

#! Obtain the fine time grid
tsf = sim.getMinorT()

figure(1)
#! Plot the default output, i.e. the states
plot(tsf,array(sim.output())[:,0])
xlabel("t")
ylabel("x")

#! Plot the controls
plot(ts[:-1],array(sim.input(CONTROLSIMULATOR_U))[:,0],'o') # Sampled on the coarse grid
plot(tsf[:-1],array(sim.getMinorU())[:,0],'.')           # Sampled on the fine grid 
legend(('x','u (coarse)','u (fine)'))

show()


#! Custom output function
#! =======================

h=SXFunction({'NUM': CONTROL_DAE_NUM_IN, CONTROL_DAE_T: t, CONTROL_DAE_YDOT: [dx,dv], CONTROL_DAE_Y: [x,v], CONTROL_DAE_P: [k,c,m], CONTROL_DAE_U: [u] ,CONTROL_DAE_Y_MAJOR: [x0,v0]},[x0,u])
h.init()

sim=ControlSimulator(f,h,ts)
sim.setOption("integrator",CVodesIntegrator)

#! Each control interval will be subdived in 8
sim.setOption("nf",8) 
sim.init()
sim.input(CONTROLSIMULATOR_X0).set([0,0])
sim.input(CONTROLSIMULATOR_P).set([1,0.1,1])
#! Our 9 control intervals have the following prescribed values for u:
sim.input(CONTROLSIMULATOR_U).set([0,-0.2,0,0.5,0,0,0,0.2,-0.8]) 
sim.evaluate()

figure(1)

plot(tsf,array(sim.output(1)),'x') 
plot(tsf,array(sim.output(0)),'*') 
legend(('x','u (coarse)','u (fine)','u (output)','x0 (output)'),loc='lower left')
show()

#! Working with interpolation
#! ===========================

f=SXFunction({'NUM': CONTROL_DAE_NUM_IN, CONTROL_DAE_T: t, CONTROL_DAE_YDOT: [dx,dv], CONTROL_DAE_Y: [x,v], CONTROL_DAE_P: [k,c,m], CONTROL_DAE_U_INTERP: [u] ,CONTROL_DAE_Y_MAJOR: [x0,v0]},[rhs])
f.init()

ui = ssym("ui")

h=SXFunction({'NUM': CONTROL_DAE_NUM_IN, CONTROL_DAE_T: t, CONTROL_DAE_YDOT: [dx,dv], CONTROL_DAE_Y: [x,v], CONTROL_DAE_P: [k,c,m], CONTROL_DAE_U_INTERP: [ui]  ,CONTROL_DAE_Y_MAJOR: [x0,v0]},[x,ui])
h.init()

sim=ControlSimulator(f,h,ts)
sim.setOption("integrator",CVodesIntegrator)

#! Each control interval will be subdived in 8
sim.setOption("control_interpolation","linear")
sim.setOption("control_endpoint",True)
sim.setOption("nf",8) 
sim.init()
sim.input(CONTROLSIMULATOR_X0).set([0,0])
sim.input(CONTROLSIMULATOR_P).set([1,0.1,1])
#! CONTROLSIMULATOR_U is larger, it has a value at the end of the last control interval, such that interpolation can happen
sim.input(CONTROLSIMULATOR_U).set([0,-0.2,0,0.5,0,0,0,0.2,-0.8,0]) 
sim.evaluate()

#! Obtain the fine time grid
tsf = sim.getMinorT()

figure(2)
#! Plot the default output, i.e. the states
plot(tsf,array(sim.output())[:,0])
xlabel("t")
ylabel("x")

#! Plot the controls
plot(ts,array(sim.input(CONTROLSIMULATOR_U))[:,0],'o') # Sampled on the coarse grid
plot(tsf,sim.output(1),'-')           # Sampled on the fine grid 
legend(('x','u (coarse)','u (fine)'))

show()

