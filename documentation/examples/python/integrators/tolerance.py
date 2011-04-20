#! Integrator tolerances
#! =====================
from casadi import *
from numpy import *
from pylab import *

t=SX("t")

x=SX("x") 
dx=SX("dx")

f=SXFunction({'NUM': ODE_NUM_IN, ODE_T: t, ODE_Y: [x,dx]},[[dx,-x]])
f.init()


tend = 2*pi*3
ts = linspace(0,tend,1000)

tolerances = [-10,-5,-4,-3,-2,-1]

figure()

for tol in tolerances:
  integrator = CVodesIntegrator(f)
  integrator.setOption("reltol",10.0**tol)
  integrator.setOption("abstol",10.0**tol)
  integrator.init()

  sim=Simulator(integrator,ts)
  sim.init()
  sim.input(SIMULATOR_X0).set([1,0])
  sim.evaluate()

  plot(ts,sim.output()[:,0],label="tol = 1e%d" % tol)

legend( loc='upper left')
xlabel("Time [s]")
ylabel("State x [-]")
show()


tolerances = logspace(-15,1,500)
endresult=[]

for tol in tolerances:
  integrator = CVodesIntegrator(f)
  integrator.setOption("reltol",tol)
  integrator.setOption("abstol",tol)
  integrator.init()
  integrator.input(INTEGRATOR_TF).set(tend)
  integrator.input(INTEGRATOR_X0).set([1,0])
  integrator.evaluate()
  endresult.append(integrator.output()[0])
  
figure()
loglog(tolerances,(array(endresult)-1),'b',label="Positive error")
loglog(tolerances,-(array(endresult)-1),'r',label="Negative error")
xlabel("Integrator relative tolerance")
ylabel("Error at the end of integration time")
legend(loc='upper left')
show()
