from casadi import *

x=SX.sym("x")
rx = SX.sym("rx")
tstart = SX.sym("tstart")
tend = SX.sym("tend")
f = SXFunction(daeIn(x=x),daeOut(ode=x))
f.init()

g = SXFunction(rdaeIn(x=x,rx=rx),rdaeOut(ode=1))
g.init()

integrator = CollocationIntegrator(f,g)
integrator.setOption("implicit_solver",KinsolSolver)
integrator.init()

integrator.setInput(1.1,"x0")

integrator.evaluate()

print integrator.getOutput("rxf")
