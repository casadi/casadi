# -*- coding: utf-8 -*-
from numpy import *
import matplotlib.pyplot as plt

# CasADi
from casadi import *

# End time
tf = 1.0

# Variables
t = ssym("t")
x = ssym("x")
z = ssym("z")
p = ssym("p")
q = ssym("q")

# Differential equation input argument
ffcn_in = SXMatrixVector(DAE_NUM_IN)
ffcn_in[DAE_T] = t
ffcn_in[DAE_Y] = vertcat((x,z))
ffcn_in[DAE_P] = vertcat((p,q))

# Differential equation output argument
ffcn_out = [vertcat((-p*x*x*z, \
                      q*q - z*z + 0.1*x))]

# Differential equation
ffcn = SXFunction(ffcn_in,ffcn_out)

# Create integrator
integrator = AcadoIntegrator(ffcn)
integrator.setOption("time_dependence",False)
integrator.setOption("num_algebraic",1)
integrator.setOption("num_grid_points",100)
integrator.setOption("tf",tf)
integrator.init()

# Initial conditions
xz0 = array([1.0, 1.000000])
pq0 = array([1.0, 1.0])
integrator.setInput(xz0, INTEGRATOR_X0)
integrator.setInput(pq0, INTEGRATOR_P)

# Seeds
integrator.setFwdSeed([0.,0.], INTEGRATOR_X0)
integrator.setFwdSeed([1.,0.], INTEGRATOR_P)
integrator.setAdjSeed([1.,0.], INTEGRATOR_XF)

# Integrate with forward and adjoint sensitivities
integrator.evaluate(1,0)
#integrator.evaluate(1,1) # NOTE ACADO does not support adjoint mode AD using interfaced functions

# Result
print "final state =           ", integrator.output(INTEGRATOR_XF).data(), " (XF)"
print "forward sensitivities = ", integrator.fwdSens(INTEGRATOR_XF).data(), " (XF)"
#print "adjoint sensitivities = ", integrator.adjSens(INTEGRATOR_X0).data(), " (X0), ", integrator.adjSens(INTEGRATOR_P).data(), " (P)"

# Create a simulator
tgrid = numpy.linspace(0,tf,100)
simulator = Simulator(integrator,tgrid)
simulator.init()
simulator.setInput(xz0, INTEGRATOR_X0)
simulator.setInput(pq0, INTEGRATOR_P)
simulator.evaluate()

plt.clf()
plt.plot(tgrid,simulator.output())
plt.legend(('differential state', 'algebraic state'))
plt.grid(True)
plt.show()

print "Script finished"

