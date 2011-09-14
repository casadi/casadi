# -*- coding: utf-8 -*-
from numpy import *
import matplotlib.pyplot as plt

# CasADi
from casadi import *

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
integrator.init()

# Initial conditions
xz0 = array([1.0, 1.000000])
integrator.setInput(xz0, INTEGRATOR_X0)
pq0 = array([1.0, 1.0])
integrator.setInput(pq0, INTEGRATOR_P)

# Integrate
integrator.evaluate()



print "Script finished"

