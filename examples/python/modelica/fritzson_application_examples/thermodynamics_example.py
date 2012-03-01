# -*- coding: utf-8 -*-
import os
import sys
import numpy as NP
from numpy import *
import matplotlib.pyplot as plt
import zipfile
import time

try:
  # JModelica
  from pymodelica import compile_jmu
  from pyjmi import JMUModel
  import pymodelica
except:
  print "This example requires a jmodelica installation"
  sys.exit(0)

# CasADi
from casadi import *

# Matplotlib interactive mode
#plt.ion()

# Compile Modelica code to XML
def comp(name):
  curr_dir = os.path.dirname(os.path.abspath(__file__))
  jmu_name = compile_jmu(name, curr_dir+"/thermodynamics_example.mo",'modelica','ipopt',{'generate_xml_equations':True, 'generate_fmi_me_xml':False})
  modname = name.replace('.','_')
  sfile = zipfile.ZipFile(curr_dir+'/'+modname+'.jmu','r')
  mfile = sfile.extract('modelDescription.xml','.')
  os.remove(modname+'.jmu')
  os.rename('modelDescription.xml',modname+'.xml')

# Compile the simplemost example (conservation of mass in control volume)
comp("BasicVolumeMassConservation")

# Read a model from XML
ocp = FlatOCP('BasicVolumeMassConservation.xml')

# Set options
ocp.setOption("make_explicit",True)

# Parse the model
ocp.init()

# Create an integrator
integrator = CVodesIntegrator(ocp.daeFcn())

# Output function
m = ocp.variable("m").var()
P = ocp.variable("P").var()
output_fcn_out = ocp.substituteDependents([m,P])
output_fcn_in = ocp.daeArg()
output_fcn = SXFunction(output_fcn_in,output_fcn_out)

# Create a simulator
grid = NP.linspace(0,1,100)
simulator = Simulator(integrator,output_fcn,grid)
simulator.init()

# Pass initial conditions
x0 = getStart(ocp.xd())
simulator.setInput(x0,INTEGRATOR_X0)

# Simulate
simulator.evaluate()
integrator.printStats()

# Plot
plt.figure(1)
plt.subplot(1,2,1)
plt.plot(grid,simulator.output())
plt.xlabel("t")
plt.ylabel("m(t)")
plt.title("c.f. Fritzson figure 15-6 (left)")

plt.subplot(1,2,2)
plt.plot(grid,simulator.output(1))
plt.xlabel("t")
plt.ylabel("P(t)")
plt.title("c.f. Fritzson figure 15-6 (right)")
plt.draw()

# Compile the next example (conservation of energy in control volume)
comp("BasicVolumeEnergyConservation")

# Allocate a parser and load the xml
ocp = FlatOCP('BasicVolumeEnergyConservation.xml')

# Set options
ocp.setOption("make_explicit",True)

# Parse the model
ocp.init()

# Create an integrator
integrator = CVodesIntegrator(ocp.daeFcn())

# Output function
T = ocp.variable("T").var()
output_fcn_out = ocp.substituteDependents([T])
output_fcn_in = ocp.daeArg()
output_fcn = SXFunction(output_fcn_in,output_fcn_out)

# Create a simulator
grid = NP.linspace(0,10,100)
simulator = Simulator(integrator,output_fcn,grid)
simulator.init()

# Pass initial conditions
x0 = getStart(ocp.xd())
simulator.setInput(x0,INTEGRATOR_X0)

# Simulate
simulator.evaluate()
integrator.printStats()

# Plot
plt.figure(2)
plt.plot(grid,simulator.output())
plt.xlabel("t")
plt.ylabel("T(t)")
plt.title("c.f. Fritzson figure 15-9")
plt.draw()

# Compile the next example (Heat transfer and work)
comp("BasicVolumeTest")

# Allocate a parser and load the xml
ocp = FlatOCP('BasicVolumeTest.xml')

# Set options
ocp.setOption("make_explicit",True)

# Parse the model
ocp.init()

# Create an integrator
integrator = CVodesIntegrator(ocp.daeFcn())

# Output function
T = ocp.variable("T").var()
U = ocp.variable("U").var()
V = ocp.variable("V").var()
output_fcn_out = ocp.substituteDependents([T,U,V])
output_fcn_in = ocp.daeArg()
output_fcn = SXFunction(output_fcn_in,output_fcn_out)

# Create a simulator
grid = NP.linspace(0,2,100)
simulator = Simulator(integrator,output_fcn,grid)
simulator.init()

# Pass initial conditions
x0 = getStart(ocp.xd())
simulator.setInput(x0,INTEGRATOR_X0)

# Simulate
simulator.evaluate()
integrator.printStats()

# Plot
plt.figure(3)
p1, = plt.plot(grid,simulator.output(0))
p2, = plt.plot(grid,simulator.output(1))
plt.xlabel("t")
plt.ylabel("T(t)")
plt.legend([p2, p1], ["T", "U"])
plt.title("c.f. Fritzson figure 15-14")

plt.figure(4)
plt.plot(grid,simulator.output(2))
plt.xlabel("t")
plt.ylabel("V(t)")
plt.title("Approximation of V")
plt.draw()

# Compile the next example (conservation of energy in control volume)
comp("CtrlFlowSystem")

# Allocate a parser and load the xml
ocp = FlatOCP('CtrlFlowSystem.xml')

# Set options
ocp.setOption("make_explicit",True)

# Parse the model
ocp.init()

# Print the ocp
print ocp

# The problem has no differential states, so instead of integrating, we just solve for mdot...


plt.show()












