#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             K.U. Leuven. All rights reserved.
#     Copyright (C) 2011-2014 Greg Horn
#
#     CasADi is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 3 of the License, or (at your option) any later version.
#
#     CasADi is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public
#     License along with CasADi; if not, write to the Free Software
#     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
#
# -*- coding: utf-8 -*-
import os
import sys
import numpy as NP
from numpy import *
import matplotlib.pyplot as plt
import zipfile
import time
import shutil

try:
  # JModelica
  from pymodelica import compile_jmu
  from pyjmi import JMUModel
  import pymodelica
  use_precompiled = False
except:
  print "No jmodelica installation, falling back to precompiled XML-files"
  use_precompiled = True

# CasADi
from casadi import *

# Matplotlib interactive mode
#plt.ion()

# Compile Modelica code to XML
def comp(name):
  curr_dir = os.path.dirname(os.path.abspath(__file__))
  if use_precompiled:
    shutil.copy(curr_dir + '/precompiled_' + name + '.xml', name + '.xml')
  else:
    jmu_name = compile_jmu(name, curr_dir+"/thermodynamics_example.mo",'modelica','ipopt',{'generate_xml_equations':True, 'generate_fmi_me_xml':False})
    modname = name.replace('.','_')
    sfile = zipfile.ZipFile(curr_dir+'/'+modname+'.jmu','r')
    mfile = sfile.extract('modelDescription.xml','.')
    os.remove(modname+'.jmu')
    os.rename('modelDescription.xml',modname+'.xml')

# Compile the simplemost example (conservation of mass in control volume)
comp("BasicVolumeMassConservation")

# Read a model from XML
ocp = SymbolicOCP()
ocp.parseFMI('BasicVolumeMassConservation.xml')

# Transform into an explicit ODE
ocp.makeExplicit()

# Inputs to the integrator
dae_fcn_in = daeIn(
  t = ocp.t,
  x = ocp.x,
  p = vertcat((ocp.pi,ocp.p))
)

# Create an integrator
dae = SXFunction(dae_fcn_in,daeOut(ode=ocp.ode(ocp.x)))
integrator = Integrator("cvodes", dae)

# Output function
output_fcn_out = [ocp.beq("m"),ocp.beq("P")]
output_fcn_in = daeIn(
  t=ocp.t,
  x = ocp.x,
  z = ocp.z,
  p = vertcat((ocp.pi,ocp.p,ocp.u))
)
output_fcn = SXFunction(output_fcn_in,output_fcn_out)

# Create a simulator
grid = NP.linspace(0,1,100)
simulator = Simulator(integrator,output_fcn,grid)
simulator.init()

# Pass initial conditions
x0 = ocp.start(ocp.x)
simulator.setInput(x0,"x0")

# Simulate
simulator.evaluate()
integrator.printStats()

# Plot
plt.figure(1)
plt.subplot(1,2,1)
plt.plot(grid,simulator.getOutput().T)
plt.xlabel("t")
plt.ylabel("m(t)")
plt.title("c.f. Fritzson figure 15-6 (left)")

plt.subplot(1,2,2)
plt.plot(grid,simulator.getOutput(1).T)
plt.xlabel("t")
plt.ylabel("P(t)")
plt.title("c.f. Fritzson figure 15-6 (right)")
plt.draw()

# Compile the next example (conservation of energy in control volume)
comp("BasicVolumeEnergyConservation")

# Allocate a parser and load the xml
ocp = SymbolicOCP()
ocp.parseFMI('BasicVolumeEnergyConservation.xml')

# Transform into an explicit ODE
ocp.makeExplicit()

# Eliminate the independent variables
ocp.eliminateIndependentParameters()

# Inputs to the integrator
dae_fcn_in = daeIn(
  t = ocp.t,
  x = ocp.x,
  p = vertcat((ocp.pi,ocp.p))
)

# Create an integrator
dae = SXFunction(dae_fcn_in,daeOut(ode=ocp.ode(ocp.x)))
integrator = Integrator("cvodes", dae)

# Output function
output_fcn_out = [ocp.beq("T")]
output_fcn_in = daeIn(
  t=ocp.t,
  x = ocp.x,
  z = ocp.z,
  p = vertcat((ocp.pi,ocp.p,ocp.u))
)
output_fcn = SXFunction(output_fcn_in,output_fcn_out)

# Create a simulator
grid = NP.linspace(0,10,100)
simulator = Simulator(integrator,output_fcn,grid)
simulator.init()

# Pass initial conditions
x0 = ocp.start(ocp.x)
simulator.setInput(x0,"x0")

# Simulate
simulator.evaluate()
integrator.printStats()

# Plot
plt.figure(2)
plt.plot(grid,simulator.getOutput().T)
plt.xlabel("t")
plt.ylabel("T(t)")
plt.title("c.f. Fritzson figure 15-9")
plt.draw()

# Compile the next example (Heat transfer and work)
comp("BasicVolumeTest")

# Allocate a parser and load the xml
ocp = SymbolicOCP()
ocp.parseFMI('BasicVolumeTest.xml')

# Eliminate the dependent variables
ocp.eliminateIndependentParameters()

# Transform into an explicit ODE
ocp.makeExplicit()

# Inputs to the integrator
dae_fcn_in = daeIn(
  t = ocp.t,
  x = ocp.x,
  p = vertcat((ocp.pi,ocp.p))
)

# Create an integrator
dae = SXFunction(dae_fcn_in,daeOut(ode=ocp.ode(ocp.x)))
integrator = Integrator("cvodes", dae)

# Output function
output_fcn_out = [ocp.beq("T"),ocp.beq("U"),ocp.beq("V")]
output_fcn_in = daeIn(
  t=ocp.t,
  x = ocp.x,
  z = ocp.z,
  p = vertcat((ocp.pi,ocp.p,ocp.u))
)
output_fcn = SXFunction(output_fcn_in,output_fcn_out)

# Create a simulator
grid = NP.linspace(0,2,100)
simulator = Simulator(integrator,output_fcn,grid)
simulator.init()

# Pass initial conditions
x0 = ocp.start(ocp.x)
simulator.setInput(x0,"x0")

# Simulate
simulator.evaluate()
integrator.printStats()

# Plot
plt.figure(3)
p1, = plt.plot(grid,simulator.getOutput(0).T)
p2, = plt.plot(grid,simulator.getOutput(1).T)
plt.xlabel("t")
plt.ylabel("T(t)")
plt.legend([p2, p1], ["T", "U"])
plt.title("c.f. Fritzson figure 15-14")

plt.figure(4)
plt.plot(grid,simulator.getOutput(2).T)
plt.xlabel("t")
plt.ylabel("V(t)")
plt.title("Approximation of V")
plt.draw()

# Compile the next example (conservation of energy in control volume)
comp("CtrlFlowSystem")

# Allocate a parser and load the xml
ocp = SymbolicOCP()
ocp.parseFMI('CtrlFlowSystem.xml')

# Transform into a semi-explicit ODE
ocp.makeSemiExplicit()

# Print the ocp
print ocp

# The problem has no differential states, so instead of integrating, we just solve for mdot...


plt.show()
