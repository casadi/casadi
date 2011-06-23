# -*- coding: utf-8 -*-
import os
import sys
import numpy as NP
from numpy import *
import matplotlib.pyplot as plt
import zipfile
import time

# JModelica
from jmodelica.jmi import compile_jmu
from jmodelica.jmi import JMUModel
import jmodelica

# CasADi
from casadi import *

# Compile Modelica code to XML
curr_dir = os.path.dirname(os.path.abspath(__file__))
jmu_name = compile_jmu("BasicVolume", curr_dir+"/fritzson_15_2.mo",'modelica','ipopt',{'generate_xml_equations':True, 'generate_fmi_me_xml':False})
sfile = zipfile.ZipFile(curr_dir+'/BasicVolume.jmu','r')
mfile = sfile.extract('modelDescription.xml','.')
os.remove('BasicVolume.jmu')
os.rename('modelDescription.xml','BasicVolume.xml')

# Allocate a parser and load the xml
parser = FMIParser('BasicVolume.xml')

# Obtain the symbolic representation of the OCP
ocp = parser.parse()

# Create functions
ocp.createFunctions()

# Create an integrator
integrator = CVodesIntegrator(ocp.oderhs_)

# Get the variables
vv = ocp.getVariables()

# Output function
output_expression = [[ocp.getExplicit(vv.m.var())],[ocp.getExplicit(vv.P.var())]]
output_fcn_in = DAE_NUM_IN * [[]]
output_fcn_in[DAE_T] = [ocp.t_]
output_fcn_in[DAE_Y] = var(ocp.x_)
output_fcn_in[DAE_YDOT] = der(ocp.x_)
output_fcn_in[DAE_P] = var(ocp.p_)
output_fcn = SXFunction(output_fcn_in,output_expression)

# Create a simulator
grid = NP.linspace(0,1,100)
simulator = Simulator(integrator,output_fcn,grid)
simulator.init()

# Pass initial conditions
x0 = getStart(ocp.x_)
simulator.setInput(x0,INTEGRATOR_X0)

# Simulate
simulator.evaluate()

# Plot
plt.subplot(1,2,1)
plt.plot(grid,simulator.output())
plt.xlabel("t")
plt.ylabel("m(t)")

plt.subplot(1,2,2)
plt.plot(grid,simulator.output(1))
plt.xlabel("t")
plt.ylabel("P(t)")

plt.show()










