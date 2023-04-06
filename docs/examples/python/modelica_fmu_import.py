#
#     MIT No Attribution
#
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
#
#     Permission is hereby granted, free of charge, to any person obtaining a copy of this
#     software and associated documentation files (the "Software"), to deal in the Software
#     without restriction, including without limitation the rights to use, copy, modify,
#     merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
#     permit persons to whom the Software is furnished to do so.
#
#     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#     INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#     PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# -*- coding: utf-8 -*-
from casadi import *

# Consider the following simple Modelica file
with open('rocket.mo', 'r') as f: print(f.read())

# Different tools exist to compile Modelica models like these into FMUs.
# The following code shows how we could do it using the open-source tool 
# OpenModelica/OMPython, cf. the OpenModelica User's Guide.
from OMPython import OMCSessionZMQ
omc = OMCSessionZMQ()
if omc.loadFile('rocket.mo').startswith('false'):
  raise Exception('Modelica compilation failed: {}'.format(omc.sendExpression('getErrorString()')))
fmu_file = omc.sendExpression('translateModelFMU(rocket)')
flag = omc.sendExpression('getErrorString()')
if not fmu_file.endswith('.fmu'): raise Exception('FMU generation failed: {}'.format(flag))
print("translateModelFMU warnings:\n{}".format(flag))

# Regardless of how the FMU was obtained, we must unzip it before CasADi can process it
import zipfile
with zipfile.ZipFile('rocket.fmu', 'r') as zip_ref: zip_ref.extractall('rocket_fmu_unzipped')

# Create a CasADi/DaeBuilder instance from the unzipped FMU
dae = DaeBuilder('rocket', 'rocket_fmu_unzipped')
dae.disp(True)

# Get state vector, initial conditions, bounds
x = dae.x()
lbx = dae.min(x)
ubx = dae.max(x)
x0 = dae.start(x)
print('x: ', x)
print('lbx: ', lbx)
print('ubx: ', ubx)
print('x0: ', x0)

# Get control vector, initial conditions, bounds
u = dae.u()
lbu = dae.min(u)
ubu = dae.max(u)
u0 = dae.start(u)
print('u: ', u)
print('lbu: ', lbu)
print('ubu: ', ubu)
print('u0: ', u0)

# Let's create a CasADi function for evaluating the ODE right-hand-side.
# We only need to specify the expressions that are varying:
f = dae.create('f', ['x', 'u'], ['ode'])

# This is a standard CasADi function that can be embedded into CasADi expressions
print(f)

# We can evaluate it and calculate derivatives as with other CasADi functions
ode0 = f(x0, u0)
print('ode0: ', ode0)

# Analytic first derivatives with sparsities are also available
jac = f.factory('jac_f', ['x', 'u'], ['jac:ode:x', 'jac:ode:u'])
jac_ode_x0, jac_ode_u0 = jac(x0, u0)
print('jac_ode_x0: ', jac_ode_x0)
print('jac_ode_u0: ', jac_ode_u0)

# Code generation for the functions are currently work in progress
f.generate('fmu_codegen')
