#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
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
