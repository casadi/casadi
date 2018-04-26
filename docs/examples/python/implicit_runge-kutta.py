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
from casadi import *
import numpy as N
import matplotlib.pyplot as plt

'''
Demonstration on how to construct a fixed-step implicit Runge-Kutta integrator
@author: Joel Andersson, K.U. Leuven 2013
'''

# End time
tf = 10.0

# Dimensions
nx = 3
np = 1

# Declare variables
x  = SX.sym('x', nx)  # state
p  = SX.sym('u', np)  # control

# ODE right hand side function
ode = vertcat((1 - x[1]*x[1])*x[0] - x[1] + p, \
              x[0], \
              x[0]*x[0] + x[1]*x[1] + p*p)
dae = {'x':x, 'p':p, 'ode':ode}
f = Function('f', [x, p], [ode])

# Number of finite elements
n = 100

# Size of the finite elements
h = tf/n

# Degree of interpolating polynomial
d = 4

# Choose collocation points
tau_root = [0] + collocation_points(d, 'legendre')

# Coefficients of the collocation equation
C = N.zeros((d+1,d+1))

# Coefficients of the continuity equation
D = N.zeros(d+1)

# Dimensionless time inside one control interval
tau = SX.sym('tau')

# For all collocation points
for j in range(d+1):
  # Construct Lagrange polynomials to get the polynomial basis at the collocation point
  L = 1
  for r in range(d+1):
    if r != j:
      L *= (tau-tau_root[r])/(tau_root[j]-tau_root[r])

  # Evaluate the polynomial at the final time to get the coefficients of the continuity equation
  lfcn = Function('lfcn', [tau], [L])
  D[j] = lfcn(1.0)

  # Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
  tfcn = Function('tfcn', [tau], [tangent(L,tau)])
  for r in range(d+1): C[j,r] = tfcn(tau_root[r])

# Total number of variables for one finite element
X0 = MX.sym('X0',nx)
P  = MX.sym('P',np)
V = MX.sym('V',d*nx)

# Get the state at each collocation point
X = [X0] + vertsplit(V,[r*nx for r in range(d+1)])

# Get the collocation quations (that define V)
V_eq = []
for j in range(1,d+1):
  # Expression for the state derivative at the collocation point
  xp_j = 0
  for r in range (d+1):
    xp_j += C[r,j]*X[r]

  # Append collocation equations
  f_j = f(X[j],P)
  V_eq.append(h*f_j - xp_j)

# Concatenate constraints
V_eq = vertcat(*V_eq)

# Root-finding function, implicitly defines V as a function of X0 and P
vfcn = Function('vfcn', [V, X0, P], [V_eq])

# Convert to SX to decrease overhead
vfcn_sx = vfcn.expand()

# Create a implicit function instance to solve the system of equations
ifcn = rootfinder('ifcn', 'newton', vfcn_sx)
V = ifcn(MX(),X0,P)
X = [X0 if r==0 else V[(r-1)*nx:r*nx] for r in range(d+1)]

# Get an expression for the state at the end of the finie element
XF = 0
for r in range(d+1):
  XF += D[r]*X[r]

# Get the discrete time dynamics
F = Function('F', [X0,P],[XF])

# Do this iteratively for all finite elements
X = X0
for i in range(n):
  X = F(X,P)

# Fixed-step integrator
irk_integrator = Function('irk_integrator', {'x0':X0, 'p':P, 'xf':X},
                          integrator_in(), integrator_out())

# Create a convensional integrator for reference
ref_integrator = integrator('ref_integrator', 'cvodes', dae, {'tf':tf})

# Test values
x0_val  = N.array([0,1,0])
p_val = 0.2

# Make sure that both integrators give consistent results
for F in (irk_integrator,ref_integrator):
  print('-------')
  print('Testing ' + F.name())
  print('-------')

  # Generate a new function that calculates forward and reverse directional derivatives
  dF = F.factory('dF', ['x0', 'p', 'fwd:x0', 'fwd:p', 'adj:xf'],
                       ['xf', 'fwd:xf', 'adj:x0', 'adj:p']);
  arg = {}

  # Pass arguments
  arg['x0'] = x0_val
  arg['p'] = p_val

  # Forward sensitivity analysis, first direction: seed p and x0[0]
  arg['fwd_x0'] = [1,0,0]
  arg['fwd_p'] = 1

  # Adjoint sensitivity analysis, seed xf[2]
  arg['adj_xf'] = [0,0,1]

  # Integrate
  res = dF(**arg)

  # Get the nondifferentiated results
  print('%30s = %s' % ('xf', res['xf']))

  # Get the forward sensitivities
  print('%30s = %s' % ('d(xf)/d(p)+d(xf)/d(x0[0])', res['fwd_xf']))

  # Get the adjoint sensitivities
  print('%30s = %s' % ('d(xf[2])/d(x0)', res['adj_x0']))
  print('%30s = %s' % ('d(xf[2])/d(p)', res['adj_p']))
