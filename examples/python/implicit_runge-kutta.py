#
#     This file is part of CasADi.
# 
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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
import numpy as NP
import matplotlib.pyplot as plt

"""
Demonstration on how to construct a fixed-step implicit Runge-Kutta integrator
@author: Joel Andersson, K.U. Leuven 2013
"""

tf = 10.0  # End time
ni = 100     # Number of finite elements

# Size of the finite elements
h = tf/ni

# Dimensions
nx = 3
nu = 1

# Declare variables
x  = ssym("x",nx)  # state
u  = ssym("u",nu)  # control

# ODE right hand side function
rhs = vertcat([(1 - x[1]*x[1])*x[0] - x[1] + u, \
               x[0], \
               x[0]*x[0] + x[1]*x[1] + u*u])
f = SXFunction(daeIn(x=x,p=u),daeOut(ode=rhs))
f.init()

# Legendre collocation points
legendre_points1 = [0,0.500000]
legendre_points2 = [0,0.211325,0.788675]
legendre_points3 = [0,0.112702,0.500000,0.887298]
legendre_points4 = [0,0.069432,0.330009,0.669991,0.930568]
legendre_points5 = [0,0.046910,0.230765,0.500000,0.769235,0.953090]

# Radau collocation points
radau_points1 = [0,1.000000]
radau_points2 = [0,0.333333,1.000000]
radau_points3 = [0,0.155051,0.644949,1.000000]
radau_points4 = [0,0.088588,0.409467,0.787659,1.000000]
radau_points5 = [0,0.057104,0.276843,0.583590,0.860240,1.000000]

# Choose collocation points
tau_root = legendre_points4

# Degree of interpolating polynomial
d = len(tau_root)-1

# Coefficients of the collocation equation
C = NP.zeros((d+1,d+1))

# Coefficients of the continuity equation
D = NP.zeros(d+1)

# Dimensionless time inside one control interval
tau = ssym("tau")
  
# For all collocation points
for j in range(d+1):
  # Construct Lagrange polynomials to get the polynomial basis at the collocation point
  L = 1
  for r in range(d+1):
    if r != j:
      L *= (tau-tau_root[r])/(tau_root[j]-tau_root[r])
  lfcn = SXFunction([tau],[L])
  lfcn.init()
  
  # Evaluate the polynomial at the final time to get the coefficients of the continuity equation
  lfcn.setInput(1.0)
  lfcn.evaluate()
  D[j] = lfcn.output()

  # Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
  for r in range(d+1):
    lfcn.setInput(tau_root[r])
    lfcn.setFwdSeed(1.0)
    lfcn.evaluate(1,0)
    C[j,r] = lfcn.fwdSens()

# Total number of variables for one finite element
x0_elem = MX("x0_elem",nx) # given
u_elem = MX("u_elem",nu)  # given
x_elem = MX("x_elem",d*nx) # sought

# Get the state at each collocation point
x_elem_split = [x0_elem if r==0 else x_elem[(r-1)*nx:r*nx] for r in range(d+1)]

# Constraint function for the NLP
g_elem = []

# For all collocation points
for j in range(1,d+1):
        
  # Get an expression for the state derivative at the collocation point
  xp_j = 0
  for r in range (d+1):
    xp_j += C[r,j]*x_elem_split[r]
      
  # Add collocation equations to the NLP
  [f_j] = daeOut(f.call(daeIn(x=x_elem_split[j],p=u_elem)),"ode")
  g_elem.append(h*f_j - xp_j)

# Concatenate constraints
g_elem = vertcat(g_elem)

# Create the nonlinear system of equations
gfcn_elem = MXFunction([x_elem,x0_elem,u_elem],[g_elem])

# Get an expression for the solution of the system of equations
ifcn_elem = NewtonImplicitSolver(gfcn_elem)
ifcn_elem.setOption("linear_solver",CSparse)
ifcn_elem.init()

# Get an expression for the state at the end of the finite element
[x_elem] = ifcn_elem.call([x0_elem,u_elem])
x_elem_split = [x0_elem if r==0 else x_elem[(r-1)*nx:r*nx] for r in range(d+1)]
xf_elem = 0
for r in range(d+1):
  xf_elem += D[r]*x_elem_split[r]
  
# Get the discrete time dynamics
ffcn_elem = MXFunction([x0_elem,u_elem],[xf_elem])
ffcn_elem.init()

# Integrate over one control interval
x_shot = x0_elem
for i in range(ni):
  [x_shot] = ffcn_elem.call([x_shot,u_elem])

# Fixed-step integrator
irk_integrator = MXFunction(integratorIn(x0=x0_elem,p=u_elem),integratorOut(xf=x_shot))
irk_integrator.setOption("name","irk_integrator")
irk_integrator.init()

# Create a convensional integrator for reference
ref_integrator = CVodesIntegrator(f)
ref_integrator.setOption("name","ref_integrator")
ref_integrator.setOption("tf",tf)
ref_integrator.init()

# Test values
x0_val  = NP.array([0,1,0])
u_val = 0.2

# Make sure that both integrators give consistent results
for integrator in (irk_integrator,ref_integrator):
  print "Testing ", repr(integrator)

  # Pass arguments
  integrator.setInput(x0_val,"x0")
  integrator.setInput(u_val,"p")
  
  # Integrate
  integrator.evaluate()

  # Get the results
  print "xf = ", integrator.output("xf")

# Print stats for CVodes
print "------------------"
print "CVodes statistics:"
ref_integrator.printStats()

