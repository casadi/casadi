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
import casadi
from jmodelica.compiler import OptimicaCompiler
from jmodelica.compiler import ModelicaCompiler

import matplotlib.pyplot as plt

import numpy as NP

oc = OptimicaCompiler()
oc.set_boolean_option('generate_xml_equations',True)
#oc.compile_model("VDP_Opt","VDP.mo")

#mc = ModelicaCompiler()
#mc.compile_model("VDP","VDP2.mo")

# Allocate a parser and load the xml
#parser = casadi.FMIParser('VDP_Opt.xml')

parser = casadi.FMIParser('/home/janderss/dev/OPTICON-SOFTWARE/swig_interface/test/xml_files/VDP_pack_VDP_Opt.xml')

# Dump representation to screen
print "XML representation"
print parser

# Obtain the symbolic representation of the OCP
ocp = parser.parse()

# Sort the variables according to type
ocp.sortVariables()

# Print the ocp to screen
print ocp

# Make explicit
#ocp.makeExplicit()

F_inputs = []
F_inputs += ocp.xdot
F_inputs += ocp.x
F_inputs += ocp.u
F_inputs += ocp.xa
F_inputs += ocp.p
F_inputs += ocp.t

# The DAE function
F = casadi.SXFunction([F_inputs],[ocp.dyneq])

# Number of elements
N = 100
tf = 5.

x_init = [0,1,0]

# Build variable data structure

n_x = len(ocp.x)
n_u = len(ocp.u)
n_w = len(ocp.xa)

dxs = []
xs = []
us = []
ws = []

# Add initial point for the states
x0 = []
for i in range(n_x):
   x0.append(casadi.SX('x_0_'+str(i+1)))
xs.append(x0)

for i in range(N):
   dxi = []
   xi = []
   ui = []
   wi = []
   for j in range(n_x):
       dxi.append(casadi.SX('dx_'+str(i+1)+'_'+str(j+1)))
   dxs.append(dxi)
   for j in range(n_x):
       xi.append(casadi.SX('x_'+str(i+1)+'_'+str(j+1)))
   xs.append(xi)
   for j in range(n_u):
       ui.append(casadi.SX('u_'+str(i+1)+'_'+str(j+1)))
   us.append(ui)
   for j in range(n_w):
       wi.append(casadi.SX('w_'+str(i+1)+'_'+str(j+1)))
   ws.append(wi)

xx = []
xx += xs[0]
for i in range(N):
   xx += dxs[i]
   xx += xs[i+1]
   xx += us[i]
   xx += ws[i]

# Equality constraints
g = []

for i in range(n_x):
   g.append(xs[0][i] - x_init[i])

for i in range(N):
   z = []
   z += dxs[i]
   z += xs[i+1]
   z += ws[i]
   z += us[i]
   z += ocp.p
   z += ocp.t
   g += list(F.eval([z])[0])
   for j in range(n_x):
       g.append(dxs[i][j] - (xs[i+1][j] - xs[i][j])/(tf/N))

cost = xs[N][0]

# Equality constraint residual
g_res = casadi.SXFunction([xx],[g])

# Objective function
obj = casadi.SXFunction([xx], [[cost]])

# Hessian
sigma = casadi.SX('sigma')

lam = []
Lag = sigma*cost
for i in range(len(g)):
   lam.append(casadi.SX('lambda_' + str(i)))
   Lag = Lag + g[i]*lam[i]

Lag_fcn = casadi.SXFunction([xx, lam, [sigma]],[[Lag]])
H = Lag_fcn.hess()

#H_fcn = casadi.SXFunction([xx, lam, [sigma]],[H])
H_fcn = Lag_fcn.hessian(0,0)

# Create Solver
solver = casadi.IpoptSolver(obj,g_res,H_fcn)

# Set options
solver.setOption("tol",1e-10)
#solver.setOption("derivative_test",'second-order')
#solver.setOption("max_iter",30)
#solver.setOption("hessian_approximation","limited-memory")

# Initialize
solver.init()

# Initial condition
x_initial_guess = len(xx) * [0]
solver.setInput(x_initial_guess,casadi.NLP_SOLVER_X0)

# Bounds on x
xx_lb = len(xx)*[-100]
xx_ub = len(xx)*[100]
solver.setInput(xx_lb,casadi.NLP_SOLVER_LBX)
solver.setInput(xx_ub,casadi.NLP_SOLVER_UBX)

# Bounds on the constraints
lubg = len(g)*[0]
solver.setInput(lubg,casadi.NLP_SOLVER_LBG)
solver.setInput(lubg,casadi.NLP_SOLVER_UBG)

# Solve the problem
solver.solve()

xx_opt = solver.getOutput(casadi.NLP_SOLVER_X)

# Retrieve output

dx_opt = NP.zeros((N,n_x))
x_opt = NP.zeros((N+1,n_x))
w_opt = NP.zeros((N,n_w))
u_opt = NP.zeros((N,n_u))

t = NP.linspace(0,tf,N+1)

x_opt[0,:] = xx_opt[0:3]

n_vars_per_el = n_x*2 + n_w + n_u

for i in range(N):
   offs = 3 + n_vars_per_el*i
   dx_opt[i,:] = xx_opt[offs:offs+n_x]
   offs = 3 + n_vars_per_el*i + n_x
   x_opt[i+1,:] = xx_opt[offs:offs+n_x]
   offs = 3 + n_vars_per_el*i + n_x + n_x
   w_opt[i,:] = xx_opt[offs:offs+n_w]
   offs = 3 + n_vars_per_el*i + n_x + n_x + n_w
   u_opt[i,:] = xx_opt[offs:offs+n_u]


plt.figure(1)
plt.clf()
plt.subplot(2,1,1)
plt.plot(t,x_opt[:,1])
plt.hold(True)
plt.plot(t,x_opt[:,2])
plt.grid(True)
plt.subplot(2,1,2)
plt.plot(t[1:],u_opt[:,0])
plt.grid(True)
plt.show()
