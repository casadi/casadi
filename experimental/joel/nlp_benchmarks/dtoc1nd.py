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
import numpy as NP
import matplotlib.pyplot as plt

# This is a CasADi version of dtoc1nd.mod from the cute test collection, original by Hande Y. Benson
# original
#n = 50 
#nx = 5
#ny = 10

# original
n = 1000
nx = 2
ny = 4

mu = 1.0

b = SX.nan(ny+1,nx+1)
c = SX.nan(ny+1,nx+1)
for i in range(1,ny+1): 
  for j in range(1,nx+1): 
    b[i,j] = float(i-j)/(nx+ny)
    c[i,j] = float(i+j)*mu/(nx+ny)


x = ssym("x",n-1, nx)
y = ssym("y",n, ny)

x_lb = -inf * NP.ones(x.shape)
x_ub =  inf * NP.ones(x.shape)
x_guess = NP.zeros(x.shape)

y_lb = -inf * NP.ones(y.shape)
y_ub =  inf * NP.ones(y.shape)
y_guess = NP.zeros(y.shape)

y_lb[0,:] = 0.0
y_ub[0,:] = 0.0

# vectorize
def NPvec(x):
  return NP.reshape(x,NP.prod(x.shape),1)

# All variables
v = vertcat([vec(x),vec(y)])
v_lb = NP.concatenate([NPvec(x_lb),NPvec(y_lb)])
v_ub = NP.concatenate([NPvec(x_ub),NPvec(y_ub)])
v_guess = NP.concatenate([NPvec(x_guess),NPvec(y_guess)])

print "vars ok"

# Objective function
f = 0
for t in range(1,n):
  for i in range(1,nx+1):
    f += (x[t-1,i-1] + 0.5)**4
for t in range(1,n+1):
  for i in range(1,ny+1):
    f += (y[t-1,i-1] + 0.25)**4

ffcn = SXFunction([v],[f])
print "obj ok"

# Constraint functions
cons1 = []
cons1_lb = []
cons1_ub = []
for t in range(1,n):
  ct = 0
  for k in range(ny*nx):
    ct += c[(k/nx)+1,k-nx*(k/nx)+1]*y[t-1,(k/nx)+1-1]*x[t-1,k-nx*(k/nx)+1-1]+0.5*y[t-1,1-1] + 0.25*y[t-1,2-1] - y[t+1-1,1-1] 
  for i in range(1,nx+1):
    ct += b[1,i]*x[t-1,i-1]
  cons1.append(ct)
  cons1_lb.append(0)
  cons1_ub.append(0)

print "cons1"

cons2 = []
cons2_lb = []
cons2_ub = []
for t in range(1,n):
  for j in range(2,ny):
    ctj = 0
    for k in range(ny*nx):
      ctj += c[(k/nx)+1,k-nx*(k/nx)+1]*y[t-1,(k/nx)+1-1]*x[t-1,k-nx*(k/nx)+1-1]-y[t+1-1,j-1] + 0.5*y[t-1,j-1] - 0.25*y[t-1,j-1-1] + 0.25*y[t-1,j+1-1] 
    for i in range(1,nx+1):
      ctj += b[j,i]*x[t-1,i-1]
    cons2.append(ctj)
    cons2_lb.append(0)
    cons2_ub.append(0)

print "cons2"

cons3 = []
cons3_lb = []
cons3_ub = []
for t in range(1,n):
  ct = 0
  for k in range(ny*nx):
    ct += c[(k/nx)+1,k-nx*(k/nx)+1]*y[t-1,(k/nx)+1-1]*x[t-1,k-nx*(k/nx)+1-1]+0.5*y[t-1,ny-1] - 0.25*y[t-1,ny-1-1] - y[t+1-1,ny-1]
  for k in range(1,nx+1):
    ct += b[ny,i]*x[t-1,i-1]
  cons3.append(ct)
  cons3_lb.append(0)
  cons3_ub.append(0)

print "cons3"

g = vertcat([cons1,cons2,cons3])
g_lb = NP.concatenate([cons1_lb,cons2_lb,cons3_lb])
g_ub = NP.concatenate([cons1_ub,cons2_ub,cons3_ub])

gfcn = SXFunction([v],[g])

# NLP solver
print "ipopt"
nlp_solver = IpoptSolver(ffcn,gfcn)
#nlp_solver.setOption("verbose",True)
#nlp_solver.setOption("max_iter",10)
nlp_solver.setOption("generate_hessian",True)
nlp_solver.init()
nlp_solver.setInput(v_guess, "x0")
nlp_solver.setInput(v_lb, "lbx")
nlp_solver.setInput(v_ub, "ubx")
nlp_solver.setInput(g_lb, "lbg")
nlp_solver.setInput(g_ub, "ubg")

nlp_solver.solve()




