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

# This is a CasADi version of optmass.mod from the cute test collection, original by Hande Y. Benson
#   
#   M. Gawande and J. Dunn,
#   "A Projected Newton Method in a Cartesian Product of Balls",
#   JOTA 59(1): 59-69, 1988.

n = 1000 # original
#n = 10000 # website
#n = 1000

speed = 0.01
pen = 0.335

x = ssym("x", 2, n+2)
x_lb = -inf * NP.ones(x.shape)
x_ub =  inf * NP.ones(x.shape)
x_guess = NP.zeros(x.shape)

x_lb[:,0] = x_ub[:,0] = 0

v = ssym("v", 2, n+2)
v_lb = -inf * NP.ones(v.shape)
v_ub =  inf * NP.ones(v.shape)
v_guess = NP.zeros(v.shape)

v_lb[0,0] = v_ub[0,0] = speed
v_lb[1,0] = v_ub[1,0] = 0.0

f = ssym("f", 2, n+1)
f_lb = -inf * NP.ones(f.shape)
f_ub =  inf * NP.ones(f.shape)
f_guess = NP.zeros(f.shape)

# vectorize
def NPvec(x):
  return NP.reshape(x,NP.prod(x.shape),1)

# All variables
vv = vertcat([vec(x),vec(v),vec(f)])
vv_lb = NP.concatenate([NPvec(x_lb),NPvec(v_lb),NPvec(f_lb)])
vv_ub = NP.concatenate([NPvec(x_ub),NPvec(v_ub),NPvec(f_ub)])
vv_guess = NP.concatenate([NPvec(x_guess),NPvec(v_guess),NPvec(f_guess)])

# Objective function
obj = pen*(v[1-1,n+1-1]**2+v[2-1,n+1-1]**2) - (x[1-1,n+1-1]**2+x[2-1,n+1-1]**2)
ffcn = SXFunction([vv],[obj])

# constraints
g = []
g_lb = []
g_ub = []

#cons1
for i in range(1,n+2):
  for j in [0,1]:
    g.append(x[j,i] - x[j,i-1] - v[j,i-1]/n - f[j,i-1]/(2*n**2))
    g_lb.append(0)
    g_ub.append(0)

#cons2
for i in range(1,n+2):
  for j in [0,1]:
    g.append(v[j,i] - v[j,i-1] - f[j,i-1]/n)
    g_lb.append(0)
    g_ub.append(0)
    
#cons3
for i in range(n+1):
  g.append(f[0,i]**2 + f[1,i]**2)
  g_lb.append(-inf)
  g_ub.append(1)
  
g = vertcat(g)
g_lb = NP.array(g_lb)
g_ub = NP.array(g_ub)
gfcn = SXFunction([vv],[g])

# NLP solver
nlp_solver = IpoptSolver(ffcn,gfcn)
#nlp_solver.setOption("max_iter",10)
nlp_solver.setOption("generate_hessian",True)
nlp_solver.setOption("linear_solver","ma57")
nlp_solver.init()
nlp_solver.setInput(vv_guess, "x0")
nlp_solver.setInput(vv_lb, "lbx")
nlp_solver.setInput(vv_ub, "ubx")
nlp_solver.setInput(g_lb, "lbg")
nlp_solver.setInput(g_ub, "ubg")

nlp_solver.solve()




