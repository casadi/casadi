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

# This is a CasADi version of clnlbeam.mod from the cute test collection, original by Hande Y. Benson
#   
#   H. Maurer and H.D. Mittelman,
#   "The non-linear beam via optimal control with bound state variables",
#   Optimal Control Applications and Methods 12, pp. 19-31, 1991.

#ni = 500 # original
#ni = 20000 # web
ni = 500 # cuter
alpha = 350.0
h = 1./ni

t = ssym("t", ni+1)
t_lb = -1.0 * NP.ones(ni+1)
t_ub =  1.0 * NP.ones(ni+1)
t_guess =  NP.array(list(0.05*cos(i*h) for i in range(ni+1)))

t_lb[0] = 0;  t_ub[0] = 0
t_lb[ni] = 0; t_ub[ni] = 0

x = ssym("x", ni+1)
x_lb = -0.05 * NP.ones(ni+1)
x_ub =  0.05 * NP.ones(ni+1)
x_guess =  NP.array(list(0.05*cos(i*h) for i in range(ni+1)))

x_lb[0] = 0;  x_ub[0] = 0
x_lb[ni] = 0; x_ub[ni] = 0

u = ssym("u", ni+1)
u_lb = -inf * NP.ones(ni+1)
u_ub =  inf * NP.ones(ni+1)
u_guess =  0.0 * NP.ones(ni+1)

# All variables
v = vertcat([t,x,u])
v_lb = NP.concatenate([t_lb,x_lb,u_lb])
v_ub = NP.concatenate([t_ub,x_ub,u_ub])
v_guess = NP.concatenate([t_guess,x_guess,u_guess])

# Make h, alpha symbolic once
h = SX(h)
alpha = SX(alpha)

# Objective function
f = 0
for i in range(ni):
  f += 0.5*h*(u[i+1]**2 + u[i]**2) + 0.5*alpha*h*(cos(t[i+1]) + cos(t[i]))

ffcn = SXFunction([v],[f])

# Constraint function
g = []
g_lb = []
g_ub = []
for i in range(ni):
  g.append(x[i+1] - x[i] - 0.5*h*(sin(t[i+1]) + sin(t[i])))
  g_lb.append(0)
  g_ub.append(0)

for i in range(ni):
  g.append(t[i+1] - t[i] - 0.5*h*u[i+1] - 0.5*h*u[i])
  g_lb.append(0)
  g_ub.append(0)
g = vertcat(g)
g_lb = NP.array(g_lb)
g_ub = NP.array(g_ub)

gfcn = SXFunction([v],[g])

# NLP solver
nlp_solver = IpoptSolver(ffcn,gfcn)
#nlp_solver.setOption("max_iter",10)
nlp_solver.setOption("generate_hessian",True)
nlp_solver.setOption("linear_solver","ma57")
nlp_solver.init()
nlp_solver.setInput(v_guess, "x0")
nlp_solver.setInput(v_lb, "lbx")
nlp_solver.setInput(v_ub, "ubx")
nlp_solver.setInput(g_lb, "lbg")
nlp_solver.setInput(g_ub, "ubg")

nlp_solver.solve()




