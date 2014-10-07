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

# This is a CasADi version of corkscrw.mod from the cute test collection, original by Hande Y. Benson
#   
#   Ph. Toint

t=1000 #original
#t=5000 #web
#t = 500 # cuter
xt = 10.0
mass = 0.37
tol = 0.1

h = xt/t
w = xt*(t+1)/2.

fmax = xt/t

x = ssym("x", t)
xe = vertcat((0.,x))
x_lb =  0.0 * NP.ones(t)
x_ub =   xt * NP.ones(t)
x_guess =  NP.array(list(i*h for i in range(1,t+1)))

#x_lb[0] = 0.0
#x_ub[0] = 0.0

y = ssym("y", t)
ye = vertcat((0.,y))
y_lb = -inf * NP.ones(t)
y_ub =  inf * NP.ones(t)
y_guess = 0.0 * NP.ones(t)

#y_lb[0] = 0.0
#y_ub[0] = 0.0

z = ssym("z", t)
ze = vertcat((1.,z))
z_lb = -inf * NP.ones(t)
z_ub =  inf * NP.ones(t)
z_guess = 0.0 * NP.ones(t)

#z_lb[0] = 1.0
#z_ub[0] = 1.0

vx = ssym("vx", t-1)
vxe = vertcat((0.,vx,0.))
vx_lb = -inf * NP.ones(t-1)
vx_ub =  inf * NP.ones(t-1)
vx_guess = 1.0 * NP.ones(t-1)

#vx_lb[0] = 0.0
#vx_ub[0] = 0.0
#vx_lb[t-1] = 0.0
#vx_ub[t-1] = 0.0

vy = ssym("vy", t-1)
vye = vertcat((0.,vy,0.))
vy_lb = -inf * NP.ones(t-1)
vy_ub =  inf * NP.ones(t-1)
vy_guess = 0.0 * NP.ones(t-1)

#vy_lb[0] = 0.0
#vy_ub[0] = 0.0
#vy_lb[t-1] = 0.0
#vy_ub[t-1] = 0.0

vz = ssym("vz", t-1)
vze = vertcat((0.,vz,0.))
vz_lb = -inf * NP.ones(t-1)
vz_ub =  inf * NP.ones(t-1)
vz_guess = 0.0 * NP.ones(t-1)

#vz_lb[0] = 0.0
#vz_ub[0] = 0.0
#vz_lb[t-1] = 0.0
#vz_ub[t-1] = 0.0

ux = ssym("ux", t)
ux_lb = -fmax * NP.ones(t)
ux_ub =  fmax * NP.ones(t)
ux_guess = 0.0 * NP.ones(t)

uy = ssym("uy", t)
uy_lb = -fmax * NP.ones(t)
uy_ub =  fmax * NP.ones(t)
uy_guess = 0.0 * NP.ones(t)

uz = ssym("uz", t)
uz_lb = -fmax * NP.ones(t)
uz_ub =  fmax * NP.ones(t)
uz_guess = 0.0 * NP.ones(t)

# All variables
v = vertcat([x,y,z,vx,vy,vz,ux,uy,uz])
v_lb = NP.concatenate([x_lb,y_lb,z_lb,vx_lb,vy_lb,vz_lb,ux_lb,uy_lb,uz_lb])
v_ub = NP.concatenate([x_ub,y_ub,z_ub,vx_ub,vy_ub,vz_ub,ux_ub,uy_ub,uz_ub])
v_guess = NP.concatenate([x_guess,y_guess,z_guess,vx_guess,vy_guess,vz_guess,ux_guess,uy_guess,uz_guess])

# Make xt, h, mass symbolic once
xt = SX(xt)
mass = SX(mass)
tol = SX(tol)

# Objective function
f = 0
for i in range(1,t+1):
  f += (i*h/w)*(xe[i] - xt)**2

ffcn = SXFunction([v],[f])
h = SX(h)

# Constraint functions
acx = []
acx_lb = []
acx_ub = []
for i in range(1,t+1):
  acx.append(mass*(vxe[i]-vxe[i-1])/h - ux[i-1])
  acx_lb.append(0)
  acx_ub.append(0)
  
acy = []
acy_lb = []
acy_ub = []
for i in range(1,t+1):
  acy.append(mass*(vye[i]-vye[i-1])/h - uy[i-1])
  acy_lb.append(0)
  acy_ub.append(0)
  
acz = []
acz_lb = []
acz_ub = []
for i in range(1,t+1):
  acz.append(mass*(vze[i]-vze[i-1])/h - uz[i-1])
  acz_lb.append(0)
  acz_ub.append(0)

psx = []
psx_lb = []
psx_ub = []
for i in range(1,t+1):
  psx.append((xe[i]-xe[i-1])/h - vxe[i])
  psx_lb.append(0)
  psx_ub.append(0)

psy = []
psy_lb = []
psy_ub = []
for i in range(1,t+1):
  psy.append((ye[i]-ye[i-1])/h - vye[i])
  psy_lb.append(0)
  psy_ub.append(0)

psz = []
psz_lb = []
psz_ub = []
for i in range(1,t+1):
  psz.append((ze[i]-ze[i-1])/h - vze[i])
  psz_lb.append(0)
  psz_ub.append(0)

sc = []
sc_lb = []
sc_ub = []
for i in range(1,t+1):
  sc.append((ye[i] - sin(xe[i]))**2 + (ze[i] - cos(xe[i]))**2 - tol**2)
  sc_lb.append(-inf)
  sc_ub.append(0)
  
g = vertcat([acx,acy,acz,psx,psy,psz,sc])
g_lb = NP.concatenate([acx_lb,acy_lb,acz_lb,psx_lb,psy_lb,psz_lb,sc_lb])
g_ub = NP.concatenate([acx_ub,acy_ub,acz_ub,psx_ub,psy_ub,psz_ub,sc_ub])

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




