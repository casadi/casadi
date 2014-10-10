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
from numpy import *
import numpy as NP
import matplotlib.pyplot as plt
import copy

# Excercise 2, chapter 10 from Larry Biegler's book 
# Joel Andersson, K.U. Leuven 2010

nk = 50    # Control discretization

# Time 
t = SX("t")
  
# Differential states
Ls = SX("Ls")    # mean crystal size
Nc = SX("Nc")    # number of nuclei per liter of solvent
L = SX("L")      # total length of crystals per liter of solvent
Ac = SX("Ac")    # total surface area of the crystals per liter of solvent
Vc = SX("Vc")    # total volume of the crysals per liter of solvent
Mc = SX("Mc")    # total mass of the crystals
Cc = SX("Cc")    # solute concentration
Tc = SX("Tc")    # cystillizer temperature

# Initial values
Ls_init = 0.0005
Nc_init = 0.
L_init = 0.
Ac_init = 0.
Vc_init = 0.
Mc_init = 2.0
Cc_init = 5.4
Tc_init = 75.

# Controls
Tj = SX("Tj") # jacket temperature
Tj_lb = 10.;  Tj_ub = 60.;  Tj_init = 50.

# Constants
Vs = 300. # volume of the solvent
W = 2025. # the total mass in the crysallizer
a = [-66.4309, 2.8604, -0.022579, 6.7117e-5]
b = [16.08852, -2.708263, 0.0670694, -3.5685e-4]
Kg = 0.00418
Bn = 385.
Cp = 0.4
Kc = 35.
Ke = 377.
eta1 = 1.1
eta2 = 5.72
Ls0 = 5e-4 # initial crystal size
L0 = 5e-5 # nucleate crystal size
Ws0 = 2. # weight of the seed crystals
rho = 1.58 # specific gravity of crystals
alpha = 0.2 # shape factor for area of crystals
beta = 1.2 # shape factor for volume of crystals

# Time horizon
tf = 25.

# Dependent variables
C_bar = 100.*Cc/(1.35+Cc)
Tequ = a[0] + a[1]*C_bar + a[2]*C_bar*C_bar + a[3]*C_bar*C_bar*C_bar # equilibrium temperature
Ta = b[0] + b[1]*C_bar + b[2]*C_bar*C_bar + b[3]*C_bar*C_bar*C_bar # lower bound of Tj

# degree of supercooling:
DeltaT = fmax(0,Tequ-Tc)   # Original formulation
# DeltaT = fmax(1e-8,Tequ-Tc) # "epsilon" to avoid divide by zero
# DeltaT = log(exp(1e-8)+exp(Tequ-Tc))   # Log sum exp
#epsilon = 1e-3
#DD = Tequ-Tc
#DeltaT = (sqrt( DD*DD + epsilon*epsilon )  + DD ) / 2

# Differential equations
Ls_dot = Kg * sqrt(Ls) * DeltaT**eta1
Nc_dot = Bn * DeltaT**eta2
L_dot = Nc * Ls_dot + L0 * Nc_dot
Ac_dot = 2 * alpha * Nc * Ls_dot + L0**2 * Nc_dot
Vc_dot = 3 * beta * Ac * Ls_dot + L0**3 * Nc_dot
Mc_dot = 3 * Ws0/Ls0**3 * Ls**2 * Ls_dot + rho*Vs*Vc_dot
Cc_dot = -1. / Vs * Mc_dot
Tc_dot = (Kc*Mc_dot - Ke*(Tc-Tj))/(W*Cp)

# State vector
x = [Ls, Nc, L, Ac, Vc, Mc, Cc, Tc]

# State bounds and initial guess
x_min = repeat(-inf,len(x))
x_max = repeat( inf,len(x))
xi_min = array([Ls_init,Nc_init,L_init,Ac_init,Vc_init,Mc_init,Cc_init,Tc_init])
xi_max = copy.deepcopy(xi_min)
xf_min = copy.deepcopy(x_min)
xf_max = copy.deepcopy(x_max)
x_init = copy.deepcopy(xi_min)

# Control vector
u = [Tj]
u_init = NP.array([Tj_init])
#u_min = NP.array(49.9)
#u_max = NP.array(50.1)
u_min = [Tj_lb]
u_max = [Tj_ub]

# ODE
xdot = [Ls_dot, Nc_dot, L_dot, Ac_dot, Vc_dot, Mc_dot, Cc_dot, Tc_dot]

# Right hand side of the ODE
ffcn = SXFunction([[t],x,u],[xdot])
ffcn.init()
  
# Objective function (meyer term)
mfcn = SXFunction([[t],x,u],[[-Ls]])
mfcn.init()

# Nonlinear constraint function
cfcn = SXFunction([[t],x,u],[[Tj-Ta]])
cfcn.init()

# Bounds on nonlinear constraints
nc = 1
c_lb = NP.repeat(0,  nc)
c_ub = NP.repeat(inf,nc)

# Dimensions
nx = len(x)
nu = len(u)

print "modelling done"


# Legendre collocation points
legendre_points1 = [0,0.500000]
legendre_points2 = [0,0.211325,0.788675]
legendre_points3 = [0,0.112702,0.500000,0.887298]
legendre_points4 = [0,0.069432,0.330009,0.669991,0.930568]
legendre_points5 = [0,0.046910,0.230765,0.500000,0.769235,0.953090]
legendre_points = [0,legendre_points1,legendre_points2,legendre_points3,legendre_points4,legendre_points5]

# Radau collocation points
radau_points1 = [0,1.000000]
radau_points2 = [0,0.333333,1.000000]
radau_points3 = [0,0.155051,0.644949,1.000000]
radau_points4 = [0,0.088588,0.409467,0.787659,1.000000]
radau_points5 = [0,0.057104,0.276843,0.583590,0.860240,1.000000]
radau_points = [0,radau_points1,radau_points2,radau_points3,radau_points4,radau_points5]

# Type of collocation points
LEGENDRE = 0
RADAU = 1
collocation_points = [legendre_points,radau_points]

# Degree of interpolating polynomial
deg = 3
  
# Radau collocation points
cp = RADAU

# Size of the finite elements
h = tf/nk

# Coefficients of the collocation equation
C = zeros((deg+1,deg+1))

# Coefficients of the continuity equation
D = zeros(deg+1)
  
# Collocation point
tau = SX("tau")

# All collocation time points
tau_root = collocation_points[cp][deg]
T = zeros((nk,deg+1))
for i in range(nk):
  for j in range(deg+1):
    T[i][j] = h*(i + tau_root[j])

# For all collocation points
for j in range(deg+1):
  # Construct Lagrange polynomials to get the polynomial basis at the collocation point
  L = 1
  for j2 in range(deg+1):
    if j2 != j:
      L *= (tau-tau_root[j2])/(tau_root[j]-tau_root[j2])
  lfcn = SXFunction([tau],[L])
  lfcn.init()
  
  # Evaluate the polynomial at the final time to get the coefficients of the continuity equation
  lfcn.setInput(1.0)
  lfcn.evaluate()
  D[j] = lfcn.getOutput()

  # Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
  for j2 in range(deg+1):
    lfcn.setInput(tau_root[j2])
    lfcn.setFwdSeed(1.0)
    lfcn.evaluate(1,0)
    C[j][j2] = lfcn.getFwdSens()

# Total number of variables
NX = nk*(deg+1)*nx      # Collocated states
NU = nk*nu              # Parametrized controls
NXF = nx                # Final state
NV = NX+NU+NXF

# NLP variable vector
V = MX("V",NV)

# All variables with bounds and initial guess
vars_lb = zeros(NV)
vars_ub = zeros(NV)
vars_init = zeros(NV)
offset = 0

# Get collocated states and parametrized control
X = resize(array([],dtype=MX),(nk+1,deg+1))
U = resize(array([],dtype=MX),nk)
for k in range(nk):
  # Collocated states
  for j in range(deg+1):
    # Get the expression for the state vector
    X[k][j] = V[offset:offset+nx]
    
    # Add the initial condition
    vars_init[offset:offset+nx] = x_init
    
    # Add bounds
    if k==0 and j==0:
      vars_lb[offset:offset+nx] = xi_min
      vars_ub[offset:offset+nx] = xi_max
    else:
      vars_lb[offset:offset+nx] = x_min
      vars_ub[offset:offset+nx] = x_max
    offset += nx
  
  # Parametrized controls
  U[k] = V[offset:offset+nu]
  vars_lb[offset:offset+nu] = u_min
  vars_ub[offset:offset+nu] = u_max
  vars_init[offset:offset+nu] = u_init
  offset += nu
  
# State at end time
X[nk][0] = V[offset:offset+nx]
vars_lb[offset:offset+nx] = xf_min
vars_ub[offset:offset+nx] = xf_max
vars_init[offset:offset+nx] = x_init
offset += nx

# Constraint function for the NLP
g = []
lbg = []
ubg = []

# For all finite elements
for k in range(nk):
  
  # For all collocation points
  for j in range(1,deg+1):
        
    # Get an expression for the state derivative at the collocation point
    xp_jk = 0
    for j2 in range (deg+1):
      xp_jk += C[j2][j]*X[k][j2]
      
    # Add collocation equations to the NLP
    [fk] = ffcn.call([T[k][j], X[k][j], U[k]])
    g += [h*fk - xp_jk]
    lbg.append(zeros(nx)) # equality constraints
    ubg.append(zeros(nx)) # equality constraints

    # Add nonlinear constraints, if any
    if nc>0:
      g += cfcn.call([T[k][j], X[k][j], U[k]])
      lbg.append(c_lb)
      ubg.append(c_ub)

  # Get an expression for the state at the end of the finite element
  xf_k = 0
  for j in range(deg+1):
    xf_k += D[j]*X[k][j]

  # Add continuity equation to NLP
  g += [X[k+1][0] - xf_k]
  lbg.append(zeros(nx))
  ubg.append(zeros(nx))

# Nonlinear constraint function
gfcn_nlp = MXFunction([V],[vertcat(g)])

# Objective function of the NLP
[f] = mfcn.call([T[nk-1][deg],X[nk][0],U[nk-1]])
ffcn_nlp = MXFunction([V], [f])

## ----
## SOLVE THE NLP
## ----
  
# Allocate an NLP solver
solver = IpoptSolver(ffcn_nlp,gfcn_nlp)

# Set options
#solver.setOption("tol",1e-6)
solver.setOption("expand_f",True)
solver.setOption("expand_g",True)
#solver.setOption("generate_hessian",True)
solver.setOption("hessian_approximation","limited-memory")
#solver.setOption("derivative_test","first-order")
#solver.setOption("bound_relax_factor",0.0)
solver.setOption("max_iter",50)

# initialize the solver
solver.init()
  
# Initial condition
solver.setInput(vars_init,"x0")

# Bounds on x
solver.setInput(vars_lb,"lbx")
solver.setInput(vars_ub,"ubx")

# Bounds on g
solver.setInput(NP.concatenate(lbg),"lbg")
solver.setInput(NP.concatenate(ubg),"ubg")

# Solve the problem
solver.solve()

# Print the optimal cost
print "optimal cost: ", solver.getOutput("f")[0]

# Retrieve the solution
v_opt = array(solver.getOutput("x"))

# Get the solution
vars_sol = solver.getOutput("x")

## ----
## SAVE SOLUTION TO DISK
## ----

Ls_opt = v_opt[0::(deg+1)*nx+nu]
Nc_opt = v_opt[1::(deg+1)*nx+nu]
L_opt = v_opt[2::(deg+1)*nx+nu]
Ac_opt = v_opt[3::(deg+1)*nx+nu]
Vc_opt = v_opt[4::(deg+1)*nx+nu]
Mc_opt = v_opt[5::(deg+1)*nx+nu]
Cc_opt = v_opt[6::(deg+1)*nx+nu]
Tc_opt = v_opt[7::(deg+1)*nx+nu]
Tj_opt = v_opt[(deg+1)*nx::(deg+1)*nx+nu]

tgrid = linspace(0,tf,nk+1)
tgrid_u = linspace(0,tf,nk)
  
# plot to screen
plt.figure(1)
plt.clf()
plt.plot(tgrid,Ls_opt)
plt.xlabel("Time (h)")
plt.ylabel("Mean Chrystal Size (m)")

plt.figure(2)
plt.clf()
plt.plot(tgrid_u,Tj_opt)
plt.xlabel("Time (h)")
plt.ylabel("Jacket Temperature (C)")

#plt.figure(3)
#plt.clf()
#plt.plot(t_opt,Ls_opt)


# show the plots
plt.show()
