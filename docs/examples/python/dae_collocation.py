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
"""
@author: Mario Zanon and Sebastien Gross, K.U. Leuven 2012
"""

from casadi import *
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# Collocation setup
# -----------------------------------------------------------------------------
nicp = 1        # Number of (intermediate) collocation points per control interval

xref = 0.1 # chariot reference

l = 1. #- -> crane, + -> pendulum
m = 1.
M = 1.
g = 9.81
tf = 5.0
nk = 50
ndstate = 6
nastate = 1
ninput = 1

# Degree of interpolating polynomial
deg = 4
# Radau collocation points
cp = "radau"
# Size of the finite elements
h = tf/nk/nicp

# Coefficients of the collocation equation
C = np.zeros((deg+1,deg+1))
# Coefficients of the continuity equation
D = np.zeros(deg+1)

# Collocation point
tau = SX.sym("tau")

# All collocation time points
tau_root = [0] + collocation_points(deg, cp)

T = np.zeros((nk,deg+1))
for i in range(nk):
    for j in range(deg+1):
        T[i][j] = h*(i + tau_root[j])

# For all collocation points: eq 10.4 or 10.17 in Biegler's book
# Construct Lagrange polynomials to get the polynomial basis at the collocation point
for j in range(deg+1):
    L = 1
    for j2 in range(deg+1):
        if j2 != j:
            L *= (tau-tau_root[j2])/(tau_root[j]-tau_root[j2])

    # Evaluate the polynomial at the final time to get the coefficients of the continuity equation
    lfcn = Function('lfcn', [tau],[L])
    D[j] = lfcn(1.0)

    # Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
    tfcn = Function('tfcn', [tau],[tangent(L,tau)])
    for j2 in range(deg+1):
        C[j][j2] = tfcn(tau_root[j2])

# -----------------------------------------------------------------------------
# Model setup
# -----------------------------------------------------------------------------
# Declare variables (use scalar graph)
t  = SX.sym("t")          # time
u  = SX.sym("u")          # control
xd  = SX.sym("xd",ndstate)      # differential state
xa  = SX.sym("xa",nastate)    # algebraic state
xddot  = SX.sym("xdot",ndstate) # differential state time derivative
p = SX.sym("p",0,1)      # parameters

x = SX.sym("x")
y = SX.sym("y")
w = SX.sym("w")

dx = SX.sym("dx")
dy = SX.sym("dy")
dw = SX.sym("dw")


res = vertcat(xddot[0] - dx,\
       xddot[1] - dy,\
       xddot[2] - dw,\
       m*xddot[3] + (x-w)*xa, \
       m*xddot[4] +     y*xa - g*m,\
       M*xddot[5] + (w-x)*xa +   u,\
       (x-w)*(xddot[3] - xddot[5]) + y*xddot[4] + dy*dy + (dx-dw)*(dx-dw))


xd[0] = x
xd[1] = y
xd[2] = w
xd[3] = dx
xd[4] = dy
xd[5] = dw


# System dynamics (implicit formulation)
ffcn = Function('ffcn', [t,xddot,xd,xa,u,p],[res])

# Objective function
MayerTerm = Function('mayer', [t,xd,xa,u,p],[(x-xref)*(x-xref) + (w-xref)*(w-xref) + dx*dx + dy*dy])
LagrangeTerm = Function('lagrange', [t,xd,xa,u,p],[(x-xref)*(x-xref) + (w-xref)*(w-xref)])

# Control bounds
u_min = np.array([-2])
u_max = np.array([ 2])
u_init = np.array((nk*nicp*(deg+1))*[[0.0]]) # needs to be specified for every time interval (even though it stays constant)

# Differential state bounds
#Path bounds
xD_min =  np.array([-inf, -inf, -inf, -inf, -inf, -inf])
xD_max =  np.array([ inf,  inf,  inf,  inf,  inf,  inf])
#Initial bounds
xDi_min = np.array([ 0.0,  l,  0.0,  0.0,  0.0,  0.0])
xDi_max = np.array([ 0.0,  l,  0.0,  0.0,  0.0,  0.0])
#Final bounds
xDf_min = np.array([-inf, -inf, -inf, -inf, -inf, -inf])
xDf_max = np.array([ inf,  inf,  inf,  inf,  inf,  inf])

#Initial guess for differential states
xD_init = np.array((nk*nicp*(deg+1))*[[ 0.0,  l,  0.0,  0.0,  0.0,  0.0]]) # needs to be specified for every time interval

# Algebraic state bounds and initial guess
xA_min =  np.array([-inf])
xA_max =  np.array([ inf])
xAi_min = np.array([-inf])
xAi_max = np.array([ inf])
xAf_min = np.array([-inf])
xAf_max = np.array([ inf])
xA_init = np.array((nk*nicp*(deg+1))*[[sign(l)*9.81]])

# Parameter bounds and initial guess
p_min = np.array([])
p_max = np.array([])
p_init = np.array([])




# -----------------------------------------------------------------------------
# Constraints setup
# -----------------------------------------------------------------------------
# Initial constraint
ic_min = np.array([])
ic_max = np.array([])
ic = SX()
#ic.append();       ic_min = append(ic_min, 0.);         ic_max = append(ic_max, 0.)
icfcn = Function('icfcn', [t,xd,xa,u,p],[ic])
# Path constraint
pc_min = np.array([])
pc_max = np.array([])
pc = SX()
#pc.append();       pc_min = append(pc_min, 0.);         pc_max = append(pc_max, 0.)
pcfcn = Function('pcfcn', [t,xd,xa,u,p],[pc])
# Final constraint
fc_min = np.array([])
fc_max = np.array([])
fc = SX()
#fc.append();       fc_min = append(fc_min, 0.);         fc_max = append(fc_max, 0.)
fcfcn = Function('fcfcn', [t,xd,xa,u,p],[fc])

# -----------------------------------------------------------------------------
# NLP setup
# -----------------------------------------------------------------------------
# Dimensions of the problem
nx = xd.nnz() + xa.nnz()  # total number of states        #MODIF
ndiff = xd.nnz()           # number of differential states #MODIF
nalg = xa.nnz()            # number of algebraic states
nu = u.nnz()               # number of controls
NP  = p.nnz()              # number of parameters

# Total number of variables
NXD = nicp*nk*(deg+1)*ndiff # Collocated differential states
NXA = nicp*nk*deg*nalg      # Collocated algebraic states
NU = nk*nu                  # Parametrized controls
NXF = ndiff                 # Final state (only the differential states)
NV = NXD+NXA+NU+NXF+NP

# NLP variable vector
V = MX.sym("V",NV)

# All variables with bounds and initial guess
vars_lb = np.zeros(NV)
vars_ub = np.zeros(NV)
vars_init = np.zeros(NV)
offset = 0

# Get the parameters
P = V[offset:offset+NP]
vars_init[offset:offset+NP] = p_init
vars_lb[offset:offset+NP] = p_min
vars_ub[offset:offset+NP] = p_max
offset += NP

# Get collocated states and parametrized control
XD = np.resize(np.array([],dtype=MX),(nk+1,nicp,deg+1)) # NB: same name as above
XA = np.resize(np.array([],dtype=MX),(nk,nicp,deg)) # NB: same name as above
U = np.resize(np.array([],dtype=MX),nk)
for k in range(nk):
    # Collocated states
    for i in range(nicp):
        #
        for j in range(deg+1):

            # Get the expression for the state vector
            XD[k][i][j] = V[offset:offset+ndiff]
            if j !=0:
                XA[k][i][j-1] = V[offset+ndiff:offset+ndiff+nalg]
            # Add the initial condition
            index = (deg+1)*(nicp*k+i) + j
            if k==0 and j==0 and i==0:
                vars_init[offset:offset+ndiff] = xD_init[index,:]

                vars_lb[offset:offset+ndiff] = xDi_min
                vars_ub[offset:offset+ndiff] = xDi_max
                offset += ndiff
            else:
                if j!=0:
                    vars_init[offset:offset+nx] = np.append(xD_init[index,:],xA_init[index,:])

                    vars_lb[offset:offset+nx] = np.append(xD_min,xA_min)
                    vars_ub[offset:offset+nx] = np.append(xD_max,xA_max)
                    offset += nx
                else:
                    vars_init[offset:offset+ndiff] = xD_init[index,:]

                    vars_lb[offset:offset+ndiff] = xD_min
                    vars_ub[offset:offset+ndiff] = xD_max
                    offset += ndiff

    # Parametrized controls
    U[k] = V[offset:offset+nu]
    vars_lb[offset:offset+nu] = u_min
    vars_ub[offset:offset+nu] = u_max
    vars_init[offset:offset+nu] = u_init[index,:]
    offset += nu

# State at end time
XD[nk][0][0] = V[offset:offset+ndiff]
vars_lb[offset:offset+ndiff] = xDf_min
vars_ub[offset:offset+ndiff] = xDf_max
vars_init[offset:offset+ndiff] = xD_init[-1,:]
offset += ndiff
assert(offset==NV)

# Constraint function for the NLP
g = []
lbg = []
ubg = []

# Initial constraints
ick = icfcn(0., XD[0][0][0], XA[0][0][0], U[0], P)
g += [ick]
lbg.append(ic_min)
ubg.append(ic_max)

# For all finite elements
for k in range(nk):
    for i in range(nicp):
        # For all collocation points
        for j in range(1,deg+1):
            # Get an expression for the state derivative at the collocation point
            xp_jk = 0
            for j2 in range (deg+1):
                xp_jk += C[j2][j]*XD[k][i][j2]       # get the time derivative of the differential states (eq 10.19b)

            # Add collocation equations to the NLP
            fk = ffcn(0., xp_jk/h, XD[k][i][j], XA[k][i][j-1], U[k], P)
            g += [fk[:ndiff]]                     # impose system dynamics (for the differential states (eq 10.19b))
            lbg.append(np.zeros(ndiff)) # equality constraints
            ubg.append(np.zeros(ndiff)) # equality constraints
            g += [fk[ndiff:]]                               # impose system dynamics (for the algebraic states (eq 10.19b))
            lbg.append(np.zeros(nalg)) # equality constraints
            ubg.append(np.zeros(nalg)) # equality constraints

            #  Evaluate the path constraint function
            pck = pcfcn(0., XD[k][i][j], XA[k][i][j-1], U[k], P)

            g += [pck]
            lbg.append(pc_min)
            ubg.append(pc_max)

        # Get an expression for the state at the end of the finite element
        xf_k = 0
        for j in range(deg+1):
            xf_k += D[j]*XD[k][i][j]

        # Add continuity equation to NLP
        if i==nicp-1:
#            print "a ", k, i
            g += [XD[k+1][0][0] - xf_k]
        else:
#            print "b ", k, i
            g += [XD[k][i+1][0] - xf_k]

        lbg.append(np.zeros(ndiff))
        ubg.append(np.zeros(ndiff))

# Periodicity constraints
#   none

# Final constraints (Const, dConst, ConstQ)
fck = fcfcn(0., XD[k][i][j], XA[k][i][j-1], U[k], P)
g += [fck]
lbg.append(fc_min)
ubg.append(fc_max)

# Objective function of the NLP
#Implement Mayer term
Obj = 0
obj = MayerTerm(0., XD[k][i][j], XA[k][i][j-1], U[k], P)
Obj += obj

# Implement Lagrange term
lDotAtTauRoot = C.T
lAtOne = D

ldInv = np.linalg.inv(lDotAtTauRoot[1:,1:])
ld0 = lDotAtTauRoot[1:,0]
lagrangeTerm = 0
for k in range(nk):
    for i in range(nicp):
        dQs = h*veccat(*[LagrangeTerm(0., XD[k][i][j], XA[k][i][j-1], U[k], P) \
                        for j in range(1,deg+1)])
        Qs = mtimes( ldInv, dQs)
        m = mtimes( Qs.T, lAtOne[1:])
        lagrangeTerm += m

Obj += lagrangeTerm

# NLP
nlp = {'x':V, 'f':Obj, 'g':vertcat(*g)}

## ----
## SOLVE THE NLP
## ----

# NLP solver options
opts = {}
opts["expand"] = True
opts["ipopt.max_iter"] = 1000
opts["ipopt.tol"] = 1e-4
opts["ipopt.linear_solver"] = 'ma27'

# Allocate an NLP solver
solver = nlpsol("solver", "ipopt", nlp, opts)
arg = {}

# Initial condition
arg["x0"] = vars_init

# Bounds on x
arg["lbx"] = vars_lb
arg["ubx"] = vars_ub

# Bounds on g
arg["lbg"] = np.concatenate(lbg)
arg["ubg"] = np.concatenate(ubg)

# Solve the problem
res = solver(**arg)

# Print the optimal cost
print("optimal cost: ", float(res["f"]))

# Retrieve the solution
v_opt = np.array(res["x"])


## ----
## RETRIEVE THE SOLUTION
## ----
xD_opt = np.resize(np.array([],dtype=MX),(ndiff,(deg+1)*nicp*(nk)+1))
xA_opt = np.resize(np.array([],dtype=MX),(nalg,(deg)*nicp*(nk)))
u_opt = np.resize(np.array([],dtype=MX),(nu,(deg+1)*nicp*(nk)+1))
offset = 0
offset2 = 0
offset3 = 0
offset4 = 0

for k in range(nk):
    for i in range(nicp):
        for j in range(deg+1):
            xD_opt[:,offset2] = v_opt[offset:offset+ndiff][:,0]
            offset2 += 1
            offset += ndiff
            if j!=0:
                xA_opt[:,offset4] = v_opt[offset:offset+nalg][:,0]
                offset4 += 1
                offset += nalg
    utemp = v_opt[offset:offset+nu][:,0]
    for i in range(nicp):
        for j in range(deg+1):
            u_opt[:,offset3] = utemp
            offset3 += 1
    #    u_opt += v_opt[offset:offset+nu]
    offset += nu

xD_opt[:,-1] = v_opt[offset:offset+ndiff][:,0]



# The algebraic states are not defined at the first collocation point of the finite elements:
# with the polynomials we compute them at that point
Da = np.zeros(deg)
for j in range(1,deg+1):
    # Lagrange polynomials for the algebraic states: exclude the first point
    La = 1
    for j2 in range(1,deg+1):
        if j2 != j:
            La *= (tau-tau_root[j2])/(tau_root[j]-tau_root[j2])
    lafcn = Function('lafcn', [tau], [La])
    Da[j-1] = lafcn(tau_root[0])

xA_plt = np.resize(np.array([],dtype=MX),(nalg,(deg+1)*nicp*(nk)+1))
offset4=0
offset5=0
for k in range(nk):
    for i in range(nicp):
        for j in range(deg+1):
            if j!=0:
                xA_plt[:,offset5] = xA_opt[:,offset4]
                offset4 += 1
                offset5 += 1
            else:
                xa0 = 0
                for j in range(deg):
                    xa0 += Da[j]*xA_opt[:,offset4+j]
                xA_plt[:,offset5] = xa0
                #xA_plt[:,offset5] = xA_opt[:,offset4]
                offset5 += 1

xA_plt[:,-1] = xA_plt[:,-2]




tg = np.array(tau_root)*h
for k in range(nk*nicp):
    if k == 0:
        tgrid = tg
    else:
        tgrid = np.append(tgrid,tgrid[-1]+tg)
tgrid = np.append(tgrid,tgrid[-1])
# Plot the results
plt.figure(1)
plt.clf()
plt.subplot(2,2,1)
plt.plot(tgrid,xD_opt[0,:],'--')
plt.title("x")
plt.grid
plt.subplot(2,2,2)
plt.plot(tgrid,xD_opt[1,:],'-')
plt.title("y")
plt.grid
plt.subplot(2,2,3)
plt.plot(tgrid,xD_opt[2,:],'-.')
plt.title("w")
plt.grid

plt.figure(2)
plt.clf()
plt.plot(tgrid,u_opt[0,:],'-.')
plt.title("Crane, inputs")
plt.xlabel('time')


plt.figure(3)
plt.clf()
plt.plot(tgrid,xA_plt[0,:],'-.')
plt.title("Crane, lambda")
plt.xlabel('time')
plt.grid()
plt.show()
