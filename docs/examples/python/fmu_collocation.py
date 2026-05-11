#
#     MIT No Attribution
#
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
#
#     Permission is hereby granted, free of charge, to any person obtaining a copy of this
#     software and associated documentation files (the "Software"), to deal in the Software
#     without restriction, including without limitation the rights to use, copy, modify,
#     merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
#     permit persons to whom the Software is furnished to do so.
#
#     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#     INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#     PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
#
# -*- coding: utf-8 -*-
"""
@authors: Mario Zanon and Sebastien Gross 2012, Joel Andersson 2026
"""
from casadi import *
from numpy import inf
import numpy as np
import matplotlib.pyplot as plt

# Dynamic model: 
m = 1.
M = 1.
g = 9.81

# Control
u  = SX.sym("u")

# Differential states
x = SX.sym("x")
y = SX.sym("y")
w = SX.sym("w")
dx = SX.sym("dx")
dy = SX.sym("dy")
dw = SX.sym("dw")
X = vertcat(x,y,w,dx,dy,dw)

# Algebraic variables
Z = SX.sym("Z")

# Ordinary differential equations
ddx = -(x-w)*Z/m
ddy = g - y*Z/m
ddw = ((x-w)*Z - u)/M
Xdot = vertcat(dx, dy, dw, ddx, ddy, ddw)

# Algebraic equation
Alg = (x-w)*(ddx - ddw) + y*ddy + dy*dy + (dx-dw)*(dx-dw)

# DAE rhs
ffcn = Function('ffcn', [X,Z,u],[Xdot,Alg])

# Objective function
xref = 0.1 # chariot reference
MayerTerm = Function('mayer', [X,Z,u],[(x-xref)*(x-xref) + (w-xref)*(w-xref) + dx*dx + dy*dy])
LagrangeTerm = Function('lagrange', [X,Z,u],[(x-xref)*(x-xref) + (w-xref)*(w-xref)])


# -----------------------------------------------------------------------------
# Collocation setup
# -----------------------------------------------------------------------------
l = 1. #- -> crane, + -> pendulum
tf = 5.0
nk = 50

# Degree of interpolating polynomial
deg = 4
# Radau collocation points
cp = "radau"
# Size of the finite elements
h = tf/nk

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

# Control bounds
u_min = np.array([-2])
u_max = np.array([ 2])
u_init = np.array((nk*(deg+1))*[[0.0]]) # needs to be specified for every time interval (even though it stays constant)

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
xD_init = np.array((nk*(deg+1))*[[ 0.0,  l,  0.0,  0.0,  0.0,  0.0]]) # needs to be specified for every time interval

# Algebraic state bounds and initial guess
xA_min =  np.array([-inf])
xA_max =  np.array([ inf])
xAi_min = np.array([-inf])
xAi_max = np.array([ inf])
xAf_min = np.array([-inf])
xAf_max = np.array([ inf])
xA_init = np.array((nk*(deg+1))*[[sign(l)*9.81]])

# -----------------------------------------------------------------------------
# NLP setup
# -----------------------------------------------------------------------------
# Dimensions of the problem
nx = X.nnz() + Z.nnz()  # total number of states        #MODIF
ndiff = X.nnz()           # number of differential states #MODIF
nalg = Z.nnz()            # number of algebraic states
nu = u.nnz()               # number of controls

# Total number of variables
NXD = nk*(deg+1)*ndiff # Collocated differential states
NXA = nk*deg*nalg      # Collocated algebraic states
NU = nk*nu                  # Parametrized controls
NXF = ndiff                 # Final state (only the differential states)
NV = NXD+NXA+NU+NXF

# NLP variable vector
V = MX.sym("V",NV)

# All variables with bounds and initial guess
vars_lb = np.zeros(NV)
vars_ub = np.zeros(NV)
vars_init = np.zeros(NV)
offset = 0

# Get collocated states and parametrized control
XD = np.resize(np.array([],dtype=MX),(nk+1,deg+1)) # NB: same name as above
XA = np.resize(np.array([],dtype=MX),(nk,deg)) # NB: same name as above
U = np.resize(np.array([],dtype=MX),nk)
for k in range(nk):
    # Collocated states
    for j in range(deg+1):

        # Get the expression for the state vector
        XD[k][j] = V[offset:offset+ndiff]
        if j !=0:
            XA[k][j-1] = V[offset+ndiff:offset+ndiff+nalg]
        # Add the initial condition
        index = (deg+1)*k + j
        if k==0 and j==0:
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
XD[nk][0] = V[offset:offset+ndiff]
vars_lb[offset:offset+ndiff] = xDf_min
vars_ub[offset:offset+ndiff] = xDf_max
vars_init[offset:offset+ndiff] = xD_init[-1,:]
offset += ndiff
assert offset==NV

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
            xp_jk += C[j2][j]*XD[k][j2]       # get the time derivative of the differential states (eq 10.19b)

        # Add collocation equations to the NLP
        [Xdotk, Algk] = ffcn(XD[k][j], XA[k][j-1], U[k])
        # impose system dynamics (for the differential states (eq 10.19b))
        g += [Xdotk - xp_jk/h]
        lbg.append(np.zeros(ndiff)) # equality constraints
        ubg.append(np.zeros(ndiff)) # equality constraints
        # impose system dynamics (for the algebraic states (eq 10.19b))
        g += [Algk]
        lbg.append(np.zeros(nalg)) # equality constraints
        ubg.append(np.zeros(nalg)) # equality constraints

    # Get an expression for the state at the end of the finite element
    xf_k = 0
    for j in range(deg+1):
        xf_k += D[j]*XD[k][j]

    # Add continuity equation to NLP
    g += [XD[k+1][0] - xf_k]
    lbg.append(np.zeros(ndiff))
    ubg.append(np.zeros(ndiff))

# Objective function of the NLP
#Implement Mayer term
Obj = 0
obj = MayerTerm(XD[k][j], XA[k][j-1], U[k])
Obj += obj

# Implement Lagrange term
lDotAtTauRoot = C.T
lAtOne = D

ldInv = np.linalg.inv(lDotAtTauRoot[1:,1:])
ld0 = lDotAtTauRoot[1:,0]
lagrangeTerm = 0
for k in range(nk):
    dQs = h*veccat(*[LagrangeTerm(XD[k][j], XA[k][j-1], U[k]) \
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
#opts["ipopt.linear_solver"] = 'ma27'

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
xD_opt = np.resize(np.array([],dtype=MX),(ndiff,(deg+1)*(nk)+1))
xA_opt = np.resize(np.array([],dtype=MX),(nalg,(deg)*(nk)))
u_opt = np.resize(np.array([],dtype=MX),(nu,(deg+1)*(nk)+1))
offset = 0
offset2 = 0
offset3 = 0
offset4 = 0

for k in range(nk):
    for j in range(deg+1):
        xD_opt[:,offset2] = v_opt[offset:offset+ndiff][:,0]
        offset2 += 1
        offset += ndiff
        if j!=0:
            xA_opt[:,offset4] = v_opt[offset:offset+nalg][:,0]
            offset4 += 1
            offset += nalg
    utemp = v_opt[offset:offset+nu][:,0]
    for j in range(deg+1):
        u_opt[:,offset3] = utemp
        offset3 += 1
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

xA_plt = np.resize(np.array([],dtype=MX),(nalg,(deg+1)*(nk)+1))
offset4=0
offset5=0
for k in range(nk):
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
for k in range(nk):
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
