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
"""
Created on Mon Jan 16 10:57:20 2012

@author: mzanon
"""

from casadi import *
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# Collocation setup
# -----------------------------------------------------------------------------
nicp = 1        # Number of (intermediate) collocation points per control interval
nk = 10         # Control discretization
tf = 2.5        # End time

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
deg = 4
# Radau collocation points
cp = RADAU
# Size of the finite elements
h = tf/nk/nicp

# Coefficients of the collocation equation
C = np.zeros((deg+1,deg+1))b
# Coefficients of the continuity equation
D = np.zeros(deg+1)

# Collocation point
tau = ssym("tau")
  
# All collocation time points
tau_root = collocation_points[cp][deg]
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
    lfcn = SXFunction([tau],[L])
    lfcn.init()
    # Evaluate the polynomial at the final time to get the coefficients of the continuity equation
    lfcn.setInput(1.0)
    lfcn.evaluate()
    D[j] = lfcn.output()
    # Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
    for j2 in range(deg+1):
        lfcn.setInput(tau_root[j2])
        lfcn.setFwdSeed(1.0)
        lfcn.evaluate(1,0)
        C[j][j2] = lfcn.fwdSens()



# -----------------------------------------------------------------------------
# Model setup
# -----------------------------------------------------------------------------
# Declare variables (use scalar graph)
t  = ssym("t")          # time
u  = ssym("u")          # control
xd  = ssym("xd",3)      # differential state
xa  = ssym("xa",0,1)    # algebraic state
xddot  = ssym("xdot",3) # differential state time derivative
p = ssym("xa",0,1)      # parameters

# ODE right hand side function
rhs = [(1 - xd[1]*xd[1])*xd[0] - xd[1] + u, \
       xd[0], \
       xd[0]*xd[0] + xd[1]*xd[1] + u*u]
# System dynamics (implicit formulation)
ffcn = SXFunction([t,xddot,xd,xa,u,p],[xddot - rhs])

# Objective function (meyer term)
mfcn = SXFunction([t,xd,xa,u,p],[xd[2]])

# Control bounds
u_min = np.array([-0.75])
u_max = np.array([ 1.0])
u_init = np.array((nk*nicp*(deg+1))*[[0.0]]) # needs to be specified for every time interval (even though it stays constant)

# Differential state bounds and initial guess
xD_min =  np.array([-inf, -inf, -inf])
xD_max =  np.array([ inf,  inf,  inf])
xDi_min = np.array([ 0.0,  1.0,  0.0])
xDi_max = np.array([ 0.0,  1.0,  0.0])
xDf_min = np.array([ 0.0,  0.0, -inf])
xDf_max = np.array([ 0.0,  0.0,  inf])
xD_init = np.array((nk*nicp*(deg+1))*[[ 0.0,  0.0,  0.0]]) # needs to be specified for every time interval

# Algebraic state bounds and initial guess
xA_min =  np.array([])
xA_max =  np.array([])
xAi_min = np.array([])
xAi_max = np.array([])
xAf_min = np.array([])
xAf_max = np.array([])
xA_init = np.array((nk*nicp*(deg+1))*[[]])

# Parameter bounds and initial guess
p_min = np.array([])
p_max = np.array([])
p_init = np.array([])

# Initialize functions
ffcn.init()
mfcn.init()


# -----------------------------------------------------------------------------
# Constraints setup
# -----------------------------------------------------------------------------
# Initial constraint
ic_min = np.array([])
ic_max = np.array([])
ic = SXMatrix()
#ic.append();       ic_min = append(ic_min, 0.);         ic_max = append(ic_max, 0.)
icfcn = SXFunction([t,xd,xa,u,p],[ic])
# Path constraint
pc_min = np.array([])
pc_max = np.array([])
pc = SXMatrix()
#pc.append();       pc_min = append(pc_min, 0.);         pc_max = append(pc_max, 0.)
pcfcn = SXFunction([t,xd,xa,u,p],[pc])
# Final constraint
fc_min = np.array([])
fc_max = np.array([])
fc = SXMatrix()
#fc.append();       fc_min = append(fc_min, 0.);         fc_max = append(fc_max, 0.)
fcfcn = SXFunction([t,xd,xa,u,p],[fc])

# Initialize the functions
icfcn.init()
pcfcn.init()
fcfcn.init()

# -----------------------------------------------------------------------------
# NLP setup
# -----------------------------------------------------------------------------
# Dimensions of the problem
nx = xd.size()          # total number of states
ndiff = nx              # number of differential states
nalg = xa.size()        # number of algebraic states
nu = u.size()           # number of controls
NP  = p.size()          # number of parameters

# Total number of variables
NXD = nicp*nk*(deg+1)*ndiff # Collocated differential states
NXA = nicp*nk*deg*nalg      # Collocated algebraic states
NU = nk*nu                  # Parametrized controls
NXF = ndiff                 # Final state (only the differential states)
NV = NXD+NXA+NU+NXF+NP

# NLP variable vector
V = msym("V",NV)
  
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
[ick] = icfcn.call([0., XD[0][0][0], XA[0][0][0], U[0], P])
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
            [fk] = ffcn.call([0., xp_jk/h, XD[k][i][j], XA[k][i][j-1], U[k], P])
            g += [fk[:ndiff]]                     # impose system dynamics (for the differential states (eq 10.19b))
            lbg.append(np.zeros(ndiff)) # equality constraints
            ubg.append(np.zeros(ndiff)) # equality constraints
            g += [fk[ndiff:]]                               # impose system dynamics (for the algebraic states (eq 10.19b))
            lbg.append(np.zeros(nalg)) # equality constraints
            ubg.append(np.zeros(nalg)) # equality constraints
            
            #  Evaluate the path constraint function
            [pck] = pcfcn.call([0., XD[k][i][j], XA[k][i][j-1], U[k], P])
            
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


# Final constraints (Const, dConst, ConstQ)
[fck] = fcfcn.call([0., XD[k][i][j], XA[k][i][j-1], U[k], P])
g += [fck]
lbg.append(fc_min)
ubg.append(fc_max)


# Nonlinear constraint function
gfcn = MXFunction([V],[vertcat(g)])

# Objective function of the NLP
# Regularization
Obj = 0
#for k in range(nk):
#    for i in range(nicp):
#        # For all collocation points
#        for j in range(1,deg+1):
#            [obj] = mfcn2.call([XD[k][i][j],U[k], P])
#            Obj += obj
#[obj] = mfcn2.call([XD[nk][0][0],zeros(9), P])
#Obj += obj
# Energy
[obj] = mfcn.call([0., XD[k][i][j], XA[k][i][j-1], U[k], P])
Obj += obj

ofcn = MXFunction([V], [Obj])


## ----
## SOLVE THE NLP
## ----
  
#assert(1==0)
# Allocate an NLP solver
solver = IpoptSolver(ofcn,gfcn)

# Set options
solver.setOption("expand_f",True)
solver.setOption("expand_g",True)
solver.setOption("generate_hessian",True)
solver.setOption("max_iter",1000)
solver.setOption("tol",1e-4)
solver.setOption("linear_solver","ma27")

# initialize the solver
solver.init()
  
# Initial condition
solver.setInput(vars_init,NLP_X_INIT)

# Bounds on x
solver.setInput(vars_lb,NLP_LBX)
solver.setInput(vars_ub,NLP_UBX)

# Bounds on g
solver.setInput(np.concatenate(lbg),NLP_LBG)
solver.setInput(np.concatenate(ubg),NLP_UBG)

# Solve the problem
solver.solve()

# Print the optimal cost
print "optimal cost: ", float(solver.output(NLP_COST))

# Retrieve the solution
v_opt = np.array(solver.output(NLP_X_OPT))
    

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
    lafcn = SXFunction([tau],[La])
    lafcn.init()
    lafcn.setInput(tau_root[0])
    lafcn.evaluate()
    Da[j-1] = lafcn.output()

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
plt.plot(tgrid,xD_opt[0,:],'--')
plt.plot(tgrid,xD_opt[1,:],'-')
plt.plot(tgrid,u_opt[0,:],'-.')
plt.title("Van der Pol optimization")
plt.xlabel('time')
plt.legend(['x0 trajectory','x1 trajectory','u trajectory'])
plt.grid()
plt.show()

