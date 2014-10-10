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
from casadi import *
from numpy import *
from os import system
import time
import sys

CN = 1.0

# First order model parameters of velocity loops for trolley and hoisting

# Gain of the closed loop transfer function for the trolley.
# Unit: [m/V]
FOM_A1 = 0.042152000000000

# Gain of the closed loop transfer function for hoisting
# Unit: [m/V]
FOM_A2 = 0.029488576312860

# Second order model parameters of velocity loops for trolley and hoisting

SOM_A1 = CN * 0.047418203070092
SOM_TAU1 = 0.012790605943772

SOM_A2 =0.034087337273386
SOM_TAU2 = 0.024695192379264

# Define the variables
#
xT = ssym("xT")         # the trolley position
vT = ssym("vT")         # the trolley velocity
xL = ssym("xL")         # the cable length
vL = ssym("vL")         # the cable velocity
theta = ssym("theta")   # the excitation angle
omega = ssym("omega")   # the angular velocity
uT = ssym("uT")         # the input to the trolley velocity controller
uL = ssym("uL")         # the input to the cable velocity controller
duT = ssym("duT")       # cart control rate
duL = ssym("duL")       # winch control rate
w = ssym("w",6)         # disturbances

# State vector
x = vertcat((xT,vT,xL,vL,theta,omega,uT,uL))
print "x = ", x

# Controls
u = vertcat((duT,duL,w))
print "u = ", u

# Set up the MPC - Optimal control problem
t0 = 0.0
tf = 1.0
nk = 10

tau1 = SOM_TAU1        # time constant
a1   = SOM_A1          # gain

#
# XL( s ) / UL( s ) = A2 / ( s * ( tau2 * s + 1 ) )
#
tau2 = SOM_TAU2        # time constant
a2   = SOM_A2          # gain

#
# Additional parameters
#
g = 9.81       # the gravitational constant
c = 0.0        # damping constant for the motion of the pendulum
m = 1318.0     # the mass of the pendulum

#
# Define the model equations
#

# the trolley acceleration
aT = FOM_A1 * duT
       
# the cable acceleration
aL = -1.0 / tau2 * vL + a2 / tau2 * uL

# ODE right hand side
dot_xT = vT
dot_vT = aT
dot_xL = vL
dot_vL = aL
dot_theta = omega
dot_omega = 1.0 / xL * (-g * sin( theta ) - aT * cos( theta ) - 2 * vL * omega - c * omega / ( m * xL ) )
dot_uT = duT
dot_uL = duL
xdot = vertcat((dot_xT,dot_vT,dot_xL,dot_vL,dot_theta,dot_omega,dot_uT,dot_uL))
print "xdot = ", xdot

# ODE right hand side function
f = SXFunction([x,u],[xdot])
f.init()

# Create an expression for the Jacobian of rhs with respect to x
df_dx = f.jac(0,0)
makeDense(df_dx) # make the expression dense
print "df_dx = ", df_dx

# Create an expression for the Jacobian of rhs with respect to u
df_du = f.jac(1,0)
makeDense(df_du) # make the expression dense
print "df_du = ", df_du

# Create a function which evaluates f, df_dx and df_du
jac_f = SXFunction([x,u],[xdot,df_dx,df_du])
#jac_f.setOption("topological_sorting","depth-first")
jac_f.setOption("live_variables",True)
jac_f.setOption("verbose",True)
#jac_f.setOption("inplace",True)
jac_f.init()

# Generate c code
jac_f.generateCode("jac_f_casadi.c")
