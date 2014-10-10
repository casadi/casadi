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
from numpy import *
import matplotlib.pyplot as plt

# CasADi
from casadi import *

# DIFFERENTIAL STATES :
r = SX("r")              #  the length r of the cable
phi = SX("phi")          #  the angle phi
theta = SX("theta")      #  the angle theta
dr = SX("dr")            #  first  derivative of r0    with respect to t
dphi = SX("dphi")        #  first  derivative of phi   with respect to t
dtheta = SX("dtheta")    #  first  derivative of theta with respect to t
n = SX("n")              #  winding number
Psi = SX("Psi")          #  the roll angle Psi
CL = SX("CL")            #  the aerodynamic lift coefficient
W = SX("W")              #  integral over the power at the generator  ( = ENERGY )

# CONTROL :
ddr0 = SX("ddr0")        #  second derivative of r0    with respect to t
dPsi = SX("dPsi")        #  first  derivative of Psi   with respect to t
dCL = SX("dCL")          #  first  derivative of CL    with respect to t

#  PARAMETERS OF THE KITE :
mk =  850.00      #  mass of the kite               #  [ kg    ]
A =  500.00      #  effective area                 #  [ m^2   ]
V =  720.00      #  volume                         #  [ m^3   ]
cd0 =    0.04      #  aerodynamic drag coefficient   #  [       ]
                                     #  ( cd0: without downwash )
K =    0.04      #  induced drag constant          #  [       ]

#   PHYSICAL CONSTANTS :
g =    9.81      #  gravitational constant         #  [ m /s^2]
rho =    1.23      #  density of the air             #  [ kg/m^3]

#  PARAMETERS OF THE CABLE :
rhoc = 1450.00      #  density of the cable           #  [ kg/m^3]
cc = 1.00      #  frictional constant            #  [       ]
dc = 0.05614      #  diameter                       #  [ m     ]

#  PARAMETERS OF THE WIND :
w0 =   10.00      #  wind velocity at altitude h0   #  [ m/s   ]
h0 =  100.00      #  the altitude h0                #  [ m     ]
hr =    0.10      #  roughness length               #  [ m     ]

# DEFINITION OF PI :
PI = 3.1415926535897932

# CROSS AREA OF THE CABLE :
AQ      =  PI*dc*dc/4.0                                       

# THE EFECTIVE MASS' :
mc      =  rhoc*AQ*r           # mass of the cable
m       =  mk + mc     / 3.0   # effective inertial mass
m_      =  mk + mc     / 2.0   # effective gravitational mass
dm      =  (rhoc*AQ/ 3.0)*dr   # time derivative of the mass

# WIND SHEAR MODEL :
h       =  r*cos(theta)                    #  altitude of the kite                   
w       =  log(h/hr) / log(h0/hr) *w0      #  the wind at altitude h 

#  effective wind vector
we = [w*sin(theta)*cos(phi) - dr, -w*sin(phi)- r*sin(theta)*dphi, -w*cos(theta)*cos(phi) + r*dtheta]

nwep    =  pow(                we[1]*we[1] + we[2]*we[2], 0.5 )
nwe     =  pow(  we[0]*we[0] + we[1]*we[1] + we[2]*we[2], 0.5 )
eta     =  arcsin( we[0]*tan(Psi)/ nwep ) #  angle eta: cf. [2]                      

#  projection of ewe
ewep = [0.0, we[1] / nwep, we[2] / nwep]                                       

#  unit vector from the left to the right wing tip
et = [sin(Psi), (-cos(Psi)*sin(eta))*ewep[1] - (cos(Psi)*cos(eta))*ewep[2], (-cos(Psi)*sin(eta))*ewep[2] + (cos(Psi)*cos(eta))*ewep[1]]

# CONSTANTS FOR SEVERAL FORCES :
Cg      =  (V*rho-m_)*g                                       
Caer    =  (rho*A/2.0 )*nwe                                   
Cf      =  (rho*dc/8.0)*r*nwe              #  simple constants                   

#  the aerodynamic drag coefficient
CD      =  cd0 + K*CL*CL                                      

#  the gravitational force
Fg = [Cg  *  cos(theta), Cg  *  0.0, Cg  *  sin(theta)]

#  the aerodynamic force
Faer = [Caer*( CL*(we[1]*et[2]-we[2]*et[1]) + CD*we[0] ), Caer*( CL*(we[2]*et[0]-we[0]*et[2]) + CD*we[1] ), Caer*( CL*(we[0]*et[1]-we[1]*et[0]) + CD*we[2] )]

#  the frictional force
Ff = [Cf  *  cc* we[0], Cf  *  cc* we[1], Cf  *  cc* we[2]]

#  the total force
F = [Fg[0] + Faer[0], Faer[1] + Ff[1], Fg[2] + Faer[2] + Ff[2]]

#  the pseudo accelaration
a_pseudo = [
  - ddr0 + r*(dtheta*dtheta + sin(theta)*sin(theta)*dphi  *dphi) - dm/m*dr,
  - 2.0*cos(theta)/sin(theta)*dphi*dtheta - 2.0*dr/r*dphi - dm/m*dphi,
  cos(theta)*sin(theta)*dphi*dphi - 2.0*dr/r*dtheta - dm/m*dtheta
  ]

#  second derivative of s     with respect to t
ddr =  F[0]/m + a_pseudo[0]           

#  second derivative of phi   with respect to t
ddphi        =  F[1]/(m*r*sin(theta)) + a_pseudo[1]

#  second derivative of theta with respect to t
ddtheta      = -F[2]/(m*r           ) + a_pseudo[2]

#  the time derivate of the kite's winding number
dn =  (dphi*ddtheta - dtheta*ddphi     ) / (2.0*PI*(dphi*dphi    + dtheta*dtheta)   )

#  the power at the generator
power        =   m*ddr*dr

#  regularisation of controls
regularisation = 5.0e2 * ddr0 * ddr0 + 1.0e8 * dPsi * dPsi + 1.0e5 * dCL * dCL + 2.5e5 * dn * dn + 2.5e7 * ddphi   * ddphi + 2.5e7 * ddtheta * ddtheta + 2.5e6 * dtheta  * dtheta

# Differential state vector
x = [r, phi, theta, dr, dphi, dtheta, n, Psi, CL, W]

# Control
u = [ddr0, dPsi, dCL]

# The right hand side of the ACADO functions
acado_in = ACADO_FCN_NUM_IN * [[]]
acado_in[ACADO_FCN_T]  = [ ]  # Time
acado_in[ACADO_FCN_XD] =  x   # Differential states
acado_in[ACADO_FCN_XA] = [ ]  # Algebraic state
acado_in[ACADO_FCN_U]  =  u   # Control
acado_in[ACADO_FCN_P]  = [ ]  # Parameter

# ODE rhs
f = [dr, dphi, dtheta, ddr0, ddphi, ddtheta, dn, dPsi, dCL, (-power + regularisation)*1.0e-6]
ffcn = SXFunction(acado_in,[f])
f
# Objective function
mfcn = SXFunction(acado_in,[[W]])

# Path constraint function
cfcn = SXFunction(acado_in,[[ddr0, dPsi, dCL]])

# Time horizon
t0 = 0.0
tf = 18.0

# Number of shooting nodes
num_nodes = 18

# Create ACADO solver
ocp_solver = AcadoInterface(ffcn,mfcn,cfcn)

# Set options
ocp_solver.setOption("start_time",t0)
ocp_solver.setOption("final_time",tf)
ocp_solver.setOption("number_of_shooting_nodes",num_nodes)
#ocp_solver.setOption("max_num_iterations",30)
#ocp_solver.setOption("max_num_integrator_steps",10000)
#ocp_solver.setOption("dynamic_sensitivity","forward_sensitivities")
ocp_solver.setOption("kkt_tolerance",1e-5)
#ocp_solver.setOption("absolute_tolerance",1e-4)
#ocp_solver.setOption("integrator_tolerance",1e-8)
#ocp_solver.setOption("auto_init",True)
#ocp_solver.setOption("relaxation_parameter",1.5)
ocp_solver.setOption("periodic_bounds",[1, 1, 1, 1, 1, 1, 0, 1, 1, 0])

# Initialize
ocp_solver.init()

lbx0 = [-inf, -inf, -inf, -inf, -inf, -inf, 0, -inf, -inf, 0]
ubx0 = [ inf,  inf,  inf,  inf,  inf,  inf, 0,  inf,  inf, 0]
ocp_solver.setInput(lbx0, "lbx0")
ocp_solver.setInput(ubx0, "ubx0")

lbx = [-inf, -0.34, 0.85, -40.0, -inf, -inf, -0.7, -0.29, 0.1, -inf]
ubx = [ inf,  0.34, 1.45,  10.0,  inf,  inf,  0.9,  0.29, 1.5,  inf]
ocp_solver.setInput(lbx, "lbx")
ocp_solver.setInput(ubx, "ubx")
  
lbc = [-25.0, -0.065, -3.5]
ubc = [ 25.0,  0.065,  3.5]
ocp_solver.setInput(lbc, "lbc")
ocp_solver.setInput(ubc, "ubc")

# Solve
ocp_solver.solve()


