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
import matplotlib.pyplot as plt
import numpy

# Sailboat model based on
#
# [MF2011]:
# Y. Masuyama, Toichi Fukasawa
# "Tacking Simulation of Sailing Yachts with New Model of Aerodynamic
# Force Variation During Tacking Maneuver"
# Journal of Sailboat Technology, Article 2011-01
#
# Joel Andersson, UW Madison 2017
#

# Create DaeBuilder instance
dae = DaeBuilder()

# Physical constants
g = 9.81 # [m/s^2] gravity
rho = 1027. # p[kg/m^3] density of ocean water

# Sailboat model constants. Cf. Table 1 [MF2011]
L = 8.80 # [m] Length of design waterline
D = 2.02 # [m] Design draft, including fin keel
m = 4410. # [kg] Displacement
GM = 1.45 # [m] metacentric height of boat
m_x = 160.; m_y_hull = 2130.; m_y_sail = 280.; m_z = 12000. # [kg] Added masses
Ixx = 17700.; Iyy = 33100.; Izz = 17200. # [kg m^2] Moments of inertia
Jxx_hull = 7200.; Jxx_sail = 8100.; Jyy = 42400.; Jzz = 6700. # [kg m^2] Added moments of inertia
X_pVV = 3.38e-1
X_pPP = 1.40e-3
X_pVVVV = -1.84
X_pT = -1.91e-2
Y_pV = -5.35e-1
Y_pP = -5.89e-3
Y_pVPP = 7.37e-1
Y_pVVP = -5.53e-1
Y_pVVV = 3.07
Y_pP = 2.19e-1
Y_pT = -4.01e-3
K_pV = 2.80e-1
K_pP = 3.36e-3
K_pVPP = -4.07e-1
K_pVVP = 2.24e-1
K_pVVV = -1.38
K_pT = -3.53e-1
N_pV = -3.23e-2
N_pP = -1.52e-2
N_pVPP = 2.71e-4
N_pVVP = -9.06e-2
N_pVVV = -2.98e-2
N_pT = -5.89e-3
C_Xd = -3.79e-2
C_Yd = -1.80e-1
C_Kd = 9.76e-2
C_Nd = 9.74e-2

# States
U = dae.add_x('U') # Velocity along the X axis
V = dae.add_x('V') # Velocity along the Y axis
phi = dae.add_x('phi') # Roll angle
theta = dae.add_x('theta') # Yaw angle

# Auxiliary variables

# Squared boat velocity
V_B2 = U**2 + V**2

# Misc common subexpressions
phi2 = phi**2

# Calculate hydrodynamic forces
beta = dae.add_aux('beta') # Leeway angle
V_p = sin(beta)
V_p2 = V_p**2
V_p3 = V_p2*V_p
V_p4 = V_p2**2
H_fact = 0.5*rho*V_B2*L*D
dae.add_d('X_H', (X_pVV*V_p2 + X_pPP*phi2 + X_pVVVV*V_p4)*H_fact)
dae.add_d('Y_H', (Y_pV*V_p + Y_pP*phi + Y_pVPP*V_p*phi2 + Y_pVVP*V_p2*phi + Y_pVVV*V_p3)*H_fact)
dae.add_d('K_H', (K_pV*V_p + K_pP*phi + K_pVPP*V_p*phi2 + K_pVVP*V_p2*phi + K_pVVV*V_p3)*H_fact*D)
dae.add_d('N_H', (N_pV*V_p + N_pP*phi + N_pVPP*V_p*phi2 + N_pVVP*V_p2*phi + N_pVVV*V_p3)*H_fact*L)
H = dae.add_fun('H', ['phi', 'beta', 'U', 'V'], ['X_H', 'Y_H', 'K_H', 'N_H'])

# Plot it for reference
ngrid_phi = 100; ngrid_beta = 100
U0 = 5.; V0 = 5. # [m/s] Boat speed for simulation
phi0 = numpy.linspace(-pi/4, pi/4, ngrid_phi)
beta0 = numpy.linspace(-pi/4, pi/4, ngrid_beta)
PHI0,BETA0 = numpy.meshgrid(phi0, beta0)
r = H(phi=PHI0, beta=BETA0, U=U0, V=V0)
for i,c in enumerate(['X_H', 'Y_H', 'K_H', 'N_H']):
    plt.subplot(2,2,i+1)
    CS = plt.contour(PHI0*180/pi, BETA0*180/pi, log10(r[c]))
    plt.clabel(CS, inline=1, fontsize=10)
    plt.title('log10(' + c + ')')
    plt.grid(True)


plt.show()
