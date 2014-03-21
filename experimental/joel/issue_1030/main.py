from __future__ import division
from casadi import *
from numpy import *
from pylab import *
import time
import cProfile
import os
import subprocess

#===============================================================================
# Here you set if you are running og the 1.7.1 or the 1.8.0-release.
runWith171 = 0
runWith180 = 0
runWith190 = 1
calc_ic = True

if runWith190:
  Mat = SX
else:
  Mat = SXMatrix

# This is a model of a small power network, and therefore there are imaginary numbers in the model.
# This has been solved by teating an imaginary number a = a_re + j*a_im as
# a = [a_re, a_im]. In AXMatrix_functions.py a number of mathematical operations with imaginary numbers has been adopted to
# this way of defining them.
execfile("SXMatrix_functions.py")

#===============================================================================
# Constants
y12 = 1-1j
y13 = 1-1j
y23 = 1-1j

Y = array(([[y12+y13, -y12, -y13], [-y12, y12+y13, -y13], [-y13, -y23, y13 + y23]]))
Y_complex = [real(Y), imag(Y)]

Il_abs = array(([[0], [0], [-0.7774]]))
Yl = array(([[0,0,0],[0,0,0],[0,0,0+1j*0.3022]]))
Yl_complex = [real(Yl), imag(Yl)]

Ya = diag([1j,1j,0])
Ya_complex = [real(Ya), imag(Ya)]
Yb = array(([1j,0],[0,1j],[0,0]))
Yb_complex = [real(Yb), imag(Yb)]

Ylarge = dot(linalg.inv(Ya + Y + Yl), Yb)
Ylarge_complex = [real(Ylarge), imag(Ylarge)]

E_abs = array(([[1.2943], [1.1609]]))
omega_s = 2*pi*50
H = array(([[3], [3.5]]))
Dg = array(([[3], [3]]))


#===============================================================================
# Function
# Declare variables
t = Mat.sym("t")                  # Time
Pm = Mat.sym("u", 2)              # control
xdyn = Mat.sym("xd",4)            # Differential state
xalg = Mat.sym("xa",6)            # Algebraic states

delta = xdyn[0:2]
omega_delta = xdyn[2:4]

U = xalg[0:3]
theta = xalg[3:6]

U_complex = imexp(U,theta)
E_complex = imexp(E_abs,delta)

#==============================
#Algebraic equation

# Ig from Kirchoff
IlA_complex = imexp(Il_abs,theta)
IlB_complex = immul(Yl_complex,U_complex)
IL_complex = [IlA_complex[0] + IlB_complex[0], IlA_complex[1] + IlB_complex[1]]
Isys_complex = immul(Y_complex, U_complex)

Ig_complex = [-IL_complex[0]-Isys_complex[0], -IL_complex[1]-Isys_complex[1]]

# Ig from (U-E)/r. r=1
IgB1_complex = immul(Ya_complex,U_complex)
IgB2_complex = immul(Yb_complex,E_complex)

IgB_complex = [IgB1_complex[0]-IgB2_complex[0], IgB1_complex[1]-IgB2_complex[1]]

# Algebraic equation. Ig_complex = IgB_complex. (Divided into real and imaginary part.)
alg = vstack([Ig_complex[0] - IgB_complex[0] ,\
              Ig_complex[1] - IgB_complex[1]])

#===============================
# Dynamic equation

# Finds Pe to use in dynamic equations
Se = imelemmul(E_complex,imconj([Ig_complex[0][0:2], Ig_complex[1][0:2]]))
Pe = Se[0]

# Dynamic equation. (Swing equation of generators.)
ode = vstack([ omega_delta*omega_s, \
               (0.5/H)*(Pm-Pe-Dg*(omega_delta))])

#===============================
# Creates SXFunction
f = SXFunction(daeIn(x=xdyn, z=xalg, p=Pm), daeOut(ode=ode, alg=alg))
f.init()

#==============================================================================
# Initial values

# Input
Pm0 = array(([[0.3054], [0.5]]))

# Algebraic states
U_abs0  = array(([[1], [1], [1.2863]]))
theta0 = array(([[0], [-0.0617], [0.0665]]))

Z0 = [U_abs0[0,0], U_abs0[1,0], U_abs0[2,0], theta0[0,0], theta0[1,0], theta0[2,0]]

# Dynamic states
delta0 = array(([[0.2382], [0.3836]]))
omega_delta0 = zeros((2,1))

X0 = vstack([delta0, omega_delta0])

#===============================================================================
#Simulates

# Creates IdasIntegrator
I = IdasIntegrator(f)
#I.setOption("linear_solver",CSparse)
#I.setOption("linear_solver_type","user_defined")

if runWith171==1:
    I.setOption("init_z",Z0)
I.setOption("calc_ic", calc_ic)
I.init()

I.setInput(X0,INTEGRATOR_X0)
I.setInput(Pm0,INTEGRATOR_P)
if runWith180 or runWith190:
    I.input(INTEGRATOR_Z0).set(Z0)
I.evaluate()

# Simulate
ts = numpy.linspace(0,60,100) # (start, stop, points)
sim = Simulator(I,ts)
sim.init()
sim.input(INTEGRATOR_X0).set(X0)
sim.input(INTEGRATOR_P).set(Pm0)
if runWith180 or runWith190:
    sim.input(INTEGRATOR_Z0).set(Z0)
sim.evaluate()

sol = sim.output().toArray()
if runWith190:
  sol = sol.T

figure(1)
plot(ts,sol[:,2],ts,sol[:,3])

show()
close()



















