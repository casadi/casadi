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
from casadi import *
import pylab as plt

# Height and velocity of ball
h = SX.sym('h')
v = SX.sym('v')
x = vertcat(h, v)

# ODE right-hand-side
hdot = v
vdot = -9.81
xdot = vertcat(hdot, vdot)

# Event indicator, trigger when crossing zero
event_indicator = -h

# DAE problem structure, with zero-crossing output
dae = dict(x = x, ode = xdot, zero = event_indicator)

# Event transition function
post_x = vertcat(h, -0.8*v)
event_transition = Function('event_transition', dict(x = x, post_x = post_x),
                            event_in(), event_out())

# Create an integrator instance for integrating over 7s
tgrid = np.linspace(0, 7, 100)
sim = integrator('sim', 'cvodes', dae, 0, tgrid,
                 dict(event_transition = event_transition))

# Simulate with initial height of 5
x0 = [5, 0]
simres = sim(x0 = x0)

# Visualize the solution
plt.figure(1)
plt.clf()
plt.plot(tgrid, simres['xf'][0, :].T)
plt.grid()
plt.show()
