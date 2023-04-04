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
from casadi import *

# Example on how to use the DaeBuilder class
# Joel Andersson, UW Madison 2017

# Start with an empty DaeBuilder instance
dae = DaeBuilder('rocket')

# Add input expressions
a = dae.add_p('a')
b = dae.add_p('b')
u = dae.add_u('u')
h = dae.add_x('h')
v = dae.add_x('v')
m = dae.add_x('m')

# Constants
g = 9.81 # gravity

# Set ODE right-hand-side
dae.set_ode('h', v)
dae.set_ode('v', (u-a*v**2)/m-g)
dae.set_ode('m', -b*u**2)

# Specify initial conditions
dae.set_start('h', 0)
dae.set_start('v', 0)
dae.set_start('m', 1)

# Add meta information
dae.set_unit('h','m')
dae.set_unit('v','m/s')
dae.set_unit('m','kg')

# Print DAE
dae.disp(True)
