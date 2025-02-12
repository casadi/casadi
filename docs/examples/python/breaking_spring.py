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
import casadi as ca
from matplotlib import pyplot as plt
import numpy as np

# Simulating and optimizing a spring that breaks with DaeBuilder and event handling
# Joel Andersson, 2024-2025

# Start with an empty DaeBuilder instance
dae = ca.DaeBuilder('breaking_spring')

# Model variables
m = dae.add('m', 'parameter', 'tunable', dict(start = 1, description = 'Mass'))
v = dae.add('v', dict(start = -5, description = 'Velocity'))
x = dae.add('x', dict(start = -1, description = 'Displacement'))
k = dae.add('k', 'parameter', 'tunable', dict(start = 2, description = 'Spring constant'))
c = dae.add('c', 'parameter', 'tunable', dict(start = 0.1, description = 'Damping constant'))
d = dae.add('d', 'input', dict(start = 0, description = 'Disturbance'))
f = dae.add('f', dict(description = 'Spring force'))
b = dae.add('b', dict(type = 'Boolean', start = False, description = 'Is the spring broken?'))

# Equations
dae.eq(dae.der(x), v)
dae.eq(f, ca.if_else(ca.logic_not(b), -k*x + d, 0))
dae.eq(dae.der(v), (f - c*v) / m)
dae.when(x > 2, [dae.assign('b', True)])
dae.disp(True)
