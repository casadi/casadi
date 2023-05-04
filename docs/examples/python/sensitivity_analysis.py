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
# -*- coding: utf-8 -*-

from __future__ import print_function
from casadi import *
from copy import deepcopy
print('Testing sensitivity analysis in CasADi')

# All ODE and DAE integrators to be tested
DAE_integrators = ['idas','collocation']
ODE_integrators = ['cvodes','rk'] + DAE_integrators

for Integrators in (ODE_integrators,DAE_integrators):
  if Integrators==ODE_integrators: # rocket example
    print('******')
    print('Testing ODE example')

    # Time
    t = SX.sym('t')

    # Parameter
    u = SX.sym('u')

    # Differential states
    s = SX.sym('s'); v = SX.sym('v'); m = SX.sym('m')
    x = vertcat(s,v,m)

    # Constants
    alpha = 0.05 # friction
    beta = 0.1   # fuel consumption rate

    # Differential equation
    ode = vertcat(
      v,
      (u-alpha*v*v)/m,
      -beta*u*u)

    # Quadrature
    quad = v**3 + ((3-sin(t)) - u)**2

    # DAE callback function
    dae = {'t':t, 'x':x, 'p':u, 'ode':ode, 'quad':quad}

    # Time length
    tf = 0.5

    # Initial position
    x0 = [0.,0.,1.]

    # Parameter
    u0 = 0.4

  else: # Simple DAE example
    print('******')
    print('Testing DAE example')

    # Differential state
    x = SX.sym('x')

    # Algebraic variable
    z = SX.sym('z')

    # Parameter
    u = SX.sym('u')

    # Differential equation
    ode = -x + 0.5*x*x + u + 0.5*z

    # Algebraic constraint
    alg = z + exp(z) - 1.0 + x

    # Quadrature
    quad = x*x + 3.0*u*u

    # DAE callback function
    dae = {'x':x, 'z':z, 'p':u, 'ode':ode, 'alg':alg, 'quad':quad}

    # End time
    tf = 5.

    # Initial position
    x0 = 1.

    # Parameter
    u0 = 0.4

  # Integrator
  for MyIntegrator in Integrators:
    print('========')
    print('Integrator: ', MyIntegrator)
    print('========')

    # Integrator options
    opts = {}
    if MyIntegrator=='collocation':
      opts['rootfinder'] = 'kinsol'

    # Integrator
    I = integrator('I', MyIntegrator, dae, 0, tf, opts)

    # Integrate to get results
    res = I(x0=x0, p=u0)
    xf = res['xf']
    qf = res['qf']
    print('%50s' % 'Unperturbed solution:', 'xf  = ', xf, ', qf  = ', qf)

    # Perturb solution to get a finite difference approximation
    h = 0.001
    res = I(x0=x0, p=u0+h)
    fd_xf = (res['xf']-xf)/h
    fd_qf = (res['qf']-qf)/h
    print('%50s' % 'Finite difference approximation:', 'd(xf)/d(p) = ', fd_xf, ', d(qf)/d(p) = ', fd_qf)

    # Calculate once, forward
    I_fwd = I.factory('I_fwd', ['x0', 'z0', 'p', 'fwd:p'], ['fwd:xf', 'fwd:qf'])
    res = I_fwd(x0=x0, p=u0, fwd_p=1)
    fwd_xf = res['fwd_xf']
    fwd_qf = res['fwd_qf']
    print('%50s' % 'Forward sensitivities:', 'd(xf)/d(p) = ', fwd_xf, ', d(qf)/d(p) = ', fwd_qf)

    # Calculate once, adjoint
    I_adj = I.factory('I_adj', ['x0', 'z0', 'p', 'adj:qf'], ['adj:x0', 'adj:p'])
    res = I_adj(x0=x0, p=u0, adj_qf=1)
    adj_x0 = res['adj_x0']
    adj_p = res['adj_p']
    print('%50s' % 'Adjoint sensitivities:', 'd(qf)/d(x0) = ', adj_x0, ', d(qf)/d(p) = ', adj_p)

    # Perturb adjoint solution to get a finite difference approximation of the second order sensitivities
    res = I_adj(x0=x0, p=u0+h, adj_qf=1)
    fd_adj_x0 = (res['adj_x0']-adj_x0)/h
    fd_adj_p = (res['adj_p']-adj_p)/h
    print('%50s' % 'FD of adjoint sensitivities:', 'd2(qf)/d(x0)d(p) = ', fd_adj_x0, ', d2(qf)/d(p)d(p) = ', fd_adj_p)

    # Forward over adjoint to get the second order sensitivities
    I_foa = I_adj.factory('I_foa', ['x0', 'z0', 'p', 'adj_qf', 'fwd:p'], ['fwd:adj_x0', 'fwd:adj_p'])
    res = I_foa(x0=x0, p=u0, adj_qf=1, fwd_p=1)
    fwd_adj_x0 = res['fwd_adj_x0']
    fwd_adj_p = res['fwd_adj_p']
    print('%50s' % 'Forward over adjoint sensitivities:', 'd2(qf)/d(x0)d(p) = ', fwd_adj_x0, ', d2(qf)/d(p)d(p) = ', fwd_adj_p)

    # Adjoint over adjoint to get the second order sensitivities
    I_aoa = I_adj.factory('I_aoa', ['x0', 'z0', 'p', 'adj_qf', 'adj:adj_p'], ['adj:x0', 'adj:p'])
    res = I_aoa(x0=x0, p=u0, adj_qf=1, adj_adj_p=1)
    adj_x0 = res['adj_x0']
    adj_p = res['adj_p']
    print('%50s' % 'Adjoint over adjoint sensitivities:', 'd2(qf)/d(x0)d(p) = ', adj_x0, ', d2(qf)/d(p)d(p) = ', adj_p)
