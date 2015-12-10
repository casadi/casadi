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
from copy import deepcopy
print "Testing sensitivity analysis in CasADi"

# All ODE and DAE integrators to be tested
DAE_integrators = ["idas","collocation"]
ODE_integrators = ["cvodes","rk"] + DAE_integrators

for Integrators in (ODE_integrators,DAE_integrators):    
  if Integrators==ODE_integrators: # rocket example
    print "******"
    print "Testing ODE example"
    
    # Time 
    t = SX.sym("t")
      
    # Parameter
    u = SX.sym("u")

    # Differential states
    s = SX.sym("s"); v = SX.sym("v"); m = SX.sym("m")
    x = vertcat([s,v,m])

    # Constants
    alpha = 0.05 # friction
    beta = 0.1   # fuel consumption rate
      
    # Differential equation
    ode = vertcat([
      v,
      (u-alpha*v*v)/m,
      -beta*u*u])
      
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
    print "******"
    print "Testing DAE example"
    
    # Differential state
    x = SX.sym("x")
    
    # Algebraic variable
    z = SX.sym("z")
    
    # Parameter
    u = SX.sym("u")
    
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
    print "========"
    print "Integrator: ", MyIntegrator
    print "========"
    
    # Integrator options
    opts = {"tf":tf}
    if MyIntegrator=="collocation":
      opts["implicit_solver"] = "kinsol"
      opts["implicit_solver_options"] = {"linear_solver":"csparse"}

    # Integrator
    I = integrator("I", MyIntegrator, dae, opts)

    # Integrate to get results
    arg = {"x0":x0, "p":u0}
    res = I(arg)
    xf = res["xf"]
    qf = res["qf"]
    print "%50s" % "Unperturbed solution:", "xf  = ", xf, ", qf  = ", qf

    # Perturb solution to get a finite difference approximation
    h = 0.001
    arg["p"] = u0+h
    res = I(arg)
    fd_xf = (res["xf"]-xf)/h
    fd_qf = (res["qf"]-qf)/h
    print "%50s" % "Finite difference approximation:", "d(xf)/d(p) = ", fd_xf, ", d(qf)/d(p) = ", fd_qf

    # Calculate once, forward
    I_fwd = I.derivative(1, 0)
    arg = {}
    arg["der_x0"] = x0
    arg["der_p"] = u0
    arg["fwd0_x0"] = 0
    arg["fwd0_p"] = 1
    res = I_fwd(arg)
    fwd_xf = res["fwd0_xf"]
    fwd_qf = res["fwd0_qf"]
    print "%50s" % "Forward sensitivities:", "d(xf)/d(p) = ", fwd_xf, ", d(qf)/d(p) = ", fwd_qf

    # Calculate once, adjoint
    I_adj = I.derivative(0, 1)
    arg = {}
    arg["der_x0"] = x0
    arg["der_p"] = u0
    arg["adj0_xf"] = 0
    arg["adj0_qf"] = 1
    res = I_adj(arg)
    adj_x0 = res["adj0_x0"]
    adj_p = res["adj0_p"]
    print "%50s" % "Adjoint sensitivities:", "d(qf)/d(x0) = ", adj_x0, ", d(qf)/d(p) = ", adj_p

    # Perturb adjoint solution to get a finite difference approximation of the second order sensitivities
    arg["der_p"] = u0+h
    res = I_adj(arg)
    fd_adj_x0 = (res["adj0_x0"]-adj_x0)/h
    fd_adj_p = (res["adj0_p"]-adj_p)/h
    print "%50s" % "FD of adjoint sensitivities:", "d2(qf)/d(x0)d(p) = ", fd_adj_x0, ", d2(qf)/d(p)d(p) = ", fd_adj_p

    # Forward over adjoint to get the second order sensitivities
    I_foa = I_adj.derivative(1, 0)
    arg = {}
    arg["der_der_x0"] = x0
    arg["der_der_p"] = u0
    arg["fwd0_der_p"] = 1
    arg["der_adj0_xf"] = 0
    arg["der_adj0_qf"] = 1
    res = I_foa(arg)
    fwd_adj_x0 = res["fwd0_adj0_x0"]
    fwd_adj_p = res["fwd0_adj0_p"]
    print "%50s" % "Forward over adjoint sensitivities:", "d2(qf)/d(x0)d(p) = ", fwd_adj_x0, ", d2(qf)/d(p)d(p) = ", fwd_adj_p

    # Adjoint over adjoint to get the second order sensitivities
    I_aoa = I_adj.derivative(0, 1)
    arg = {}
    arg["der_der_x0"] = x0
    arg["der_der_p"] = u0
    arg["der_adj0_xf"] = 0
    arg["der_adj0_qf"] = 1
    arg["adj0_adj0_p"] = 1
    res = I_aoa(arg)
    adj_adj_x0 = res["adj0_der_x0"]
    adj_adj_p = res["adj0_der_p"]
    print "%50s" % "Adjoint over adjoint sensitivities:", "d2(qf)/d(x0)d(p) = ", adj_adj_x0, ", d2(qf)/d(p)d(p) = ", adj_adj_p

  
