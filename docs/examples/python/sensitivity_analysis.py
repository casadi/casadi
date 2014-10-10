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
DAE_integrators = ["idas","collocation","oldcollocation"]
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
    ffcn = SXFunction(daeIn(t=t,x=x,p=u),daeOut(ode=ode,quad=quad))

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
    ffcn = SXFunction(daeIn(x=x,z=z,p=u),daeOut(ode=ode,alg=alg,quad=quad))
    
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

    # Integrator
    I = Integrator(MyIntegrator,ffcn)
    
    # Integrator options
    I.setOption("tf",tf)
    if MyIntegrator in ("collocation","oldcollocation"):
      I.setOption("implicit_solver","kinsol")
      I.setOption("implicit_solver_options",{"linear_solver":"csparse"})
      if MyIntegrator=="oldcollocation": I.setOption("expand_f",True)
    I.init()

    # Integrate to get results
    I.setInput(x0,"x0")
    I.setInput(u0,"p")
    I.evaluate()
    xf = I.getOutput("xf")
    qf = I.getOutput("qf")
    print "%50s" % "Unperturbed solution:", "xf  = ", xf, ", qf  = ", qf

    # Perturb solution to get a finite difference approximation
    h = 0.001
    I.setInput(u0+h,"p")
    I.evaluate()
    xf_pert = I.getOutput("xf")
    qf_pert = I.getOutput("qf")
    print "%50s" % "Finite difference approximation:", "d(xf)/d(p) = ", (xf_pert-xf)/h, ", d(qf)/d(p) = ", (qf_pert-qf)/h

    # Calculate once, forward
    I_fwd = I.derivative(1,0)
    I_fwd.setInput(x0,"der_x0")
    I_fwd.setInput(u0,"der_p")
    I_fwd.setInput(0.0,"fwd0_x0")
    I_fwd.setInput(1.0,"fwd0_p")
    I_fwd.evaluate()
    fwd_xf = I_fwd.getOutput("fwd0_xf")
    fwd_qf = I_fwd.getOutput("fwd0_qf")
    print "%50s" % "Forward sensitivities:", "d(xf)/d(p) = ", fwd_xf, ", d(qf)/d(p) = ", fwd_qf

    # Calculate once, adjoint
    I_adj = I.derivative(0,1)
    I_adj.setInput(x0,"der_x0")
    I_adj.setInput(u0,"der_p")
    I_adj.setInput(0.0,"adj0_xf")
    I_adj.setInput(1.0,"adj0_qf")
    I_adj.evaluate()
    adj_x0 = I_adj.getOutput("adj0_x0")
    adj_p = I_adj.getOutput("adj0_p")
    print "%50s" % "Adjoint sensitivities:", "d(qf)/d(x0) = ", adj_x0, ", d(qf)/d(p) = ", adj_p

    # Perturb adjoint solution to get a finite difference approximation of the second order sensitivities
    I_adj.setInput(x0,"der_x0")
    I_adj.setInput(u0+h,"der_p")
    I_adj.setInput(0.0,"adj0_xf")
    I_adj.setInput(1.0,"adj0_qf")
    I_adj.evaluate()
    adj_x0_pert = I_adj.getOutput("adj0_x0")
    adj_p_pert = I_adj.getOutput("adj0_p")
    print "%50s" % "FD of adjoint sensitivities:", "d2(qf)/d(x0)d(p) = ", (adj_x0_pert-adj_x0)/h, ", d2(qf)/d(p)d(p) = ", (adj_p_pert-adj_p)/h

    # Forward over adjoint to get the second order sensitivities
    I_foa = I_adj.derivative(1,0)
    I_foa.setInput(x0,"der_der_x0")
    I_foa.setInput(u0,"der_der_p")
    I_foa.setInput(1.0,"fwd0_der_p")
    I_foa.setInput(0.0,"der_adj0_xf")
    I_foa.setInput(1.0,"der_adj0_qf")
    I_foa.evaluate()
    
    fwd_adj_x0 = I_foa.getOutput("fwd0_adj0_x0")
    fwd_adj_p = I_foa.getOutput("fwd0_adj0_p")
    print "%50s" % "Forward over adjoint sensitivities:", "d2(qf)/d(x0)d(p) = ", fwd_adj_x0, ", d2(qf)/d(p)d(p) = ", fwd_adj_p

    # Adjoint over adjoint to get the second order sensitivities
    I_aoa = I_adj.derivative(0,1)
    I_aoa.setInput(x0,"der_der_x0")
    I_aoa.setInput(u0,"der_der_p")
    I_aoa.setInput(0.0,"der_adj0_xf")
    I_aoa.setInput(1.0,"der_adj0_qf")
    I_aoa.setInput(1.0,"adj0_adj0_p")
    I_aoa.evaluate()
    adj_adj_x0 = I_aoa.getOutput("adj0_der_x0")
    adj_adj_p = I_aoa.getOutput("adj0_der_p")
    print "%50s" % "Adjoint over adjoint sensitivities:", "d2(qf)/d(x0)d(p) = ", adj_adj_x0, ", d2(qf)/d(p)d(p) = ", adj_adj_p

  
