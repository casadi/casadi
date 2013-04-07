#
#     This file is part of CasADi.
# 
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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
DAE_integrators = [CollocationIntegrator,IdasIntegrator]
ODE_integrators = DAE_integrators + [CVodesIntegrator]

for Integrators in (ODE_integrators,DAE_integrators):    
  if Integrators==ODE_integrators: # rocket example
    print "******"
    print "Testing ODE example"
    
    # Time 
    t = ssym("t")
      
    # Parameter
    u = ssym("u")

    # Differential states
    s = ssym("s"); v = ssym("v"); m = ssym("m")
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
    x = ssym("x")
    
    # Algebraic variable
    z = ssym("z")
    
    # Parameter
    u = ssym("u")
    
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
    print "Integrator: ", MyIntegrator.__name__
    print "========"

    # Integrator
    I = MyIntegrator(ffcn)
    
    # Integrator options
    I.setOption("tf",tf)
    if MyIntegrator==CollocationIntegrator:
      I.setOption("implicit_solver",KinsolSolver)
      I.setOption("implicit_solver_options",{"linear_solver":CSparse})
      I.setOption("expand_f",True)
    I.init()

    # Integrate to get results
    I.setInput(x0,"x0")
    I.setInput(u0,"p")
    I.evaluate()
    xf = deepcopy(I.output("xf"))
    qf = deepcopy(I.output("qf"))
    print "%50s" % "Unperturbed solution:", "xf  = ", xf, ", qf  = ", qf

    # Perturb solution to get a finite difference approximation
    h = 0.001
    I.setInput(u0+h,"p")
    I.evaluate()
    xf_pert = deepcopy(I.output("xf"))
    qf_pert = deepcopy(I.output("qf"))
    print "%50s" % "Finite difference approximation:", "d(xf)/d(p) = ", (xf_pert-xf)/h, ", d(qf)/d(p) = ", (qf_pert-qf)/h

    # Operator overloading approach
    I.setInput(x0,"x0")
    I.setInput(u0,"p")
    I.setFwdSeed(DMatrix.zeros(x.shape),"x0")
    I.setFwdSeed(1.0,"p")
    I.reset(1,0,0)
    I.integrate(tf)
    oo_xf = deepcopy(I.fwdSens("xf"))
    oo_qf = deepcopy(I.fwdSens("qf"))
    print "%50s" % "Forward sensitivities via OO:", "d(xf)/d(p) = ", oo_xf, ", d(qf)/d(p) = ", oo_qf

    # Calculate once, forward
    I_fwd = I.derivative(1,0)
    
    I_fwd.setInput(x0,INTEGRATOR_X0)
    I_fwd.setInput(u0,INTEGRATOR_P)
    I_fwd.setInput(DMatrix.zeros(x.shape),INTEGRATOR_NUM_IN+INTEGRATOR_X0)
    I_fwd.setInput(1.0,INTEGRATOR_NUM_IN+INTEGRATOR_P)
    I_fwd.evaluate()
    fwd_xf = deepcopy(I_fwd.output(INTEGRATOR_NUM_OUT+INTEGRATOR_XF))
    fwd_qf = deepcopy(I_fwd.output(INTEGRATOR_NUM_OUT+INTEGRATOR_QF))
    print "%50s" % "Forward sensitivities:", "d(xf)/d(p) = ", fwd_xf, ", d(qf)/d(p) = ", fwd_qf

    # Calculate once, adjoint
    I_adj = I.derivative(0,1)
    I_adj.setInput(x0,INTEGRATOR_X0)
    I_adj.setInput(u0,INTEGRATOR_P)
    I_adj.setInput(DMatrix.zeros(x.shape),INTEGRATOR_NUM_IN+INTEGRATOR_XF)
    I_adj.setInput(1.0,INTEGRATOR_NUM_IN+INTEGRATOR_QF)
    I_adj.evaluate()
    adj_x0 = deepcopy(I_adj.output(INTEGRATOR_NUM_OUT+INTEGRATOR_X0))
    adj_p = deepcopy(I_adj.output(INTEGRATOR_NUM_OUT+INTEGRATOR_P))
    print "%50s" % "Adjoint sensitivities:", "d(qf)/d(x0) = ", adj_x0, ", d(qf)/d(p) = ", adj_p

    # Perturb adjoint solution to get a finite difference approximation of the second order sensitivities
    I_adj.setInput(x0,INTEGRATOR_X0)
    I_adj.setInput(u0+h,INTEGRATOR_P)
    I_adj.setInput(DMatrix.zeros(x.shape),INTEGRATOR_NUM_IN+INTEGRATOR_XF)
    I_adj.setInput(1.0,INTEGRATOR_NUM_IN+INTEGRATOR_QF)
    I_adj.evaluate()
    adj_x0_pert = deepcopy(I_adj.output(INTEGRATOR_NUM_OUT+INTEGRATOR_X0))
    adj_p_pert = deepcopy(I_adj.output(INTEGRATOR_NUM_OUT+INTEGRATOR_P))
    print "%50s" % "FD of adjoint sensitivities:", "d2(qf)/d(x0)d(p) = ", (adj_x0_pert-adj_x0)/h, ", d2(qf)/d(p)d(p) = ", (adj_p_pert-adj_p)/h

    # Forward over adjoint to get the second order sensitivities
    I_adj.setInput(x0,INTEGRATOR_X0)
    I_adj.setInput(u0,INTEGRATOR_P)
    I_adj.setFwdSeed(1.0,INTEGRATOR_P)
    I_adj.setInput(DMatrix.zeros(x.shape),INTEGRATOR_NUM_IN+INTEGRATOR_XF)
    I_adj.setInput(1.0,INTEGRATOR_NUM_IN+INTEGRATOR_QF)
    I_adj.evaluate(1,0)
    fwd_adj_x0 = deepcopy(I_adj.fwdSens(INTEGRATOR_NUM_OUT+INTEGRATOR_X0))
    fwd_adj_p = deepcopy(I_adj.fwdSens(INTEGRATOR_NUM_OUT+INTEGRATOR_P))
    print "%50s" % "Forward over adjoint sensitivities:", "d2(qf)/d(x0)d(p) = ", fwd_adj_x0, ", d2(qf)/d(p)d(p) = ", fwd_adj_p

    # Adjoint over adjoint to get the second order sensitivities
    I_adj.setInput(x0,INTEGRATOR_X0)
    I_adj.setInput(u0,INTEGRATOR_P)
    I_adj.setInput(DMatrix.zeros(x.shape),INTEGRATOR_NUM_IN+INTEGRATOR_XF)
    I_adj.setInput(1.0,INTEGRATOR_NUM_IN+INTEGRATOR_QF)
    I_adj.setAdjSeed(1.0,INTEGRATOR_NUM_OUT+INTEGRATOR_P)
    I_adj.evaluate(0,1)
    adj_adj_x0 = deepcopy(I_adj.adjSens(INTEGRATOR_X0))
    adj_adj_p = deepcopy(I_adj.adjSens(INTEGRATOR_P))
    print "%50s" % "Adjoint over adjoint sensitivities:", "d2(qf)/d(x0)d(p) = ", adj_adj_x0, ", d2(qf)/d(p)d(p) = ", adj_adj_p

  
