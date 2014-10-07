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
from numpy import *
import matplotlib.pyplot as plt

# Construct an explicit euler integrator
def create_integrator_euler():
  u = SX("u") # control for one segment

  # Initial position
  s0 = SX("s0") # initial position
  v0 = SX("v0") # initial speed
  m0 = SX("m0") # initial mass

  t0 = SX("t0") # initial time
  tf = SX ("tf") # final time

  nj = 1000 # Number of integration steps per control segment
  dt = (tf-t0)/nj # time step
  alpha = SX(0.05)    # friction
  beta = SX(0.1)      # fuel consumption rate

  # Integrate over the interval with Euler forward
  s = s0;   v = v0;    m = m0

  dm = -dt*beta*u*u
  for j in range(nj):
    s += dt*v
    v += dt * (u-alpha*v*v)/m
    m += dm

  # State vector
  x =  [s,v,m]
  x0 = [s0,v0,m0];

  # Input to the integrator function being created
  integrator_in = INTEGRATOR_NUM_IN * [[]]
  integrator_in[INTEGRATOR_T0] = [t0]
  integrator_in[INTEGRATOR_TF] = [tf]
  integrator_in[INTEGRATOR_X0] =  x0
  integrator_in[INTEGRATOR_P]  =  [u]

  # Create a dummy state derivative vector
  xp0 = list(ssym("x",len(x)))
  integrator_in[INTEGRATOR_XP0] = xp0

  # Create the explicit Euler integrator
  integrator = SXFunction(integrator_in,[x])
  return integrator

# Construct an explicit multistep integrator from Sundials
def create_integrator_cvodes():
  # Time 
  t = SX("t")

  # Differential states
  s = SX("s"); v = SX("v"); m = SX("m")
  y = [s,v,m]

  # Control
  u = SX("u")

  # Parameters
  alpha = 0.05 # friction
  beta = 0.1 # fuel consumption rate

  # Differential equation - explicit form
  sdot = v
  vdot = (u-alpha*v*v)/m
  mdot = -beta*u*u
  rhs = [sdot,vdot,mdot]

  # Initial conditions
  y0 = [0,0,1]
  
  # Input of the ode rhs
  ffcn_in = ODE_NUM_IN * [[]]
  ffcn_in[ODE_T] = [t]
  ffcn_in[ODE_Y] =  y
  ffcn_in[ODE_P] = [u]
  
  # ODE right hand side
  ffcn = SXFunction(ffcn_in,[rhs])
  ffcn.setOption("name","ODE right hand side")

  # Explicit integrator (CVODES)
  integrator = CVodesIntegrator(ffcn)
  
  # Set options
  #integrator.setOption("exact_jacobian",True)
  #integrator.setOption("linear_multistep_method","bdf") # adams or bdf
  #integrator.setOption("nonlinear_solver_iteration","newton") # newton or functional
  integrator.setOption("fsens_err_con",True)
  integrator.setOption("abstol",1e-6)
  integrator.setOption("reltol",1e-6)
  #integrator.setOption("fsens_all_at_once",False)

  return integrator

# Create an IDAS instance (fully implicit integrator)
def create_integrator_idas():
  # Time 
  t = SX("t")

  # Differential states
  s = SX("s")
  v = SX("v")
  m = SX("m")
  y = [s,v,m]
  
  # State derivatives
  sdot = SX("sdot")
  vdot = SX("vdot")
  mdot = SX("mdot")
  ydot = [sdot,vdot,mdot]

  # Control
  u = SX("u")

  # Parameters
  alpha = 0.05 # friction
  beta = 0.1 # fuel consumption rate

  # Differential equation (fully implicit form)
  sres = v - sdot
  vres = (u-alpha*v*v)/m - vdot
  mres = -beta*u*u - mdot
  res = [sres, vres, mres]

  # Input of the DAE residual function
  ffcn_in = DAE_NUM_IN * [[]]
  ffcn_in[DAE_T] = [t]
  ffcn_in[DAE_Y] = y
  ffcn_in[DAE_YDOT] = ydot
  ffcn_in[DAE_P] = [u]

  # DAE residual function
  ffcn = SXFunction(ffcn_in,[res])
  ffcn.setOption("name","DAE residual")
  
  # Create an integrator
  integrator = IdasIntegrator(ffcn)

  # Set options
  integrator.setOption("calc_ic",True)
  integrator.setOption("is_differential",[1,1,1])
  integrator.setOption("fsens_err_con",True)
  integrator.setOption("abstol",1e-6)
  integrator.setOption("reltol",1e-6)
  integrator.setOption("steps_per_checkpoint",100)

  return integrator


# Main function

# Time length
T = 10.0

# Shooting length
nu = 20 # Number of control segments
DT = double(T)/nu

# Initial position, speed and mass
s0 = 0 # initial position
v0 = 0 # initial speed
m0 = 1 # initial mass
X0 = [s0,v0,m0]

# Create integrators
integrator_euler = create_integrator_euler()
integrator_cvodes = create_integrator_cvodes()
#integrator_idas = create_integrator_idas()

for integrator in [integrator_euler, integrator_cvodes]:
  # Enable AD

  # Initialize the integrator
  integrator.init()

  # control for all segments
  U = MX("U",nu)

  # Dummy input corresponding to the state derivative
  xdot = MX([0,0,0])

  # Integrate over all intervals
  X=MX(X0)
  T0 = MX(0) # Beginning of time interval (changed from k*DT due to probable Sundials bug)
  TF = MX(DT) # End of time interval (changed from (k+1)*DT due to probable Sundials bug)
  for k in range(nu):
    # build up a graph with function calls
    X = integrator([T0,TF,X,U[k],xdot,MX()])

  # Objective function
  F = inner_prod(U,U)

  # Terminal constraints
  G = vertcat((X[0],X[1]))
  
  # Create the NLP
  ffcn = MXFunction([U],[F]) # objective function
  gfcn = MXFunction([U],[G]) # constraint function

  # Allocate an NLP solver
  #solver = LiftedNewtonSolver(ffcn,gfcn)
  solver = IpoptSolver(ffcn,gfcn)
  
  # Set options
  solver.setOption("tol",1e-10)
  solver.setOption("hessian_approximation","limited-memory");

  # initialize the solver
  solver.init()

  # Bounds on u and initial condition
  Umin = nu * [-10] # lower bound
  solver.setInput(Umin,"lbx")

  Umax = nu * [10]  # upper bound
  solver.setInput(Umax,"ubx")

  Usol = nu * [0.4] # initial guess
  solver.setInput(Usol,"x0")

  # Bounds on g
  Gmin = Gmax = [10, 0]
  solver.setInput(Gmin,"lbg")
  solver.setInput(Gmax,"ubg")

  # Solve the problem
  solver.solve()

  # Get the solution
  xopt = solver.getOutput("x")

  # Plot the optimal trajectory
  plt.figure()
  plt.plot(xopt) 
  plt.show()
  
  
