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
  xp0 = create_symbolic("x",len(x))
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

  alpha = 0.05 # friction
  beta = 0.1 # fuel consumption rate
  
  # Differential equation
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
  ffcn.setOption("ad_order",1)

  # Explicit integrator (CVODES)
  integrator = CVodesIntegrator(ffcn)
  
  # Set options
  #integrator.setOption("exact_jacobian",True)
  #integrator.setOption("linear_multistep_method","bdf") # adams or bdf
  #integrator.setOption("nonlinear_solver_iteration","newton") # newton or functional
  integrator.setOption("ad_order",1)
  integrator.setOption("fsens_err_con",True)
  integrator.setOption("quad_err_con",True)
  integrator.setOption("abstol",1e-6)
  integrator.setOption("reltol",1e-6)
  #integrator.setOption("fsens_all_at_once",False)

  return integrator



# Main function

# Time length
T = 10.0

# Shooting length
nu = 1000 # Number of control segments
DT = double(T)/nu

# Initial position, speed and mass
s0 = 0 # initial position
v0 = 0 # initial speed
m0 = 1 # initial mass
X0 = [s0,v0,m0]

# Create integrators
integrator_euler = create_integrator_euler()
integrator_cvodes = create_integrator_cvodes()

for integrator in [integrator_euler, integrator_cvodes]:
  # Enable AD
  integrator.setOption("ad_order",1)

  # Initialize the integrator
  integrator.init()

  # control for all segments
  U = MX("U",nu)

  # Dummy input corresponding to the state derivative
  xdot = MX([0,0,0])

  # Integrate over all intervals
  X=MX(X0)
  for k in range(nu):
    X = integrator([MX(k*DT),MX((k+1)*DT),X,U[k],xdot])  # build up a graph with function calls

  # Objective function
  F = inner_prod(U,U)

  # Terminal constraints
  G = vertcat(X[0],X[1])
  
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
  solver.setInput(Umin,NLP_LBX)

  Umax = nu * [10]  # upper bound
  solver.setInput(Umax,NLP_UBX)

  Usol = nu * [0.4] # initial guess
  solver.setInput(Usol,NLP_X_INIT)

  # Bounds on g
  Gmin = Gmax = [10, 0]
  solver.setInput(Gmin,NLP_LBG)
  solver.setInput(Gmax,NLP_UBG)

  # Solve the problem
  solver.solve()

  # Get the solution
  xopt = solver.getOutput(NLP_X_OPT)

  # Plot the optimal trajectory
  plt.figure()
  plt.plot(xopt) 
  plt.show()
  
  