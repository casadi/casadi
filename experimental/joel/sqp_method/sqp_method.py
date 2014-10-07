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
from casadi import *
from numpy import *
import matplotlib.pyplot as plt

# The algorithm is a Quasi-Newton method with damped BFGS updating to that
# assures positive definitenes of the Hessian approximation. Line search is
# carried out via backtracking until with the Armijo condition applied to
# the T1 (in Nocedal phi1) merit function is satisfied.

# Controls
u = ssym("u",10)

# Starting position
s1=0
v1=0
m1=1

# Final position
s100=SX(10)
v100=SX(0)

# Time step
dt=SX(0.1)
alpha=SX(0.05)
beta=SX(0.1)

# Constraints
vmax=1.3

# starting point
u0=0.4*ones(10)
umin = -1*ones(10)
umax = 0.5*ones(10)

# minimize fuel
f = inner_prod(u,u)/2

# Calculate acceleration for all points
a = SX(100,1,0)
for i in range(10):
  a[10*i:10*(i+1)] = u[i]

# Loop over all k to get s and v and the endpoints
s = SX(101,1,0)
s[0] = s1
v = SX(101,1,0)
v[0] = v1
m = SX(101,1,0)
m[0] = m1

for k in range(100):
  s[k+1] = s[k] + dt*v[k]
  v[k+1] = v[k] + dt / m[k] * (a[k] - alpha * v[k]**2)
  m[k+1] = m[k] - dt * beta * a[k]**2

# Equality constraint function
g = [s[100]-s100, v[100]-v100]

# Create NLP functions
ffcn = SXFunction([u],[f])
ffcn.init()
gfcn = SXFunction([u],[vertcat(g)])
gfcn.init()

# Solve with IPOPT
solver = IpoptSolver(ffcn,gfcn)
solver.setOption("hessian_approximation", "limited-memory")

#solver.setOption("verbose",True)
solver.init()

# Set bounds and initial guess
solver.setInput(umin,     "lbx")
solver.setInput(umax,     "ubx")
solver.setInput(u0,       "x0")
solver.setInput(zeros(2), "lbg")
solver.setInput(zeros(2), "ubg")

# Solve the problem
solver.solve()

# Retrieve the solution
u_opt = array(solver.getOutput("x"))

# Get values at the beginning of each finite element
tgrid = linspace(0,10,10)

# Initial guess
x = u0

# Create SQP method
sqp_solver = SQPMethod(ffcn,gfcn)

# qpOASES
sqp_solver.setOption("qp_solver",QPOasesSolver)
sqp_solver.setOption("qp_solver_options",{"printLevel" : "low"})

# OOQP
#sqp_solver.setOption("qp_solver",OOQpSolver)

sqp_solver.init()
sqp_solver.setInput(umin,     "lbx")
sqp_solver.setInput(umax,     "ubx")
sqp_solver.setInput(u0,       "x0")
sqp_solver.setInput(zeros(2), "lbg")
sqp_solver.setInput(zeros(2), "ubg")
sqp_solver.evaluate()

# Retrieve the solution
u_opt2 = array(sqp_solver.getOutput("x"))

# Plot the results
plt.figure(1)
plt.clf()
#plt.plot(tgrid,x_opt,'--')
#plt.plot(tgrid,y_opt,'-')
plt.plot(tgrid,u_opt,'-.')
plt.plot(tgrid,u_opt2,'-')
plt.title("Rocket trajectory optimization")
plt.xlabel('time')
plt.legend(['IPOPT','SQPMethod'])
plt.grid()







# Parameters in the algorithm
max_iter = 100 # maximum number of sqp iterations
toldx = 1e-12 # stopping criterion for the stepsize
tolgL = 1e-12 # stopping criterion for the lagrangian gradient
merit_mu = 0.  # current 'mu' in the T1 merit function

# Get dimensions
m = gfcn.output().size() # Number of equality constraints
n = len(x)  # Number of variables
q = 0 # Number of inequality constraints

# Initial guess for the lagrange multipliers
lambda_k = zeros(m)
lambda_x_k = zeros(n)

# Initial guess for the Hessian
Bk = eye(n)

# Jacobian function
jfcn = gfcn.jacobian()
jfcn.init()

# Allocate a QP solver
H_sparsity = sp_dense(n,n)
G_sparsity = sp_dense(n)
A_sparsity = jfcn.output().sparsity()

# qpOASES
qp_solver = QPOasesSolver(H_sparsity,A_sparsity)
qp_solver.setOption("printLevel","low")

# IPOPT
#qp_solver =  NLPQpSolver(H_sparsity,A_sparsity)
#qp_solver.setOption("nlp_solver", IpoptSolver)

# OOQP
#qp_solver = OOQpSolver(H_sparsity,A_sparsity)

qp_solver.init()

# No bounds on the control
qp_solver.setInput(-inf,"lbx")
qp_solver.setInput( inf,"ubx")

# Header
print ' k  nls | dx         gradL      eq viol    ineq viol'
k = 0

while True:
  # Evaluate the constraint function
  gfcn.setInput(x)
  gfcn.evaluate()
  gk = gfcn.getOutput()
  
  # Evaluate the Jacobian
  jfcn.setInput(x)
  jfcn.evaluate()
  Jgk = jfcn.getOutput()
  
  # Evaluate the gradient of the objective function
  ffcn.setInput(x)
  ffcn.setAdjSeed(1.0)
  ffcn.evaluate(0,1)
  fk = ffcn.getOutput()
  gfk = DMatrix(ffcn.getAdjSens())
  
  # Pass data to QP solver
  qp_solver.setInput(Bk,"h")
  qp_solver.setInput(Jgk,"a")
  qp_solver.setInput(gfk,"g")
  qp_solver.setInput(-gk,"lba")
  qp_solver.setInput(-gk,"uba")

  # Solve the QP subproblem
  qp_solver.evaluate()

  # Get the optimal solution
  p = qp_solver.getOutput("primal")
  
  # Get the dual solution for the inequalities
  lambda_hat = -qp_solver.getOutput("lambda_a")
  
  # Get the dual solution for the bounds
  lambda_x_hat = -qp_solver.getOutput("lambda_x")
  
  # Get the gradient of the Lagrangian
  gradL = ffcn.getAdjSens() - dot(trans(Jgk),lambda_hat) - lambda_x_hat
  
  ## Pass adjoint seeds to g
  #gfcn.setAdjSeed(lambda_hat)
  #gfcn.evaluate(0,1)

  # Do a line search along p
  #[tk,merit_mu,nlinsearch] = linesearch(ffun,gfun,hfun,x,p,fk,gk,hk,Jfk,Jgk,Jhk,Bk,merit_mu)
  mu = merit_mu
  
  # parameters in the algorithm
  sigma = 1.  # Bk in BDGS is always pos.def.
  rho = 0.5
  mu_safety = 1.1 # safety factor for mu (see below)
  eta = 0.0001 # text to Noc 3.4
  tau = 0.2
  max_iter = 100

  # 1-norm of the feasability violations
  feasviol = sumRows(fabs(gk))

  # Use a quadratic model of T1 to get a lower bound on mu (eq. 18.36 in Nocedal)
  mu_lb = (inner_prod(gfk,p) + sigma/2.0*dot(trans(p),dot(Bk,p)))/(1.-rho)/feasviol

  # Increase mu if it is below the lower bound
  if float(mu) < float(mu_lb):
    mu = mu_lb*mu_safety

  # Calculate T1 at x (18.27 in Nocedal)
  T1 = fk + mu*feasviol

  # Calculate the directional derivative of T1 at x (cf. 18.29 in Nocedal)
  DT1 = inner_prod(gfk,p) - mu*sumRows(fabs(gk))
  
  lsiter = 0
  alpha = 1
  while True:
    # Evaluate prospective x
    x_new = x+alpha*p
    ffcn.setInput(x_new)
    ffcn.evaluate()
    fk_new = ffcn.getOutput()

    # Evaluate gk, hk and get 1-norm of the feasability violations
    gfcn.setInput(x_new)
    gfcn.evaluate()
    gk_new = gfcn.getOutput()
    feasviol_new = sumRows(fabs(gk_new))

    # New T1 function
    T1_new = fk_new + mu*feasviol_new

    # Check Armijo condition, SQP version (18.28 in Nocedal)
    if float(T1_new) <= float(T1 + eta*alpha*DT1):
      break;

    # Backtrace
    alpha = alpha*tau
    
    # Go to next iteration
    lsiter = lsiter+1
    if lsiter >= max_iter:
      raise Exception("linesearch failed!")

  # Step size
  tk = alpha

  # Calculate the new step
  dx = p*tk
  x = x + dx
  lambda_k = tk*lambda_hat + (1-tk)*lambda_k
  lambda_x_k = tk*lambda_x_hat + (1-tk)*lambda_x_k
  k = k+1

  # Gather and print iteration information
  normdx = norm_2(dx) # step size
  normgradL = norm_2(gradL) # size of the Lagrangian gradient
  eq_viol = sumRows(fabs(gk)) # constraint violation
  ineq_viol = nan # sumRows(max(0,-hk)); % inequality constraint violation

  print "%3d %3d |%0.4e %0.4e %0.4e %0.4e" % (k,lsiter,normdx,normgradL,eq_viol,ineq_viol)

  # Check convergence on dx
  if float(normdx) < float(toldx):
    print "Convergence (small dx)"
    break
  elif float(normgradL) < float(tolgL):
    print "Convergence (small gradL)"
    break
    
  # Evaluate the constraint function
  gfcn.setInput(x)
  gfcn.evaluate()
  gk = gfcn.getOutput()
  
  # Evaluate the Jacobian
  jfcn.setInput(x)
  jfcn.evaluate()
  Jgk = jfcn.getOutput()
    
  # Evaluate the gradient of the objective function
  ffcn.setInput(x)
  ffcn.setAdjSeed(1.0)
  ffcn.evaluate(0,1)
  fk = ffcn.getOutput()
  gfk = ffcn.getAdjSens()

  # Check if maximum number of iterations reached
  if k >= max_iter:
    print "Maximum number of SQP iterations reached!"
    break

  # Complete the damped BFGS update (Procedure 18.2 in Nocedal)
  gradL_new = gfk - dot(trans(Jgk),lambda_k) - lambda_x_k
  yk = gradL_new - gradL
  Bdx = dot(Bk,dx)
  dxBdx = dot(trans(dx),Bdx)
  ydx = inner_prod(dx,yk)
  if float(ydx) >= 0.2*float(dxBdx):
    thetak = 1.
  else:
    thetak = 0.8*dxBdx/(dxBdx - ydx)
  rk = thetak*dx + (1-thetak)*Bdx # rk replaces yk to assure Bk pos.def.
  Bk = Bk - outer_prod(Bdx,Bdx)/dxBdx + outer_prod(rk,rk)/ inner_prod(rk,dx)
print "SQP algorithm terminated after %d iterations" % (k-1)

plt.show()










