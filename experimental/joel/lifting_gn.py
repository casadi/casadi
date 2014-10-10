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
from casadi import tools
from matplotlib import pylab as plt
import numpy as NP

# Example 5.3 in Albersmeyer paper

# Automatic initialization
manual_init = True

# Use the Gauss-Newton method
gauss_newton = False

# QP-solver
if False:
  QpSolverClass = OOQpSolver
  qp_solver_options = {}
else:
  QpSolverClass = QPOasesSolver
  qp_solver_options = {"printLevel" : "none"}

# Initial condition
x0_test = [0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.20, 0.30]

#for (i,x0) in enumerate(x0_test[:]):
#for (i,x0) in enumerate([0.02, 0.03, 0.04, 0.05, 0.06]):
#for (i,x0) in enumerate([0.07, 0.08, 0.09, 0.10, 0.20, 0.30]):
for (i,x0) in enumerate([0.08]):

  plt.figure(i+1)
  plt.clf()
  leg = []

  # Bounds on the state
  x_min = -1
  x_max =  1
  xf_min = 0
  xf_max = 0

  # Control
  u_min = -1
  u_max =  1

  # End time
  T = 3.

  # Control discretization
  nk = 30

  # Time step
  dT = T/nk

  # Discretized control
  u = ssym("u",nk)

  # Initial guess for u
  u_guess = DMatrix.zeros(nk)
  u_min = u_min*DMatrix.ones(nk)
  u_max = u_max*DMatrix.ones(nk)

  # Lifted variables
  L = SX()

  # Objective terms
  F = SX()

  # Get an expression for the state that the final time
  x = SX(x0)

  for k in range(nk):
    # Get new value for X
    x = x + dT*(x*(x+1)+u[k])
    
    # Append terms to objective function
    F.append(u[k])
    F.append(x)
    
    # Lift x
    L.append(x)

  # Bounds on G
  G = x
  g_min = xf_min
  g_max = xf_max

  if gauss_newton:
    # Objective function (GN)
    F1 = SXFunction([u],[F])
    
  else:
    # Objective function (SQP)
    F1 = SXFunction([u],[inner_prod(F,F)])

  # Constraint function
  F2 = SXFunction([u],[G])

  # Solve with ipopt
  #nlp_solver = SQPMethod(F1,F2)
  #nlp_solver.setOption("qp_solver",QPOasesSolver)
  #nlp_solver.setOption("qp_solver_options",{"printLevel":"none"})
  nlp_solver = IpoptSolver(F1,F2)
  nlp_solver.init()
  nlp_solver.setInput(u_guess,"x0")
  nlp_solver.setInput(u_min,"lbx")
  nlp_solver.setInput(u_max,"ubx")
  nlp_solver.setInput(g_min,"lbg")
  nlp_solver.setInput(g_max,"ubg")
  nlp_solver.solve()

  #raise Exception("a")

  # Lifting function
  ifcn = SXFunction([u],[L])

  # Problem formulation ends
  # Everything below should go into a lifted newton SQP solver class

  # Options
  tol = 1e-6     # Stopping tolerance
  max_iter = 30  # Maximum number of iterations

  # Extract the free variable and expressions for F and xdef
  u = F1.inputExpr(0)
  f1 = F1.outputExpr(0)
  f2 = F2.outputExpr(0)
  xdef = ifcn.outputExpr(0)

  ## Lifted variables
  x = ssym("x",xdef.size())

  # Substitute in the lifted variables x into the expressions for xdef, F1 and F2
  ex = SXVector([f1,f2])
  substituteInPlace(x, xdef, ex, True)
  [f1,f2] = ex
  
  if gauss_newton: # if Gauss-Newton no multipliers needed
    mux = SX()
    mug = SX()
    
  else: # If SQP, get the gradient of the lagrangian now
    
    # Derivatives of lifted variables
    xdot = ssym("xdot",x.size())
    
    # Lagrange multipliers
    mux = ssym("mux",u.size1())
    mug = ssym("mug",f2.size1())

    # Gradient of the Lagrangian
    xu = vertcat((u,x))
    lgrad = gradient(f1 - inner_prod(mug,f2) + inner_prod(xdot,xdef),xu)

    # Gradient of the Lagrangian
    f1 = lgrad[:u.size1(),0] # + mux # What about the mux term?

    # Definition of xdot
    xdotdef = lgrad[u.size1():,0]
    
    # Reverse direction of x
    xdot[:,0] = SX(list(reversed(list(xdot))))
    xdotdef[:,0] = SX(list(reversed(list(xdotdef))))
    
    # Append to xdef and x
    x.append(xdot)
    xdef.append(xdotdef)
    
  # Residual function G
  G = SXFunction([u,x,mux,mug],[xdef-x,f1,f2])
  G.init()

  # Difference vector d
  d = ssym("d",xdef.size1())

  # Substitute out the x from the zdef
  z = xdef-d
  ex = SXVector([f1,f2])
  substituteInPlace(x, z, ex, False)
  [f1,f2] = ex

  # Modified function Z
  Z = SXFunction([u,d,mux,mug],[z,f1,f2])
  Z.init()

  # Matrix A and B in lifted Newton
  A = Z.jac(0,0)
  B1 = Z.jac(0,1)
  B2 = Z.jac(0,2)
  AB  = SXFunction([u,d,mux,mug],[A,B1,B2])
  AB.init()

  # Variables
  u_k = u_guess
  x_k = DMatrix.zeros(x.shape)
  d_k = DMatrix.zeros(x.shape)
  mux_k = DMatrix.zeros(mux.shape)
  mug_k = DMatrix.zeros(mug.shape)
  dmux_k = DMatrix.zeros(mux.shape)
  dmug_k = DMatrix.zeros(mug.shape)
  f1_k = DMatrix.nan(f1.shape)
  f2_k = DMatrix.nan(f2.shape)

  if manual_init:
    # Initialize node values manually
    G.setInput(u_k,0)
    G.setInput(x_k,1)
    G.setInput(mux_k,2)
    G.setInput(mug_k,3)
    G.evaluate()
    G.getOutput(d_k,0)
    G.getOutput(f1_k,1) # mux is zero (initial multiplier guess)
    G.getOutput(f2_k,2)
  else:
    # Initialize x0 by function evaluation
    Z.setInput(u_k,0)
    Z.setInput(d_k,1)
    Z.setInput(mux_k,2)
    Z.setInput(mug_k,3)
    Z.evaluate()
    Z.getOutput(x_k,0)
    Z.getOutput(f1_k,1) # mux is zero (initial multiplier guess)
    Z.getOutput(f2_k,2)
    
  # Zero seeds
  u0seed = DMatrix.zeros(u.shape)
  d0seed = DMatrix.zeros(d.shape)
  mux0seed = DMatrix.zeros(mux.shape)
  mug0seed = DMatrix.zeros(mug.shape)

  # Iterate
  k = 0
  while True:
    
    # Get A_k and Bk
    AB.setInput(u_k,0)
    AB.setInput(d_k,1)
    AB.setInput(mux_k,2)
    AB.setInput(mug_k,3)
    AB.evaluate()
    A_k = AB.getOutput(0)
    B1_k = AB.getOutput(1) # NOTE: # mux dissappears (constant term)
    B2_k = AB.getOutput(2)
    
    # Get a_k and b_k
    Z.setInput(u_k,0)
    Z.setInput(d_k,1)
    Z.setInput(mux_k,2)
    Z.setInput(mug_k,3)
    Z.setFwdSeed(u0seed,0)
    Z.setFwdSeed(d_k,1)
    Z.setFwdSeed(mux0seed,2)
    Z.setFwdSeed(mug0seed,3)

    Z.evaluate(1,0)
    #Z.getOutput(x_k,0)
    Z.getOutput(f1_k,1)
    Z.getOutput(f2_k,2)
    a_k = -Z.getFwdSens(0)
    b1_k = f1_k-Z.getFwdSens(1) # mux disappears from Z (constant term)
    b2_k = f2_k-Z.getFwdSens(2)

    if gauss_newton:
      # Gauss-Newton Hessian
      H = mul(trans(B1_k),B1_k)
      g = mul(trans(B1_k),b1_k)
      A = B2_k
      a = b2_k
    else:
      # Exact Hessian
      H = B1_k
      g = b1_k # +/- mux_k here?
      A = B2_k
      a = b2_k

    if k==0:
      # Allocate a QP solver
      qp_solver = QpSolverClass(H.sparsity(),A.sparsity())
      qp_solver.setOption(qp_solver_options)
      qp_solver.init()

    # Formulate the QP
    qp_solver.setInput(H,"h")
    qp_solver.setInput(g,"g")
    qp_solver.setInput(A,"a")
    qp_solver.setInput(u_min-u_k,"lbx")
    qp_solver.setInput(u_max-u_k,"ubx")
    qp_solver.setInput(g_min-a,"lba")
    qp_solver.setInput(g_max-a,"uba")

    # Solve the QP
    qp_solver.evaluate()

    # Get the primal solution
    du_k = qp_solver.getOutput("primal")
    
    # Get the dual solution
    if not gauss_newton:
      qp_solver.getOutput(dmux_k,"lambda_x")
      qp_solver.getOutput(dmug_k,"lambda_a")
      dmux_k = -dmux_k
      dmug_k = -dmug_k
    
    # Calculate the step in x
    Z.setFwdSeed(du_k,0)
    Z.setFwdSeed(d0seed,1) # could the a_k term be moved here?
    Z.setFwdSeed(dmux_k,2)
    Z.setFwdSeed(dmug_k,3)
    Z.evaluate(1,0)
    dx_k = Z.getFwdSens(0)
        
    # Take a full step
    u_k =     u_k +   du_k
    x_k =     x_k +    a_k + dx_k
    mug_k = mug_k + dmug_k
    mux_k = mux_k + dmux_k

    # Call algorithm 2 to obtain new d_k and fk
    G.setInput(u_k,0)
    G.setInput(x_k,1)
    G.setInput(mux_k,2)
    G.setInput(mug_k,3)
    G.evaluate()
    G.getOutput(d_k,0)
    G.getOutput(f1_k,1) # mux?
    G.getOutput(f2_k,2)

    # Norm of residual error
    norm_res = sqrt(sum(i*i for i in d_k))

    # Norm of step size
    step_du_k = sum(i*i for i in du_k)
    step_dmug_k = sum(i*i for i in dmug_k)
    norm_step = sqrt(step_du_k + step_dmug_k) # add mux

    # Norm of constraint violation
    viol_umax = sum(max(i,0)**2 for i in u_k-u_max)
    viol_umin = sum(max(i,0)**2 for i in u_min-u_k)
    viol_gmax = sum(max(i,0)**2 for i in f2_k-g_max)
    viol_gmin = sum(max(i,0)**2 for i in g_min-f2_k)
    norm_viol = sqrt(viol_umax + viol_umin + viol_gmax + viol_gmin)

    # Print progress (including the header every 10 rows)
    if k % 10 == 0:
      print " %4s" % "iter", " %20s" % "norm_res", " %20s" % "norm_step", " %20s" % "norm_viol"
    print   " %4d" %  k,     " %20e" %  norm_res,  " %20e" %  norm_step,  " %20e" %  norm_viol
    
    # Check if stopping criteria is satisfied
    if norm_viol + norm_res  + norm_step < tol:
      print "Convergens achieved!"
      break
    
    # Increase iteration count
    k = k+1
    
    # Check if number of iterations have been reached
    if k >= max_iter:
      print "Maximum number of iterations (", max_iter, ") reached"
      break

    # Plot the progress
    plotx = vertcat([x0,x_k[:nk]])
    plott = NP.linspace(0,1,plotx.size1())
    plt.plot(plott,plotx,'*-')
    leg.append(str(k))

  plt.axis([0,1,-x0,2*x0])
  plt.legend(leg)
  plt.grid(True)

  
print "u_k = ", u_k
<<<<<<< HEAD
print "nlp_solver.getOutput("x_opt") = ", nlp_solver.getOutput("x_opt")
=======
print "nlp_solver.getOutput(NLP_SOLVER_X) = ", nlp_solver.getOutput(NLP_SOLVER_X)
>>>>>>> issue_566

plt.show()


