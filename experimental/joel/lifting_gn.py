from casadi import *
from casadi import tools
from matplotlib import pylab as plt
import numpy as NP

# Example 5.3 in Albersmeyer paper
# Solve F(u) = u**16 - 2 == 0

# Initial condition
x0_test = [0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.20, 0.30]

# Automatic initialization
manual_init = True

# Use the Gauss-Newton method
gauss_newton = False

#for (i,x0) in enumerate(x0_test[:]):
for (i,x0) in enumerate([0.02]):

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
  L = []

  # Objective terms
  F = []

  # Get an expression for the state that the final time
  x = SXMatrix(x0)

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

  # Lifting function
  ifcn = SXFunction([u],[vertcat(L)])

  # Problem formulation ends
  # Everything below should go into a lifted newton SQP solver class

  # Options
  TOL = 1e-6     # Stopping tolerance
  max_iter = 30  # Maximum number of iterations

  # Extract the free variable and expressions for F and xdef
  u = F1.inputSX()
  f1 = F1.outputSX()
  f2 = F2.outputSX()
  xdef = ifcn.outputSX()

  if gauss_newton: # if Gauss-Newton no multipliers needed
    mux = SXMatrix()
    mug = SXMatrix()
    
  else: # If SQP, get the gradient of the lagrangian now
    
    # Lagrange multipliers
    mux = ssym("mux",u.size())
    mug = ssym("mug",f2.size())
    
    # Lagrange function
    lag = f1 + mul(trans(mux),u) + mul(trans(mug),f2)

    # Gradient of the Lagrangian
    f1 = jacobian(lag,u)

  ## Lifted variables
  x = ssym("x",xdef.size())

  # Substitute in the lifted variables x into the expressions for xdef, F1 and F2
  ex = SXMatrixVector([f1,f2])
  substituteInPlace(x, xdef, ex, True, False)
  [f1,f2] = ex

  # Residual function G
  G = SXFunction([u,x,mux,mug],[xdef-x,f1,f2])
  G.init()

  # Difference vector d
  d = ssym("d",xdef.size())

  # Substitute out the x from the zdef
  z = xdef-d
  ex = SXMatrixVector([f1,f2])
  substituteInPlace(x, z, ex, False, False)
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
  dmux_k = DMatrix.zeros(mux.shape)
  mug_k = DMatrix.zeros(mug.shape)
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
    G.getOutput(f1_k,1)
    G.getOutput(f2_k,2)
  else:
    # Initialize x0 by function evaluation
    Z.setInput(u_k,0)
    Z.setInput(d_k,1)
    Z.setInput(mux_k,2)
    Z.setInput(mug_k,3)
    Z.evaluate()
    Z.getOutput(x_k,0)
    Z.getOutput(f1_k,1)
    Z.getOutput(f2_k,2)
    
  # Print header
  print " %4s" % "iter", " %20s" % "norm_f1_k", " %20s" % "norm_f2_k", " %20s" % "norm_d_k", " %20s" % "norm_du_k", " %20s" % "feas_viol"

  # No seed in the u direction
  useed = DMatrix.zeros(u.shape)

  # Iterate
  k = 0
  while True:
    
    # Get A_k and Bk
    AB.setInput(u_k,0)
    AB.setInput(d_k,1)
    AB.setInput(mux_k,2)
    AB.setInput(mug_k,3)
    AB.evaluate()
    A_k = AB.output(0)
    B1_k = AB.output(1)
    B2_k = AB.output(2)
    
    # Get a_k and b_k
    Z.setInput(u_k,0)
    Z.setInput(d_k,1)
    Z.setInput(mux_k,2)
    Z.setInput(mug_k,3)
    Z.setFwdSeed(useed,0)
    Z.setFwdSeed(d_k,1)
    Z.evaluate(1,0)
    #Z.getOutput(x_k,0)
    Z.getOutput(f1_k,1)
    Z.getOutput(f2_k,2)
    a_k = -Z.fwdSens(0)
    b1_k = f1_k-Z.fwdSens(1)
    b2_k = f2_k-Z.fwdSens(2)

    if gauss_newton:
      # Gauss-Newton Hessian
      H = mul(trans(B1_k),B1_k)
      g = mul(trans(B1_k),b1_k)
      A = B2_k
      a = b2_k
    else:
      # Exact Hessian
      H = B1_k
      g = b1_k-mul(B2_k,mug_k)
      A = B2_k
      a = b2_k
        
    if k==0:
      #H.printDense()
      #g.printDense()
      #A.printDense()
      #a.printDense()
      
      #print "mug_k = ", mug_k
      #print "b1_k = ", b1_k
      #print g_min-a, ",", g_max-a
      
      # Allocate a QP solver
      qp_solver = OOQPSolver(H.sparsity(),A.sparsity())
      qp_solver.init()

    # Formulate the QP
    qp_solver.setInput(H,QP_H)
    qp_solver.setInput(g,QP_G)
    qp_solver.setInput(A,QP_A)
    qp_solver.setInput(u_min,QP_LBX)
    qp_solver.setInput(u_max,QP_UBX)
    qp_solver.setInput(g_min-a,QP_LBA)
    qp_solver.setInput(g_max-a,QP_UBA)

    # Solve the QP
    qp_solver.evaluate()

    # Get the primal solution
    du_k = qp_solver.output(QP_PRIMAL)
    
    print "du_k = ", du_k
    
    # Get the dual solution
    if not gauss_newton:
      dmux_k = qp_solver.output(QP_DUAL_X)
      dmug_k = qp_solver.output(QP_DUAL_A)
      
      print "dmux_k = ", dmux_k
      print "dmug_k = ", dmug_k
    
    # Perform the full Newton step
    x_k = x_k + a_k + mul(A_k,du_k)
    u_k = u_k + du_k
    #mux_k = mux_k + dmux_k
    #mug_k = mug_k + dmug_k
    mux_k = dmux_k
    mug_k = dmug_k
    
    print "x_k = ", x_k
    
    # Call algorithm 2 to obtain new d_k and fk
    G.setInput(u_k,0)
    G.setInput(x_k,1)
    G.setInput(mux_k,2)
    G.setInput(mug_k,3)
    G.evaluate()
    G.getOutput(d_k,0)
    G.getOutput(f1_k,1)
    G.getOutput(f2_k,2)
    
    # Get error
    norm_f1_k = float(norm_2(f1_k))
    norm_f2_k = float(norm_2(f2_k))
    norm_d_k = float(norm_2(d_k))
    norm_du_k = float(norm_2(du_k))
    
    # Constraint violation
    viol_umax = float(norm_2(max(u_k-u_max,0)))
    viol_umin = float(norm_2(max(u_min-u_k,0)))
    viol_xmax = float(norm_2(max(f2_k-g_max,0)))
    viol_xmin = float(norm_2(max(g_min-f2_k,0)))
    viol_u = float(norm_2([viol_umin,viol_umax]))
    viol_x = float(norm_2([viol_xmin,viol_xmax]))
    feas_viol = float(norm_2([viol_u,viol_x]))
    
    # Print
    print " %4d" % k, " %20e" % norm_f1_k, " %20e" % norm_f2_k, " %20e" % norm_d_k, " %20e" % norm_du_k, " %20e" % feas_viol
    
    # Check if stopping criteria achieved
    if feas_viol + norm_d_k  + norm_du_k < TOL:
      print "Convergens achieved!"
      break
    
    # Increase iteration count
    k = k+1
    
    # Check if number of iterations have been reached
    if k >= max_iter:
      print "Maximum number of iterations (", max_iter, ") reached"
      break

    plotx = vertcat([x0,x_k])
    plott = NP.linspace(0,1,plotx.size1())
    plt.plot(plott,plotx,'*-')
    leg.append(str(k))

  plt.axis([0,1,-x0,2*x0])
  plt.legend(leg)
  plt.grid(True)

plt.show()


