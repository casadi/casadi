from casadi import *
from casadi import tools
from matplotlib import pylab as plt
import numpy as NP

# Example 5.3 in Albersmeyer paper
# Solve F(u) = u**16 - 2 == 0

# Initial condition
x0_test = [0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.20, 0.30]

# Automatic initialization
auto_init = False

for (i,x0) in enumerate(x0_test[:]):
#for (i,x0) in enumerate([0.02]):

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

  # Constraint function
  G = []

  # Get an expression for the state that the final time
  x = SXMatrix(x0)

  # Lift the initial conditions
  #L.append(x)

  for k in range(nk):
    # Get new value for X
    x = x + dT*(x*(x+1)+u[k])
    
    # Append terms to objective function
    F.append(u[k])
    F.append(x)
    
    # Append terms to constraint function
    #G.append(x)
    
    # Lift x
    #if k==0: 
    L.append(x)

  # Bounds on G
  #g_min = x_min*DMatrix.ones(nk)
  #g_max = x_max*DMatrix.ones(nk)
  #g_min[-1] = xf_min
  #g_max[-1] = xf_max
  G.append(x)
  g_min = xf_min
  g_max = xf_max

  # Objective function (GN)
  F1 = SXFunction([u],[F])

  # Constraint function
  F2 = SXFunction([u],[vertcat(G)])

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

  ## Lifted variables
  x = ssym("x",xdef.size())

  # Substitute in the lifted variables x into the expressions for xdef, F1 and F2
  ex = SXMatrixVector([f1,f2])
  substituteInPlace(x, xdef, ex, True, False)
  [f1,f2] = ex

  # Residual function G
  G = SXFunction([u,x],[xdef-x,f1,f2])
  G.init()

  # Difference vector d
  d = ssym("d",xdef.size())

  # Substitute out the x from the zdef
  z = xdef-d
  ex = SXMatrixVector([f1,f2])
  substituteInPlace(x, z, ex, False, False)
  [f1,f2] = ex

  # Modified function Z
  Z = SXFunction([u,d],[z,f1,f2])
  Z.init()

  # Matrix A and B in lifted Newton
  A = Z.jac(0,0)
  B1 = Z.jac(0,1)
  B2 = Z.jac(0,2)
  AB  = SXFunction([u,d],[A,B1,B2])
  AB.init()

  # Variables
  uk = u_guess
  dk = DMatrix.zeros(xdef.size())
  xk = DMatrix.zeros(xdef.size())
  f1k = DMatrix.nan(f1.shape)
  f2k = DMatrix.nan(f2.shape)

  if auto_init:
    # Initialize x0 by function evaluation
    Z.setInput(uk,0)
    Z.setInput(dk,1)
    Z.evaluate()
    Z.getOutput(xk,0)
    Z.getOutput(f1k,1)
    Z.getOutput(f2k,2)
  else:
    # Initialize node values manually
    G.setInput(uk,0)
    G.setInput(xk,1)
    G.evaluate()
    G.getOutput(dk,0)
    G.getOutput(f1k,1)
    G.getOutput(f2k,2)
    
  #print float(norm_2(f1k))
  #print float(norm_2(vertcat((max(f2k-g_max,0),max(g_min-f2k,0)))))
  #print float(norm_2(max(f2k-g_max,0)+max(g_min-f2k)))
  #print max(f2k-g_max,0)+max(g_min-f2k,0)
  


  # Print header
  print " %4s" % "iter", " %20s" % "norm_f1k", " %20s" % "norm_f2k", " %20s" % "norm_dk", " %20s" % "norm_du", " %20s" % "feas_viol"

  # No seed in the u direction
  useed = DMatrix.zeros(u.shape)

  # Iterate
  k = 0
  while True:
    
    # Get Ak and Bk
    AB.setInput(uk,0)
    AB.setInput(dk,1)
    AB.evaluate()
    Ak = AB.output(0)
    B1k = AB.output(1)
    B2k = AB.output(2)
    
    # Get ak and bk
    Z.setInput(uk,0)
    Z.setInput(dk,1)
    Z.setFwdSeed(useed,0)
    Z.setFwdSeed(dk,1)
    Z.evaluate(1,0)
    #Z.getOutput(xk,0)
    Z.getOutput(f1k,1)
    Z.getOutput(f2k,2)
    ak = -Z.fwdSens(0)
    b1k = f1k-Z.fwdSens(1)
    b2k = f2k-Z.fwdSens(2)

    # Gauss-Newton Hessian and linear term
    H = mul(trans(B1k),B1k)
    g = mul(trans(B1k),b1k)
    A = B2k
    
    
    if k==0:
      # Allocate a QP solver
      qp_solver = OOQPSolver(H.sparsity(),A.sparsity())
      qp_solver.init()
    
    qp_solver.setInput(H,QP_H)
    qp_solver.setInput(g,QP_G)
    qp_solver.setInput(A,QP_A)
    qp_solver.setInput(u_min,QP_LBX)
    qp_solver.setInput(u_max,QP_UBX)
    qp_solver.setInput(g_min-b2k,QP_LBA)
    qp_solver.setInput(g_max-b2k,QP_UBA)
    qp_solver.evaluate()
    du = qp_solver.output(QP_PRIMAL)
    
    # Perform the Newton step
    xk = xk + ak + mul(Ak,du)
    uk = uk + du
    
    # Call algorithm 2 to obtain new dk and fk
    G.setInput(uk,0)
    G.setInput(xk,1)
    G.evaluate()
    G.getOutput(dk,0)
    G.getOutput(f1k,1)
    G.getOutput(f2k,2)
    
    # Get error
    norm_f1k = float(norm_2(f1k))
    norm_f2k = float(norm_2(f2k))
    norm_dk = float(norm_2(dk))
    norm_du = float(norm_2(du))
    
    # Constraint violation
    viol_umax = float(norm_2(max(uk-u_max,0)))
    viol_umin = float(norm_2(max(u_min-uk,0)))
    viol_xmax = float(norm_2(max(f2k-g_max,0)))
    viol_xmin = float(norm_2(max(g_min-f2k,0)))
    viol_u = float(norm_2([viol_umin,viol_umax]))
    viol_x = float(norm_2([viol_xmin,viol_xmax]))
    feas_viol = float(norm_2([viol_u,viol_x]))
    
    # Print
    print " %4d" % k, " %20e" % norm_f1k, " %20e" % norm_f2k, " %20e" % norm_dk, " %20e" % norm_du, " %20e" % feas_viol
    
    # Check if stopping criteria achieved
    if feas_viol + norm_dk  + norm_du < TOL:
      print "Convergens achieved!"
      break
    
    # Increase iteration count
    k = k+1
    
    # Check if number of iterations have been reached
    if k >= max_iter:
      print "Maximum number of iterations (", max_iter, ") reached"
      break

    plotx = vertcat([x0,xk])
    plott = NP.linspace(0,1,plotx.size1())
    plt.plot(plott,plotx,'*-')
    leg.append(str(k))

  plt.axis([0,1,-x0,2*x0])
  plt.legend(leg)
  plt.grid(True)

plt.show()


