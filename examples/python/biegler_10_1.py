# -*- coding: utf-8 -*-
from casadi import *
from numpy import *
import matplotlib.pyplot as plt

# Excercise 1, chapter 10 from Larry Biegler's book
print "program started"

# Test with different number of elements
for N in range(1,11):
  print "N = ", N
  
  # Degree of interpolating polynomial
  K = 2
  
  # Legrandre roots
  tau_root = [0., 0.211325, 0.788675]

  # Radau roots (K=3)
  #tau_root = [0, 0.155051, 0.644949, 1]

  # Time
  t = SX("t")
  
  # Differential equation
  z = SX("z")
  F = SXFunction([[z]],[[z*z - 2*z + 1]]);
  F.setOption("name","dz/dt")
  
  z0 = -3
  
  # Analytic solution
  z_analytic = SXFunction([[t]], [[(4*t-3)/(3*t+1)]]);
  z_analytic.setOption("name","analytic solution");
  
  # Collocation point
  tau = SX("tau")

  # Step size
  h = 1.0/N
  
  # Lagrange polynomials
  l = []
  for j in range(K+1):
    L = 1
    for k in range(K+1):
      if(k != j):
        L *= (tau-tau_root[k])/(tau_root[j]-tau_root[k])

    print "l(", j, ") = ", L

    f = SXFunction([[tau]],[[L]])
    f.setOption("ad_order",1)
    f.setOption("name", "l(" + str(j) + ")")
    
    # initialize
    f.init()
    l.append(f)
  
  # Get the coefficients of the continuity equation
  D = []
  for j in range(K+1):
    l[j].setInput(1.)
    l[j].evaluate()
    res = l[j].getOutput()
    D.append(res[0])

  print "D = ", D

  # Get the coefficients of the collocation equation using AD
  C = []
  for j in range(K+1):
    Cj = []
    for k in range(K+1):
      l[j].setInput(tau_root[k])
      l[j].setFwdSeed(1.0)
      l[j].evaluate(1,0)
      sens = l[j].getFwdSens()
      Cj.append(sens[0])
    C.append(Cj)
  
  print "C = ", C
  
  # Collocated states
  Z = []
  for i in range(N):
    Zi = []
    for j in range(K+1):
      Zi.append( SX("Z_"+str(i)+"_"+str(j)))
    Z.append(Zi)
  print "Z = ", Z
    
  # State at final time
  ZF = SX("ZF")
  
  # All variables
  x = []
  for i in range(N):
    for j in range(K+1):
      x.append(Z[i][j])
  print "x = ", x
  
  # Construct the "NLP"
  g = []
  for i in range(N):
    for k in range(1,K+1):
      # Add collocation equations to NLP
      rhs = 0
      for j in range(K+1):
        rhs += Z[i][j]*C[j][k]
      FF = F.eval([[Z[i][k]]])
      g.append(h*FF[0][0]- rhs)

    # Add continuity equation to NLP
    rhs = 0
    for j in range(K+1):
      rhs += D[j]*Z[i][j]

    if(i<N-1):
      g.append(Z[i+1][0] - rhs)

  print "g = ", g
  
  # Constraint function
  gfcn = SXFunction([x],[g])

  # Dummy objective function
  obj = SXFunction([x], [[x[0]*x[0]]])
  
  ## ----
  ## SOLVE THE NLP
  ## ----
  
  # Allocate an NLP solver
  solver = IpoptSolver(obj,gfcn)

  # Set options
  solver.setOption("tol",1e-10);
  solver.setOption("hessian_approximation","limited-memory");

  # initialize the solver
  solver.init();

  # Initial condition
  xinit = len(x) * [0]
  solver.setInput(xinit,NLP_X_INIT)

  # Bounds on x
  lbx = len(x)*[-100]
  ubx = len(x)*[100]
  lbx[0] = ubx[0] = z0
  solver.setInput(lbx,NLP_LBX)
  solver.setInput(ubx,NLP_UBX)
  
  # Bounds on the constraints
  lubg = len(g)*[0]
  solver.setInput(lubg,NLP_LBG)
  solver.setInput(lubg,NLP_UBG)
  
  # Solve the problem
  solver.solve()
  
  ## Print the time points
  t_opt = N*(K+1) * [0]
  for i in range(N):
    for j in range(K+1):
      t_opt[j + (K+1)*i] = h*(i + tau_root[j])
  
  print "time points: ", t_opt

  # Print the optimal cost
  print "optimal cost: ", tuple(solver.getOutput(NLP_COST))

  # Print the optimal solution
  xopt = tuple(solver.getOutput(NLP_X_OPT))
  print "optimal solution: ", xopt
 
  # plot to screen
  plt.plot(t_opt,xopt)

# show the plots
plt.show()
  