#
#     MIT No Attribution
#
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
#
#     Permission is hereby granted, free of charge, to any person obtaining a copy of this
#     software and associated documentation files (the "Software"), to deal in the Software
#     without restriction, including without limitation the rights to use, copy, modify,
#     merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
#     permit persons to whom the Software is furnished to do so.
#
#     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#     INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#     PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
#
# -*- coding: utf-8 -*-
from casadi import *
from numpy import *
import matplotlib.pyplot as plt

# Excercise 1, chapter 10 from Larry Biegler's book
print("program started")

# Test with different number of elements
for N in range(1,11):
  print("N = ", N)

  # Degree of interpolating polynomial
  K = 2

  # Legrandre roots
  tau_root = [0., 0.211325, 0.788675]

  # Radau roots (K=3)
  #tau_root = [0, 0.155051, 0.644949, 1]

  # Time
  t = SX.sym("t")

  # Differential equation
  z = SX.sym("z")
  F = Function("dz_dt", [z],[z*z - 2*z + 1])

  z0 = -3

  # Analytic solution
  z_analytic = Function("z_analytic", [t], [(4*t-3)/(3*t+1)])

  # Collocation point
  tau = SX.sym("tau")

  # Step size
  h = 1.0/N

  # Get the coefficients of the continuity and collocation equations
  D = DM.zeros(K+1)
  C = DM.zeros(K+1,K+1)
  for j in range(K+1):
    # Lagrange polynomial
    L = 1
    for k in range(K+1):
      if(k != j):
        L *= (tau-tau_root[k])/(tau_root[j]-tau_root[k])

    # Evaluate at end for coefficients of continuity equation
    lfcn = Function("lfcn", [tau],[L])
    D[j] = lfcn(1.)

    # Differentiate and evaluate at collocation points
    tfcn = Function("tfcn", [tau],[tangent(L,tau)])
    for k in range(K+1): C[j,k] = tfcn(tau_root[k])
  print("C = ", C)
  print("D = ", D)

  # Collocated states
  Z = SX.sym("Z",N,K+1)

  # Construct the NLP
  x = vec(Z.T)
  g = []
  for i in range(N):
    for k in range(1,K+1):
      # Add collocation equations to NLP
      rhs = 0
      for j in range(K+1):
        rhs += Z[i,j]*C[j,k]
      FF = F(Z[i,k])
      g.append(h*FF-rhs)

    # Add continuity equation to NLP
    rhs = 0
    for j in range(K+1):
      rhs += D[j]*Z[i,j]

    if(i<N-1):
      g.append(Z[i+1,0] - rhs)

  g = vertcat(*g)

  print("g = ", g)

  # NLP
  nlp = {'x':x, 'f':x[0]**2, 'g':g}

  ## ----
  ## SOLVE THE NLP
  ## ----

  # NLP solver options
  opts = {"ipopt.tol" : 1e-10}

  # Allocate an NLP solver and buffer
  solver = nlpsol("solver", "ipopt", nlp, opts)
  arg = {}

  # Initial condition
  arg["x0"] = x.nnz() * [0]

  # Bounds on x
  lbx = x.nnz()*[-100]
  ubx = x.nnz()*[100]
  lbx[0] = ubx[0] = z0
  arg["lbx"] = lbx
  arg["ubx"] = ubx

  # Bounds on the constraints
  arg["lbg"] = 0
  arg["ubg"] = 0

  # Solve the problem
  res = solver(**arg)

  ## Print the time points
  t_opt = N*(K+1) * [0]
  for i in range(N):
    for j in range(K+1):
      t_opt[j + (K+1)*i] = h*(i + tau_root[j])

  print("time points: ", t_opt)

  # Print the optimal cost
  print("optimal cost: ", float(res["f"]))

  # Print the optimal solution
  xopt = res["x"].nonzeros()
  print("optimal solution: ", xopt)

  # plot to screen
  plt.plot(t_opt,xopt)

# show the plots
plt.show()
