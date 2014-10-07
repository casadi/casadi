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
  t = SX.sym("t")
  
  # Differential equation
  z = SX.sym("z")
  F = SXFunction([z],[z*z - 2*z + 1])
  F.setOption("name","dz/dt")
  F.init()
  
  z0 = -3
  
  # Analytic solution
  z_analytic = SXFunction([t], [(4*t-3)/(3*t+1)])
  z_analytic.setOption("name","analytic solution")
  
  # Collocation point
  tau = SX.sym("tau")

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

    f = SXFunction([tau],[L])
    f.setOption("name", "l(" + str(j) + ")")
    
    # initialize
    f.init()
    l.append(f)
  
  # Get the coefficients of the continuity equation
  D = DMatrix.zeros(K+1)
  for j in range(K+1):
    l[j].setInput(1.)
    l[j].evaluate()
    D[j] = l[j].getOutput()
  print "D = ", D

  # Get the coefficients of the collocation equation using AD
  C = DMatrix.zeros(K+1,K+1)
  for j in range(K+1):
    tfcn = l[j].tangent()
    tfcn.init()
    for k in range(K+1):
      tfcn.setInput(tau_root[k])
      tfcn.evaluate()
      C[j,k] = tfcn.getOutput()
  print "C = ", C
  
  # Collocated states
  Z = SX.sym("Z",N,K+1)
    
  # Construct the NLP
  x = vec(Z.T)
  g = SX()
  for i in range(N):
    for k in range(1,K+1):
      # Add collocation equations to NLP
      rhs = 0
      for j in range(K+1):
        rhs += Z[i,j]*C[j,k]
      [FF] = F([Z[i,k]])
      g.append(h*FF-rhs)

    # Add continuity equation to NLP
    rhs = 0
    for j in range(K+1):
      rhs += D[j]*Z[i,j]

    if(i<N-1):
      g.append(Z[i+1,0] - rhs)

  print "g = ", g

  # NLP
  nlp = SXFunction(nlpIn(x=x),nlpOut(f=x[0]**2,g=g))

  ## ----
  ## SOLVE THE NLP
  ## ----
  
  # Allocate an NLP solver
  solver = NlpSolver("ipopt", nlp)

  # Set options
  solver.setOption("tol",1e-10)

  # initialize the solver
  solver.init()

  # Initial condition
  xinit = x.size() * [0]
  solver.setInput(xinit,"x0")

  # Bounds on x
  lbx = x.size()*[-100]
  ubx = x.size()*[100]
  lbx[0] = ubx[0] = z0
  solver.setInput(lbx,"lbx")
  solver.setInput(ubx,"ubx")
  
  # Bounds on the constraints
  lubg = g.size()*[0]
  solver.setInput(lubg,"lbg")
  solver.setInput(lubg,"ubg")
  
  # Solve the problem
  solver.evaluate()
  
  ## Print the time points
  t_opt = N*(K+1) * [0]
  for i in range(N):
    for j in range(K+1):
      t_opt[j + (K+1)*i] = h*(i + tau_root[j])
  
  print "time points: ", t_opt

  # Print the optimal cost
  print "optimal cost: ", float(solver.getOutput("f"))

  # Print the optimal solution
  xopt = solver.output("x").data()
  print "optimal solution: ", xopt
 
  # plot to screen
  plt.plot(t_opt,xopt)

# show the plots
plt.show()
  
