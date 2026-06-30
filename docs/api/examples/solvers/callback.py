#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             KU Leuven. All rights reserved.
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
# Callback
# =====================

import casadi as ca
from numpy import *

# In this example, we will demonstrate callback functionality for Ipopt.
# Note that you need the fix https://github.com/casadi/casadi/wiki/enableIpoptCallback before this works
#
# We start with constructing the rosenbrock problem

x=ca.SX.sym("x")
y=ca.SX.sym("y")

f = (1-x)**2+100*(y-x**2)**2
nlp={'x':ca.vertcat(x,y), 'f':f,'g':x+y}
fcn = ca.Function('f', [x, y], [f])

import matplotlib
if "Agg" not in matplotlib.get_backend():
  matplotlib.interactive(True)

from pylab import figure, subplot, contourf, colorbar, draw, show, plot, title

import time

class MyCallback(ca.Callback):
  def __init__(self, name, nx, ng, np, opts={}):
    ca.Callback.__init__(self)

    self.nx = nx
    self.ng = ng
    self.np = np

    figure(1)

    x_,y_ = mgrid[-1:1.5:0.01,-1:1.5:0.01]
    z_ = ca.DM.zeros(x_.shape)

    for i in range(x_.shape[0]):
      for j in range(x_.shape[1]):
        z_[i,j] = fcn(x_[i,j],y_[i,j])
    contourf(x_,y_,z_)
    colorbar()
    title('Iterations of Rosenbrock')
    draw()

    self.x_sols = []
    self.y_sols = []

    # Initialize internal objects
    self.construct(name, opts)

  def get_n_in(self): return ca.nlpsol_n_out()
  def get_n_out(self): return 1
  def get_name_in(self, i): return ca.nlpsol_out(i)
  def get_name_out(self, i): return "ret"

  def get_sparsity_in(self, i):
    n = ca.nlpsol_out(i)
    if n=='f':
      return ca.Sparsity. scalar()
    elif n in ('x', 'lam_x'):
      return ca.Sparsity.dense(self.nx)
    elif n in ('g', 'lam_g'):
      return ca.Sparsity.dense(self.ng)
    else:
      return ca.Sparsity(0,0)
  def eval(self, arg):
    # Create dictionary
    darg = {}
    for (i,s) in enumerate(ca.nlpsol_out()): darg[s] = arg[i]

    sol = darg['x']
    self.x_sols.append(float(sol[0]))
    self.y_sols.append(float(sol[1]))

    if hasattr(self,'lines'):
      if "template" not in matplotlib.get_backend(): # Broken for template: https://github.com/matplotlib/matplotlib/issues/8516/
        self.lines[0].set_data(self.x_sols,self.y_sols)

    else:
      self.lines = plot(self.x_sols,self.y_sols,'or-')

    draw()
    time.sleep(0.25)

    return [0]

mycallback = MyCallback('mycallback', 2, 1, 0)
opts = {}
opts['iteration_callback'] = mycallback
opts['ipopt.tol'] = 1e-8
opts['ipopt.max_iter'] = 50
solver = ca.nlpsol('solver', 'ipopt', nlp, opts)
sol = solver(lbx=-10, ubx=10, lbg=-10, ubg=10)

# By setting matplotlib interactivity off, we can inspect the figure at ease

matplotlib.interactive(False)
show()
