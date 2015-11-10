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
#! Callback
#! =====================
from casadi import *
from numpy import *

#! In this example, we will demonstrate callback functionality for Ipopt.
#! Note that you need the fix https://github.com/casadi/casadi/wiki/enableIpoptCallback before this works
#!
#! We start with constructing the rosenbrock problem
x=SX.sym("x")
y=SX.sym("y")

f = (1-x)**2+100*(y-x**2)**2
nlp={'x':vertcat((x,y)), 'f':f,'g':x+y}
fcn = Function('f', [x, y], [f])
 
# #! Simple callback
# #! ===============
# #! First we demonstrate a callback that does nothing more than printing some information on the screen

# it = 0

# @pycallback
# def simplecallback(f):
#   global it
#   print "====Hey, I'm an iteration===="
#   print "X_OPT = ", f.getOutput("x")
#   print f.getStat("iteration")
#   it = it + 1
#   if it > 5:
#     print "5 Iterations, that is quite enough!"
#     return 1 # With this statement you can halt the iterations
#   else:
#     return 0

# opts = {}
# opts["iteration_callback"] = simplecallback
# opts["tol"] = 1e-8
# opts["max_iter"] = 20
# solver = Function.nlpsol("solver", "ipopt", nlp, opts)
# solver.setInput([-10]*2,"lbx")
# solver.setInput([10]*2,"ubx")
# solver.setInput([-10],"lbg")
# solver.setInput([10],"ubg")
# solver.evaluate()

#! Matplotlib callback
#! ===================
#! Now let's do some useful visualisations
#! We create a callable python class, i.e. one that has a __call__ member

from pylab import figure, subplot, contourf, colorbar, draw, show, plot, title

import time
import matplotlib
matplotlib.interactive(True)

class MyCallback(Callback):
  def __init__(self, name, nx, ng, np, opts={}):
    Callback.__init__(self)

    self.nx = nx
    self.ng = ng
    self.np = np

    opts['input_scheme'] = Function.nlpsol_out()
    opts['output_scheme'] = ['ret']

    figure(1)
    subplot(111)
    
    x_,y_ = mgrid[-1:1.5:0.01,-1:1.5:0.01]
    z_ = DM.zeros(x_.shape)
    
    for i in range(x_.shape[0]):
      for j in range(x_.shape[1]):
        [z_[i,j]] = fcn([x_[i,j],y_[i,j]])
    contourf(x_,y_,z_)
    colorbar()
    title('Iterations of Rosenbrock')
    draw()
    
    self.x_sols = []
    self.y_sols = []

    # Initialize internal objects
    self.construct(name, opts)

  def get_n_in(self): return Function.nlpsol_n_out()
  def get_n_out(self): return 1


  def get_input_shape(self, i):
    n = Function.nlpsol_out(i)
    if n=='f':
      return (1,1)
    elif n in ('x', 'lam_x'):
      return (self.nx, 1)
    elif n in ('g', 'lam_g'):
      return (self.ng, 1)
    else:
      return (0, 0)
  def eval(self, arg):
    # Create dictionary
    darg = {}
    for (i,s) in enumerate(Function.nlpsol_out()): darg[s] = arg[i]

    sol = darg['x']
    self.x_sols.append(float(sol[0]))
    self.y_sols.append(float(sol[1]))
    subplot(111)
    if hasattr(self,'lines'):
      self.lines[0].set_xdata(self.x_sols)
      self.lines[0].set_ydata(self.y_sols)
    else:
      self.lines = plot(self.x_sols,self.y_sols,'or-')

    draw()
    time.sleep(0.25)

    return [0]

mycallback = MyCallback('mycallback', 2, 1, 0)
opts = {}
opts["iteration_callback"] = mycallback
opts["tol"] = 1e-8
opts["max_iter"] = 50
solver = Function.nlpsol("solver", "ipopt", nlp, opts)
solver.setInput([-10]*2,"lbx")
solver.setInput([10]*2,"ubx")
solver.setInput([-10],"lbg")
solver.setInput([10],"ubg")
solver.evaluate()

#! By setting matplotlib interactivity off, we can inspect the figure at ease
matplotlib.interactive(False)
show()
