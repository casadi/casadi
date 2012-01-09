#! Callback
#! =====================
from casadi import *
from numpy import *

#! In this example, we will demonstrate callback functionality for Ipopt.
#! Note that you need the fix https://sourceforge.net/apps/trac/casadi/wiki/enableIpoptCallback  before this works
#!
#! We start with constructing the rosenbrock problem
x=SX("x")
y=SX("y")

f=SXFunction([[x,y]],[(1-x)**2+100*(y-x**2)**2])
g=SXFunction([[x,y]],[x+y])
    
#! Simple callback
#! ===============
#! First we demonstrate a callback that does nothing more than printing some information on the screen

#! We create a callable python class, i.e. oen that has a __call__ member:

class MyCallback:
  def __init__(self):
    self.iter = 0 
  def __call__(self,f,*args):
    print "====Hey, I'm an iteration===="
    print "X_OPT = ", f.input(NLP_X_OPT)
    print f.getStats()
    self.iter = self.iter + 1
    if self.iter > 5:
      print "5 Iterations, that is quite enough!"
      f.output(0).set(1) # With this statement you can halt the iterations

mycallback = MyCallback()

#! We create a casadi function out of this callable object.
#! The sparsities given here as input must match the sparsities of the outputs of our NLP Solver
c = PyFunction( mycallback, [sp_dense(2,1),sp_dense(1,1),sp_dense(1,1),sp_dense(2,1),sp_dense(1,1)], [sp_dense(1,1)] )
c.init()


solver = IpoptSolver(f,g)
solver.setOption("iteration_callback",c)
solver.setOption("tol",1e-8)
solver.setOption("max_iter",20)
solver.init()
solver.input(NLP_LBX).set([-10]*2)
solver.input(NLP_UBX).set([10]*2)
solver.input(NLP_LBG).set([-10])
solver.input(NLP_UBG).set([10])
solver.solve()

#! Matplotlib callback
#! ===================
#! Now let's do some useful visualisations

from pylab import figure, subplot, contourf, colorbar, draw, show, plot, title

import time
import matplotlib
matplotlib.interactive(True)
  
class MyCallback:
  def __init__(self):
    figure(1)
    subplot(111)
    
    x_,y_ = mgrid[-1:1.5:0.01,-1:1.5:0.01]
    z_ = DMatrix.zeros(x_.shape)
    
    for i in range(x_.shape[0]):
      for j in range(x_.shape[1]):
        f.input().set([x_[i,j],y_[i,j]])
        f.evaluate()
        z_[i,j] = float(f.output())
    contourf(x_,y_,z_)
    colorbar()
    title('Iterations of Rosenbrock')
    draw()
    
    self.x_sols = []
    self.y_sols = []
    
  def __call__(self,f,*args):
    sol = f.input(NLP_X_OPT)
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
    
mycallback = MyCallback()

#! We create a casadi function out of this callable object.
#! The sparsities given here as input must match the sparsities of the outputs of our NLP Solver
c = PyFunction( mycallback, [sp_dense(2,1),sp_dense(1,1),sp_dense(1,1),sp_dense(2,1),sp_dense(1,1)], [sp_dense(1,1)] )
c.init()


solver = IpoptSolver(f,g)
solver.setOption("iteration_callback",c)
solver.setOption("tol",1e-8)
solver.setOption("max_iter",50)
solver.init()
solver.input(NLP_LBX).set([-10]*2)
solver.input(NLP_UBX).set([10]*2)
solver.input(NLP_LBG).set([-10])
solver.input(NLP_UBG).set([10])
solver.solve()

#! By setting matplotlib interactivity off, we can inspect the figure at ease
matplotlib.interactive(False)
show()
