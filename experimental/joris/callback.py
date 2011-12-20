from casadi import *
import numpy as NP

x = ssym("x")
y = ssym("y")
z = ssym("z")
v = vertcat([x,y,z])

f = SXFunction([v],[x**2 + 100*z**2])
g = SXFunction([v],[z + (1-x)**2 - y])

solv = IpoptSolver(f,g)

class Log:
  def __init__(self):
    self.iter = 0 
  def __call__(self,f,*args):
    print "====Hey, I'm an iteration===="
    print "X_OPT = ", f.input(NLP_X_OPT)
    print f.getStats()
    self.iter = self.iter + 1
    if self.iter > 5:
      print "This is quite enough."
      f.output(0).set(1)

log = Log()

c = PyFunction( log, [sp_dense(3,1),sp_dense(1,1),sp_dense(1,1),sp_dense(3,1)], [sp_dense(1,1)] )
c.init()

solv.setOption("iteration_callback",c)

solv.init()

solv.setInput([2.5,3.0,0.75],NLP_X_INIT)
solv.setInput(0,NLP_UBG)
solv.setInput(0,NLP_LBG)
solv.solve()

print solv.output(NLP_X_OPT)

