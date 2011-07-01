#! Exact Hessian
#! =====================
from casadi import *
from numpy import *
import casadi as c

#! We will investigate the use of an exact Hessian with the help of the Rosenbrock function
x=SX("x")
y=SX("y")
    
obj = (1-x)**2+100*(y-x**2)**2
#! We choose to add a single constraint
constr = x**2+y**2

f=SXFunction([[x,y]],[obj])
g=SXFunction([[x,y]],[constr])
solver = IpoptSolver(f,g)
    
#! We need the hessian of the lagrangian.
#! A problem with n decision variables and m constraints gives us a hessian of size n x n
  
sigma=SX("sigma")  # A scalar factor
lambd=SX("lambd")  # Multipier of the problem, shape m x 1.


h=SXFunction([[x,y],[lambd],[sigma]],[sigma*hessian(obj,[x,y])+lambd*hessian(constr,[x,y])])
   
#! We solve the problem with an exact hessian
solver = IpoptSolver(f,g,h)
solver.init()
solver.input(NLP_LBX).set([-10]*2)
solver.input(NLP_UBX).set([10]*2)
solver.input(NLP_LBG).set([0])
solver.input(NLP_UBG).set([1])
solver.solve()

for sol in array(solver.output()):
  print "%.15f" % sol

#! To compare the behaviour of convergence, we solve the same problem without exact hessian
solver = IpoptSolver(f,g)
solver.init()
solver.input(NLP_LBX).set([-10]*2)
solver.input(NLP_UBX).set([10]*2)
solver.input(NLP_LBG).set([0])
solver.input(NLP_UBG).set([1])
solver.solve()

for sol in array(solver.output()):
  print "%.15f" % sol

