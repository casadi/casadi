#! CasADi tutorial
#! ==================
#! This tutorial file explains the interface with IPOPT
#! Ipopt solves problems of the form:
#!
#!
#! Minimize     f(x)
#! x in R^n
#! s.t         g_L <= g(x) <= g_U
#! x_L <= x    <= x_U
from numpy import *
import numpy as n
from casadi import *
#! Let's solve a simple scalar non-linear program:

x = SX("x")

y = SX("y")

#f  = SXFunction([x,y], tan(x)-1) 

#g=SXFunction([x,y],y)

#print f.eval(x)
#! Quadractic program
#! ------------------
#! (Ab)using Ipopt to do a simple quadratic problem
#!
#! 
P = n.eye(5)
A = n.diag([2.0,1,4,10,-2])
q = [1,0,0,0,0]
b = [1,1,1,1,1]

X = MX("x",5,1)
P = MX(DMatrix(P))
q = MX(DMatrix(q))
A = MX(DMatrix(A))

# Objective function
F = 0.5*mul(mul(trans(X),P),X) + mul(trans(q),X)

f = MXFunction([X],[F])
f.init()
f.setInput([1,1,1,1,1])
f.evaluate()
#! Test the objective for some value of x:
print f.output().toArray()

# constraint function
g = MXFunction([X],[X+X])
g.init()

solver = IpoptSolver(f,g)
solver.printOptions()

solver.init()

#! The default lower an upper bound on the optimizations variables is zero.
#! Change them to unbouded as follows:
solver.setInput([-100,-100,-100,-100,-100],NLP_LBX)
solver.setInput([100, 100, 100, 100, 100],NLP_UBX)

#! Inequality constraints.
#! The lower bound is also necessary, although Ipopt itself does not seem to require it
solver.setInput(b,NLP_UBG)
solver.setInput([-100,-100,-100,-100,-100],NLP_LBG)


solver.solve()
print solver.output(NLP_X_OPT)
#! Nested optimization
