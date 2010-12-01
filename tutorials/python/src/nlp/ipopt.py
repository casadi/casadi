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
from casadi import *
#! Let's solve a simple scalar non-linear program:

x = SX("x")

y = SX("y")

f  = SXFunction([x,y], tan(x)-1) 

g=SXFunction([x,y],y)

print f.eval(x)
#! Quadractic program
#! ------------------
#! (Ab)using Ipopt to do a simple quadratic problem
#!
#! 
P = eye(5)
A = diag([2.0,1,4,10,-2])
q = [1,0,0,0,0]
b = [1,1,1,1,1]

X = MX("x",5,1)
P = MX(list(P.ravel()),5,5)
q = MX(q,5,1)
A = MX(list(A.ravel()),5,5)

# Objective function
F = 0.5*prod(prod(trans(X),P),X) + prod(trans(q),X)

f = MXFunction(X,F)
f.setOption("ad_order",1)
f.init()
f.setInput([1,1,1,1,1])
f.evaluate()
#! Test the objective for some value of x:
print f.getOutput()

# constraint function
g = MXFunction(X,prod(A,X))
g.setOption("ad_order",1)
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
print solver.getOutput(NLP_X_OPT)
#! Nested optimization
