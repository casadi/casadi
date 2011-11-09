#! IpoptSolver
#! =====================
from casadi import *
from numpy import *
import casadi as c

#! In this example, we will solve a few optimization problems with increasing complexity
#!
#! Scalar unconstrained problem
#! ============================
#$ $\text{min}_x \quad \quad (x-1)^2$ \\
#$ subject to $-10 \le x \le 10$ \\
#$
#$ with x scalar

x=SX("x")
f=SXFunction([x],[[(x-1)**2]])

solver = IpoptSolver(f)
solver.init()
solver.input(NLP_LBX).set([-10])
solver.input(NLP_UBX).set([10])
solver.solve()

#! The solution is obviously 1:
print solver.output()
assert(abs(solver.output()[0]-1)<1e-9)

#! Constrained problem
#! ============================
#$ $\text{min}_x \quad \quad (x-1)^T.(x-1)$ \\
#$ subject to $-10 \le x \le 10$ \\
#$ subject to $0 \le x_1 + x_2 \le 1$ \\
#$ subject to $ x_0 = 2$ \\
#$
#$ with $x \in \mathbf{R}^n$

n = 5

x=ssym("x",n)
#! Note how we do not distinguish between equalities and inequalities here
f=SXFunction([x],[[c.dot((x-1).T,x-1)]])
g=SXFunction([x],[[x[1]+x[2],x[0]]])

solver = IpoptSolver(f,g)
solver.init()
solver.input(NLP_LBX).set([-10]*n)
solver.input(NLP_UBX).set([10]*n)
#$  $ 2 \le x_0 \le 2$ is not really as bad it looks. Ipopt will recognise this situation as an equality constraint. 
solver.input(NLP_LBG).set([0,2])
solver.input(NLP_UBG).set([1,2])
solver.solve()

#! The solution is obviously [2,0.5,0.5,1,1]:
print solver.output()
for (i,e) in zip(range(n),[2,0.5,0.5,1,1]):
  assert(abs(solver.output()[i]-e)<1e-7)


#! Problem with parameters
#! ============================
#$ $\text{min}_x \quad \quad (x-a)^2$ \\
#$ subject to $-10 \le x \le 10$ \\
#$
#$ with x scalar

x=SX("x")
a=SX("a")
a_ = 2
f=SXFunction([x,a],[[(x-a)**2]])

solver = IpoptSolver(f)
solver.init()
solver.input(NLP_LBX).set([-10])
solver.input(NLP_UBX).set([10])
solver.input(NLP_P).set([a_])
solver.solve()

#! The solution is obviously a:
print solver.output()
assert(abs(solver.output()[0]-a_)<1e-9)

#! The parameter can change inbetween two solve calls:
solver.input(NLP_P).set([2*a_])
solver.solve()

#! The solution is obviously 2*a:
print solver.output()
assert(abs(solver.output()[0]-2*a_)<1e-9)

