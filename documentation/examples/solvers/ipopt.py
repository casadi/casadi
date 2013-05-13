#
#     This file is part of CasADi.
# 
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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
#! IpoptSolver
#! =====================
from casadi import *
from numpy import *

#! In this example, we will solve a few optimization problems with increasing complexity
#!
#! Scalar unconstrained problem
#! ============================
#$ $\text{min}_x \quad \quad (x-1)^2$ \\
#$ subject to $-10 \le x \le 10$ \\
#$
#$ with x scalar

x=ssym("x")
nlp=SXFunction(nlIn(x=x),nlOut(f=(x-1)**2))

solver = IpoptSolver(nlp)
solver.init()
solver.input("lbx").set([-10])
solver.input("ubx").set([10])
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
nlp=SXFunction(nlIn(x=x),nlOut(f=mul((x-1).T,x-1),g=vertcat([x[1]+x[2],x[0]])))

solver = IpoptSolver(nlp)
solver.init()
solver.input("lbx").set([-10]*n)
solver.input("ubx").set([10]*n)
#$  $ 2 \le x_0 \le 2$ is not really as bad it looks. Ipopt will recognise this situation as an equality constraint. 
solver.input("lbg").set([0,2])
solver.input("ubg").set([1,2])
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

x=ssym("x")
a=ssym("a")
a_ = 2
nlp=SXFunction(nlIn(x=x,p=a),nlOut(f=(x-a)**2))

solver = IpoptSolver(nlp)
solver.init()
solver.input("lbx").set([-10])
solver.input("ubx").set([10])
solver.input("p").set([a_])
solver.solve()

#! The solution is obviously a:
print solver.output()
assert(abs(solver.output()[0]-a_)<1e-9)

#! The parameter can change inbetween two solve calls:
solver.input("p").set([2*a_])
solver.solve()

#! The solution is obviously 2*a:
print solver.output()
assert(abs(solver.output()[0]-2*a_)<1e-9)

