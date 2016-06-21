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
#! nlpsol
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

x=SX.sym('x')
nlp = {'x':x, 'f':(x-1)**2}

solver = nlpsol('solver', 'ipopt', nlp)
sol = solver(lbx=-10, ubx=10)

#! The solution is obviously 1:
print(sol['x'])
assert(abs(sol['x']-1)<1e-9)

#! Constrained problem
#! ============================
#$ $\text{min}_x \quad \quad (x-1)^T.(x-1)$ \\
#$ subject to $-10 \le x \le 10$ \\
#$ subject to $0 \le x_1 + x_2 \le 1$ \\
#$ subject to $ x_0 = 2$ \\
#$
#$ with $x \in \mathbf{R}^n$

n = 5

x=SX.sym('x',n)
#! Note how we do not distinguish between equalities and inequalities here
nlp = {'x':x, 'f':mtimes((x-1).T,x-1), 'g':vertcat(x[1]+x[2],x[0])}

solver = nlpsol('solver', 'ipopt', nlp)
sol = solver(lbx=-10, ubx=10, lbg=[0,2], ubg=[1,2])
#$ $ 2 \le x_0 \le 2$ is not really as bad it looks. 
#$ Ipopt will recognise this situation as an equality constraint.

#! The solution is obviously [2,0.5,0.5,1,1]:
print(sol['x'])
for (i,e) in zip(list(range(n)),[2,0.5,0.5,1,1]):
  assert(abs(sol['x'][i]-e)<1e-7)


#! Problem with parameters
#! ============================
#$ $\text{min}_x \quad \quad (x-a)^2$ \\
#$ subject to $-10 \le x \le 10$ \\
#$
#$ with x scalar

x=SX.sym('x')
a=SX.sym('a')
a_ = 2
nlp={'x':x, 'p':a, 'f':(x-a)**2}

solver = nlpsol('solver', 'ipopt', nlp)
sol = solver(lbx=-10, ubx=10, p=a_)

#! The solution is obviously a:
print(sol['x'])
assert(abs(sol['x']-a_)<1e-9)

#! The parameter can change inbetween two solve calls:
sol = solver(lbx=-10, ubx=10, p=2*a_)

#! The solution is obviously 2*a:
print(sol['x'])
assert(abs(sol['x']-2*a_)<1e-9)
