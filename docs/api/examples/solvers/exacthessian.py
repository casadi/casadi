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
#! Exact Hessian
#! =====================
from casadi import *
from numpy import *
import casadi as c

#! We will investigate the use of an exact Hessian with the help of the Rosenbrock function
x=SX.sym('x')
y=SX.sym('y')
obj = (1-x)**2+100*(y-x**2)**2
constr = x**2+y**2
nlp={'x':vertcat(x,y), 'f':obj, 'g':constr}

#! We solve the problem with an exact Hessian (default)
solver = nlpsol('solver', 'ipopt', nlp)
sol = solver(lbx=-10, ubx=10, lbg=0, ubg=1)
print('Optimal solution (exact Hessian): %s' % sol['x'])

#! Same problem but with limited memory BFSG
solver = nlpsol('solver', 'ipopt', nlp, {'ipopt.hessian_approximation':'limited-memory'})
sol = solver(lbx=-10, ubx=10, lbg=0, ubg=1)
print('Optimal solution (BFGS): %s' % sol['x'])
