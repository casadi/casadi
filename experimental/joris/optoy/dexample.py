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
from pylab import *
from casadi import *
from doptoy import *

# VDP regulate to zero
x = state()
y = state()
q = state()

u = control()

x.dot = (1-y**2)*x-y+u
y.dot = x
q.dot = x**2+y**2

ocp(q.end,[u>=-1,u<=1,q.start==0,x.start==1,y.start==0],T=10,N=20)

print x.sol

plot(x.sol)
plot(y.sol)
plot(u.sol)

show()


# VDP time-optimal 
x = state()
y = state()
q = state()

u = control()

T = var(lb=0,init=10)

x.dot = (1-y**2)*x-y+u
y.dot = x
q.dot = x**2+y**2

ocp(T,[u>=-1,u<=1,q.start==0,x.start==1,y.start==0,x.end==0,y.end==0],T=T,N=20)

print T.sol, x.sol
plot(x.sol)
plot(y.sol)
plot(u.sol)

show()
