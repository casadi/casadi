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
from casadi import *
from casadi.tools import *

x = ssym("[x y]",2)
dx = ssym("[dx dy]",2)
u = ssym("[u v]",2)
du = ssym("[du dv]",2)
ddx = ssym("[ddx ddy]",2)

q = vertcat([x,u])
dq = vertcat([dx,du])
ddq = ssym("[ddx,ddy,ddu,ddv]",4)

[m,g,L] = ssym("[m,g,L]",3)

length = L
gravity=g
T = 0.5*m*mul(dx.T,dx)     # Kinetic    energy
V = m*g*x[1]               # Potential  energy
C = 0.5*(mul(x.T,x)-L**2)  # Constraint energy

lambd = ssym("lambda")     # Constraint multiplier

L = T - V - lambd*C        # Lagrange function




Ldx = jacobian(L,dx)       # partial dL/dx

 # A function to take a time-derivative
dt = lambda L: mul(jacobian(L,x),dx) + mul(jacobian(L,dx),ddx)  


print "L = ", L

print "Some intermediary results:"
print "Ldx = diff(L,dx)   =   ", Ldx
print "diff(L,x)          =   ", jacobian(L,x)
print "diff(Ldx,t)        =   diff(Ldx,x)*dx + diff(Ldx,dx)*ddx"
print "diff(Ldx,x)        =   ", jacobian(Ldx,x)
print "diff(Ldx,dx)       =   ", jacobian(Ldx,dx)

print ""


print "Finally Lagrange provides us with these equations:"

U = dt(Ldx) - jacobian(L,x).T

print U

print 
print "But we want a first-order ODE, so we will introduce helper variables " ,u


U = vertcat([substitute(U,vertcat([dx,ddx]),vertcat([u,du])),u-dx])

g = C

print " DAE in fully implicit form:"
print "  0 = U(y,y',z,t)"
print "  0 = g(y,z,t)"
print ""
print " U: ", U
print " g: ", g

print ""
print ""

print "DAE index is a concept defined for semi-explicit DAE's:"
print "  y' = f(y,z,t)"
print "  0  = g(y,z,t)"

print ""
print "  We can find f as jacobian(U,y')^-1*U(y,0,z,t)"
print "   U(y,0,z,t) : " , substitute(U,dq,[0,0,0,0])
print "   jacobian(U,y') : " ,    jacobian(U,dq)
print "   jacobian(U,y')^-1 : " , inv(jacobian(U,dq))
print ""
f =  -mul(inv(jacobian(U,dq)),substitute(U,dq,[0,0,0,0]))


print "This DAE has a y' = f(y,z,t) part:"
print "f : " , f
print "g : " , C
print "y : ", q
print "z : ", lambd


print ""
print ""
print "Will determine index of DAE"
print "Start out with g as before"

g = C


 # A function to take a time-derivative
dt = lambda L: mul(jacobian(L,q),dq) + mul(jacobian(L,dq),ddq)  


for i in range(1,4):
  print " Is it index-%d ?" %  i
  print " ==================="
  print "  Is J=jacobian(g,z) invertible?"
  
  J = jacobian(g,lambd)
  print "  J =", J
  print ""
  if J.size()==1:
    print "  Yes. You are finished"
    break
  else:
    print "  No. Find a new g:    g ->    diff(g)  & substitute y' = f(y,z,t) "

    print "  This is diff(g): ", dt(g)
    g =  substitute(dt(g),dq,f)
    print "  This is after substitution: ", g
    
  print ""
  
  
print ""
sys = vertcat([U,g])
print "So we end up with this %s residual: " % (str(sys.shape))

print sys

print "But wait, there's too much symbols in there. This is not a closed system!"
print "We are screwed"

print "Check that there's no second derivatives left: ", jacobian(sys,ddq).size()==0

dae = SXFunction({'NUM': DAE_NUM_IN, DAE_Y: vertcat([q,lambd]), DAE_YDOT: vertcat([dq,SX("dlambda")]), DAE_P: vertcat([length,m,gravity])},[sys])
dae.init()

p_ = [5,1,10]
x_ = [1,2,0,0,0]
dx_ = [0,0,0,0,0]

J = SXFunction

dae.setInput(x_,"y")
dae.setInput(dx_,"ydot")
dae.setInput(p_,"p")

dae.evaluate()
print "res @ inital consitions: ", dae.getOutput()

integr = IdasIntegrator(dae)
#integr.setOption('is_differential',[1]*4 + [0])
integr.init()
integr.setInput(x_,"x0")
integr.setInput(p_,"p")
integr.setInput(dx_,"xp0")
integr.evaluate()

print integr.getOutput()

