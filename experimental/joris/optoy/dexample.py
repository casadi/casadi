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
