from casadi import *
from optoy import *

# Simple unconstrained rosenbrock problem
x = var()
y = var()

cost = minimize((1-x)**2+100*(y-x**2)**2)
print "cost = ", cost
print "sol = ", x.sol, y.sol

# Bounds
x.lb = 2
y = var(ub=6)

print "cost = ", minimize((1-x)**2+100*(y-x**2)**2)
print "sol = ", x.sol, y.sol

# Constraints
x = var()
y = var()

print "cost = ", minimize((1-x)**2+100*(y-x**2)**2,[x**2+y**2<=1, x+y>=0])
print "sol = ", x.sol, y.sol

# Matrix symbols
xy = var(2)

print "cost = ", minimize((1-xy[0])**2+100*(xy[1]-xy[0]**2)**2)
print "sol = ", xy.sol

# Adding a parameter
p = par()
p.value = 100

print "cost = ", minimize((1-x)**2+p*(y-x**2)**2)
print "sol = ", x.sol, y.sol

