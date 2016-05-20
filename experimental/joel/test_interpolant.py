from casadi import *

# 1D
grid = [[0, 1, 2]]
values = [0, 1, 2]
print "grid = ", grid
print "values = ", values

F = interpolant('F', 'linear', grid, values)
print 'F(2.4) = ', F(2.4)
print 'F(1.4) = ', F(1.4)
print 'F(0.4) = ', F(0.4)
print 'F(-.6) = ', F(-.6)

# 2D
grid = [[0, 1, 2], [0, 1, 2]]
values = [0, 1, 2, 10, 11, 12, 20, 21, 22]
print "grid = ", grid
print "values = ", values

F = interpolant('F', 'linear', grid, values)
print 'F([2.4, 0.5]) = ', F([2.4, 0.5])
print 'F([1.4, 0.5]) = ', F([1.4, 0.5])
print 'F([0.4, 0.5]) = ', F([0.4, 0.5])
print 'F([-.6, 0.5]) = ', F([-.6, 0.5])
print 'F([-.6, 1.5]) = ', F([-.6, 1.5])
print 'F([-.6, 2.5]) = ', F([-.6, 2.5])
print 'F([-.6, 3.5]) = ', F([-.6, 3.5])

# Differentiate
gF = F.jacobian();
print 'gF([2.4, 0.5]) = ', gF([2.4, 0.5])
