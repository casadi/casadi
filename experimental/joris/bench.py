from casadi import *
from time import time
from numpy import *
import casadi as c
import numpy as n


N = 5000

# We work with a sparse diaginal matrix with 5000 elements
x = DMatrix(sp_diag(N),5)

x_s = x.toCsr_matrix()

t = time()
dummy = n.dot(x_s.T,x_s)
print "Scipy pure = %.4f s" % (time()-t)

t = time()
dummy = c.dot(x.T,x)
print "CasADi pure = %.4f s" % (time()-t)

X = MX("X",x.sparsity())

t = time()
f = MXFunction([X],[c.prod(X.T,X)])
f.init()
print "CasADi MX wrapped init overhead = %.4f s" % (time()-t)
f.input().set(x)

t = time()
f.evaluate()
print "CasADi MX wrapped = %.4f s" % (time()-t)


