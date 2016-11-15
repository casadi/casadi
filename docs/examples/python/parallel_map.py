import numpy as np
import time
print("importing casadi...")
from casadi import *


# number of inputs to evaluate in parallel
N = 300


# dummy input
dummyInput = np.linspace(0.0, 2.0*np.pi, N)


# make a dummy function that's moderately expensive to evaluate
print("creating dummy function....")
x = SX.sym('x')
y = x
for k in range(100000):
    y = sin(y)
f0 = Function('f', [x], [y])


# evaluate it serially, the old-fasioned way
X = MX.sym('x',N)
Y = vertcat(*[f0(X[k]) for k in range(N)])
fNaiveParallel = Function('fParallel', [X], [Y])

print("evaluating naive parallel function...")
t0 = time.time()
outNaive = fNaiveParallel(dummyInput)
t1 = time.time()
print("evaluated naive parallel function in %.3f seconds" % (t1 - t0))


# evaluate it using new serial map construct
fMap = f0.map(N)

print("evaluating serial map function...")
t0 = time.time()
outMap = fMap(dummyInput)
t1 = time.time()
print("evaluated serial map function in %.3f seconds" % (t1 - t0))
# the following has different shaped outputs, so it's commented out
#print outNaive == outMap


# evaluate it using new parallel map construct
fMap = f0.map(N, "openmp")

print("evaluating parallel map function...")
t0 = time.time()
outMap = fMap(dummyInput)
t1 = time.time()
print("evaluated parallel map function in %.3f seconds" % (t1 - t0))
# the following has different shaped outputs, so it's commented out
#print outNaive == outMap
