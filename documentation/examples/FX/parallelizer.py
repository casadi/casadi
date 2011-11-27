#! Demonstration of Parallelizer
#! =============================
from casadi import *

n = 12

#! Construct a functon that is expensive to evaluate
x = msym("x",100,100)
y = msym("y",100,1)
z = x

for i in range(20):
  z = mul(x,z)


f = MXFunction([x,y],[z])
f.init()

#! Evaluate this function ten times in parallel
p = Parallelizer([f]*n)
p.setOption("parallelization","openmp")
p.init()

#! Note that the parallelizer MX input/output interface is a repitition of our function's I/O interface
assert(p.getNumInputs() == n*f.getNumInputs())

p.evaluate(1,1)

print p.getStats()

