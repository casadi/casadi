from casadi import *
from numpy import *
from os import system
import time

x = symbolic("x",7,7)
f = det(x)
x = vec(x)
x0 = [random.rand() for xi in x.data()]

fcn = SXFunction([x],[[f]])
fcn.init()

# adjoint
gf = fcn.grad()

gfcn = SXFunction([x],[gf])
gfcn.init()

srcname = "grad_det.c"
gfcn.generateCode(srcname)

objname_no_opt = "grad_det_no_opt.so"
print "Compiling without optimization: ", objname_no_opt
system("gcc -shared " + srcname + " -o " + objname_no_opt)

objname_O3_opt = "grad_det_O3_opt.so"
print "Compiling with O3 optimization: ", objname_O3_opt
system("gcc -shared -O3 " + srcname + " -o " + objname_O3_opt)

# Read function
efcn_no_opt = ExternalFunction(objname_no_opt)
efcn_O3_opt = ExternalFunction(objname_O3_opt)
efcn_no_opt.init()
efcn_O3_opt.init()

for f in [gfcn,efcn_no_opt,efcn_O3_opt]:
  print "f.__class__ = ", f.__class__
  f.setInput(x0)
  t1 = time.time()
  nrep = 10000
  for r in range(nrep):
    f.evaluate()
  t2 = time.time()
  print "result = ", f.output().data()
  dt = (t2-t1)/nrep
  print "time = ", dt*1e3, " ms"
  
  num_op = len(gfcn.algorithm())
  print "number of elementary operations: ", num_op
  print "time per elementary operations: ", dt/num_op*1e9, " ns"
  
  
