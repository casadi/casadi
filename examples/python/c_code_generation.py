from casadi import *
from numpy import *
from os import system


x = SX("x"); y = SX("y")
fcn = SXFunction([[x,y]],[[x+y*sin(x)]])
fcn.init()

srcname = "fcn.c"
fcn.generateCode(srcname)

objname = "fcn.so"

print "Compiling..."
cmd = "gcc " + "-shared" + " -O3 "+ srcname + " -o " + objname
system(cmd)

# Read function
efcn = ExternalFunction(objname)
efcn.init()

for f in [fcn,efcn]:
  print "f.__class__ = ", f.__class__
  f.setInput([2,3])
  f.evaluate()
  print "result = ", f.output()
  
