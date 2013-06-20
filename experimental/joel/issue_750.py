from casadi import *

x = msym("x",2,1)

y = x[0]+x[1]

f=MXFunction([x],[y])
f.setOption("ad_mode","forward")

print "These should be the same"
for mode in [True,False]:
  f.setOption("numeric_jacobian",mode)
  f.init()

  Jf=f.jacobian(0,0)
  Jf.init()
  Jf.input().set([1.2,2.3])
  Jf.evaluate()

  print Jf.output()
