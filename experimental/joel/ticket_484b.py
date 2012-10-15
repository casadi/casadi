from casadi import *

x = MX("x")
f = MXFunction([x],[x])
f.setOption("ad_mode","reverse")
#f.setOption("verbose",True)
f.init()
J = f.jacobian()
#J.setOption("verbose",True)
f.init()
J.init()
J.evaluate()
  
