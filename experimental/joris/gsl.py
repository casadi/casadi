from casadi import *

x=SX("x")
y=SX("y")

f = SXFunction({'NUM': DAE_NUM_IN, DAE_Y:[x,y]},[[y,-x]])
f.init()

I = GslIntegrator(f)
I.setOption('tf',0.1)
I.init()
I.input(INTEGRATOR_X0).set([1,0])
I.evaluate()

print I.output()
