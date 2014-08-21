import casadi_core.*

x = SX.sym('x')

f = SXFunction({x},{cos(x)})
f.init()

f.setInput(3,0)

f.evaluate()


f.getOutput()

res = f.getOutput()-DMatrix(cos(3))
assert(res.isZero())
