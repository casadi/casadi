from casadi import *

x = msym("x",2,1)

fun = MXFunction([x],[diag(x[[1,0]])])
fun.init()

for f,sym,Function in [(fun,msym,MXFunction),]:
  f.init()
  print f

  print Function
  r = DMatrix([0,1])
  q = sym("q",2)

  f.setInput(r)
  f.setFwdSeed(r)
  f.evaluate(1,0)
  f.getFwdSens().printDense()

  _,[[fwdsens]],_ = f.eval([r],[[q]],[])
  vf = Function([q],[fwdsens])
  vf.init()
  vf.setInput(r)
  vf.evaluate()
  vf.getOutput().printDense()

