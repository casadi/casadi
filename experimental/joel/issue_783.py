from casadi import *

x = msym("x",2,1)

fun = MXFunction([x],[diag(x[[1,0]])])
fun.init()

for f,sym,Function,X in [(fun,msym,MXFunction,MX),]:
  f.init()
  print Function
  
  a = sym("a",sp_diag(2))

  _,_,[[f1]] = f.eval([X.ones(2)],[],[[a]])

  vf = Function([a],[f1])
  vf.init()
  print vf

  vf.setInput(1.0)
  vf.setAdjSeed([0,1])
  vf.evaluate(0,1)
  vf.getAdjSens().printDense()


  a2 = sym("a2",2)
  _,_,[[f2]] = vf.eval([X.ones(sp_diag(2))],[],[[a2]])

  vf2 = Function([a2],[f2])
  vf2.init()

  vf2.setInput([0,1])
  vf2.evaluate()
  vf2.getOutput().printDense()

