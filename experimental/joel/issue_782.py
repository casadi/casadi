from casadi import *
from copy import copy

x = msym("x",2,1)

z=copy(x)
z[[1,0]]*=x

fun = MXFunction([x],[z])
fun.init()

for f,sym,Function,X in [(fun,msym,MXFunction,MX),]:
  f.init()
  print Function

  aseeds = sym("a",2)

  _,_,[[adjsens]] = f.eval([X.ones(2)],[],[[aseeds]])

  vf = Function([aseeds],[adjsens])
  vf.init()

  print vf
  
  vf.setInput([1,1])
  vf.setAdjSeed([0,1])
  vf.evaluate(0,1)
  print vf.getAdjSens()
  
  _,_,[[adjsens2]] = vf.eval([X.ones(2)],[],[[X([0,1])]])

  vf2 = Function([],[adjsens2])
  vf2.init()
  vf2.evaluate()

  print vf2.getOutput()
