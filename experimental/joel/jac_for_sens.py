from casadi import *
from numpy import *
a = ssym("a",5)
b = ssym("b",5)

for jac_for_sens in (True,False):
  print "jac_for_sens = ", jac_for_sens, ":"
  f = SXFunction([a,b],[sqrt(b-sin(a)),inner_prod(a,b),outer_prod(a,b)])
  f.setOption("jac_for_sens",jac_for_sens)
  f.init()
  #print f.jacSparsityOld(0,0)

  f.setInput([1,2,3,4,5],0)
  f.setInput([10,20,30,40,50],1)

  f.setFwdSeed([1,0,0,0,1],0)
  f.setFwdSeed([0,0,0,0,0],1)

  f.setAdjSeed([0,0,0,0,0],0)
  f.setAdjSeed(1,1)
  f.setAdjSeed(DMatrix(5,5,0),2)

  f.evaluate(1,1)

  print "f.fwdSens(0) = ", f.fwdSens(0).data()
  print "f.fwdSens(1) = ", f.fwdSens(1).data()
  print "f.fwdSens(2) = ", f.fwdSens(2).data()
  print "f.adjSens(0) = ", f.adjSens(0).data()
  print "f.adjSens(1) = ", f.adjSens(1).data()

  print "--"