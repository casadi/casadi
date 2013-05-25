
from casadi import *
from copy import copy

x = msym("x",2)


w = copy(x)
w[[0,0,1,1]]+=2

wr = copy(x)
wr[[0,1]]+=2

f1 = MXFunction([x],[w])
f1.init()

f2 = MXFunction([x],[wr])
f2.init()

num = DMatrix([1.1,1.3])
fseed_ = [1,0]

for f in f1,f2:
  f.input().set(num)
  f.fwdSeed().set(fseed_)
  f.evaluate(1,0)

print "These should be the same: "
print f1.output(), f2.output()

print "These should all be the same: "
print f1.fwdSens(), f2.fwdSens()

inputs = [msym("x",2)]
fseeds = [[msym("fseed",2)]]

for f in f1,f2:
  res,fwdsens,_ = f.eval(inputs,fseeds,[])
  fs = MXFunction(inputs+fseeds[0],fwdsens[0])
  fs.init()

  fs.input(0).set(num)
  fs.input(1).set(fseed_)

  fs.evaluate()

  print fs.output()
