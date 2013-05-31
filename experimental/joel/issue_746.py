from casadi import *

x = msym("x",2)

f1 = MXFunction([x],[x[[1,0]]])
f1.init()

f2 = MXFunction([x],[vertcat([x[1],x[0]])])
f2.init()

num = DMatrix([1.1,1.3])
fseed_ = [1,0]

for f in f1,f2:
  f.setInput(num)
  f.setFwdSeed(fseed_)
  f.evaluate(1,0)

print "These should be the same: "
print f1.output(), f2.output()

print "These should be the same: "
print f1.fwdSens(), f2.fwdSens()

inputs = [msym("x",2)]
fseeds = [[msym("fseed",2)]]
aseeds = [[msym("aseed",2)]]


print "These should be the same: "
for f in f1,f2:
  res,fwdsens,_ = f.eval(inputs,fseeds,[])
  fs = MXFunction(inputs+fseeds[0],fwdsens[0])
  fs.init()

  fs.setInput(num,0)
  fs.setInput(fseed_,1)

  fs.evaluate()

  print fs.output()
