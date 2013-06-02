from casadi import *

in0 = [0.1,0.7,1.3]
in1 = [7.1,2.9]

x = ssym("x",3,1)
X = msym("x",3,1)

y = ssym("y",2,1)
Y = msym("y",2,1)

f = SXFunction([x,y],[y,x*y[0]])
F = MXFunction([X,Y],[Y,X*Y[0]])

f.init()
F.init()

f.setInput(in0,0)
F.setInput(in0,0)

f.setInput(in1,1)
F.setInput(in1,1)

f.setFwdSeed(0,0)
F.setFwdSeed(0,0)
f.setFwdSeed([1,0],1)
F.setFwdSeed([1,0],1)

f.setAdjSeed(0,0)
F.setAdjSeed(0,0)
f.setAdjSeed([1,0,0],1)
F.setAdjSeed([1,0,0],1)

f.evaluate(1,1)
F.evaluate(1,1)

print f.getAdjSens(1),F.getAdjSens(1) 
print f.getFwdSens(1),F.getFwdSens(1) 


