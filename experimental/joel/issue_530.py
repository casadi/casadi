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

f.input(0).set(in0)
F.input(0).set(in0)

f.input(1).set(in1)
F.input(1).set(in1)

f.adjSeed(1).set([1,0,0])
F.adjSeed(1).set([1,0,0])

f.evaluate(0,1)
F.evaluate(0,1)

print f.adjSens(1),F.adjSens(1) 


