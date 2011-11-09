#! SXFunction constructors
#! =======================
from casadi import *

x = SX("x") # A scalar symbolic
y = ssym("y",2,1) # A matrix symbolic


ins = [x,y] # function inputs
outs = [x,y,[[x,x],[x,x]],y*x,0]

print outs

f = SXFunction(ins,outs)
f.init()

#! f now has two inputs and a 4 outputs:
print f.getNumInputs()
print f.getNumOutputs()

#! The outputs has the following string representation.
#! Note how all elements of out have been converted to SXMatrix by
#! automatic typecasting functionality

for i in range(3):
  print f.outputSX(i)
