from casadi import *

x = ssym("x")
y = ssym("y")
f=SXFunction([y,x],[x-arcsin(y)])
f.init()
solver=NewtonImplicitSolver(f)
solver.setOption("linear_solver",CSparse)
solver.init()
solver.setAdjSeed(1.0)
solver.setInput(0.2)
solver.setOutput(0.1)
solver.evaluate(0,1)

print "sin(0.2) = ", sin(0.2)
print "cos(0.2) = ", cos(0.2)
print "y = ", solver.output(0)
print solver.adjSens(0)
