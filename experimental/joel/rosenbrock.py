from casadi import *
import numpy as NP
import matplotlib.pyplot as plt

x = ssym("x")
y = ssym("y")
z = ssym("z")
v = vertcat([x,y,z])

f = SXFunction([v],[x**2 + 100*z**2])
g = SXFunction([v],[z + (1-x)**2 - y])

#solv = IpoptSolver(f,g)

solv = SQPMethod(f,g)
solv.setOption("qp_solver",IpoptQPSolver)

solv.init()

solv.setInput([2.5,3.0,0.75],NLP_X_INIT)
solv.setInput(0,NLP_UBG)
solv.setInput(0,NLP_LBG)
solv.solve()

print solv.output(NLP_X_OPT)


