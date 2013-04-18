from casadi import *
x = ssym("x",2)
f = SXFunction([x],[x[0]*(3+x[1]**2)])
s = IpoptSolver(f)
s.init()
s.setInput([1.,    -10.],"lbx")
s.setInput([1.,     10.],"ubx")
s.evaluate()
print s.output("lam_x") # [0,0] # wrong!!!

# Change upper bound from an equality to an inequality constraint
s.setInput([1.001,  10.],"ubx")
s.evaluate()
print s.output("lam_x") # [-3,0] # correct!!!
