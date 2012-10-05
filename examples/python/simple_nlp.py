from casadi import *

# Declare variables
x = ssym("x",2)

# Form the NLP objective
f = SXFunction([x],[x[0]**2 + x[1]**2])

# Form the NLP constraints
g = SXFunction([x],[x[0]+x[1]-10])

# Allocate a solver
#solver = IpoptSolver(f,g)
#solver = WorhpSolver(f,g)
solver = SQPMethod(f,g); solver.setOption("qp_solver",QPOasesSolver)
solver.init()

# Set constraint bounds
solver.setInput(0.,NLP_LBG)

# Solve the NLP
solver.evaluate()

# Print solution
print "-----"
print "objective at solution = ", solver.output(NLP_COST)
print "primal solution = ", solver.output(NLP_X_OPT)
print "dual solution (x) = ", solver.output(NLP_LAMBDA_X)
print "dual solution (g) = ", solver.output(NLP_LAMBDA_G)

