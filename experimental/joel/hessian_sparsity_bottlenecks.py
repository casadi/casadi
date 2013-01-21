from casadi import *

# Create a function with a Hessian with one dense row, one dense column and a diagonal
n = 100000
x = ssym("x",n)
f = log(x[0])*inner_prod(x,x)
fcn = SXFunction([x],[f])
fcn.setOption("verbose",True)
fcn.init()

# Generate Hessian
hfcn = fcn.hessian()

# Visualize the sparsity pattern
#hfcn.init()
#hfcn.output().sparsity().spy()
