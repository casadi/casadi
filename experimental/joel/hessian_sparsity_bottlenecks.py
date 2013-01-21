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


# Note about the bipartitive graph defined by the Hessian sparsity pattern:
#
# First row has:
#   O(n) distance-1 neighbors
#   O(n^2) distance-2 neighbors
#   O(n^3) distance-3 neighbors
#
# Other rows have:
#   O(1) distance-1 neighbor
#   O(n) distance-2 neighbors
#   O(n^2) distance-3 neighbors

# Averarge number of distance-3 neighbors: O(n^2)   # Cf. Gebremedhin2005, StarColoringAlg1
# Averarge number of distance-2 neighbors: O(n),    # Cf. Gebremedhin2005, StarColoringAlg2
