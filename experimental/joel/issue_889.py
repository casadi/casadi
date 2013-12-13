from casadi import *
#A = DMatrix.ones(3,3) + DMatrix.ones(sp_diag(3)) # dense
#A = DMatrix.ones(sp_diag(3)) # diagonal
A = 2*DMatrix.ones(sp_diag(3)); A[1,0] = 1; A[2,0] = 1; A[2,1] = 1 # lower triangular

print "A = "
A.printDense()

print "MX"
r = msym("r",3)
sol = CSparse(A.sparsity())
sol.init()
x = sol.solve(A,r.T,True).T
f = MXFunction([r],[x])
f.init()
DMatrix.ones(f.jacSparsity()).printDense()

print "SX"
r = ssym("r",3)
x = solve(A,r)
f = SXFunction([r],[x])
f.init()
DMatrix.ones(f.jacSparsity()).printDense()

