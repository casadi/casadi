from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
from numpy.linalg import solve, norm
from numpy.random import rand
from numpy import *
from casadi import *

nrow = 5; ncol = 5
nnz = 12
  
rowind = (0, 3, 6, 8, 10, 12)
col = (0, 1, 4, 1, 2, 4, 0, 2, 0, 3, 3, 4)
s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;
val = (s, l, l, u, l, l, u, p, u, e, u, r)

# Solve with scipy
A = csr_matrix((val,col,rowind),shape=(nrow,ncol),dtype=float)
b = ones((nrow,1),dtype=float)
r_scipy = spsolve(A,b)

# Solve with SuperLU
linear_solver = SuperLU(nrow,ncol,rowind,col)
  
# Set options
linear_solver.setOption("colperm", "natural")

# Initialize
linear_solver.init()

# Pass Non-zero elements
linear_solver.setInput(val,0)

# Pass right hand side
linear_solver.setInput(nrow*[1.],1)

# Solve
linear_solver.evaluate()
  
# Print the solution
r_superlu = linear_solver.getOutput()
  
# print 
print "r_scipy = ", r_scipy
print "r_superlu = ", r_superlu
print "norm = ", norm(r_superlu-r_scipy)





#A = lil_matrix((1000, 1000))
#A[0, :100] = rand(100)
#A[1, 100:200] = A[0, :100]
#A.setdiag(rand(1000))

#A = A.tocsr()
#b = rand(1000)
#x = spsolve(A, b)

#x_ = solve(A.todense(), b)

#print "errnorm = ", norm(x-x_)


#from casadi import *
