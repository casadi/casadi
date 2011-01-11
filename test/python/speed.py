from scipy.sparse import *
from casadi import *
from random import *
import time
import numpy as N

n = 1000
m = 1000

# ~10 % non-zeros
nel = n*m/10

# tests
solvers = ["scipy","casadi"]

# build with linked lists
use_lil = True

# number of repeats
rep = 4

for sol in solvers:
  print "testing ", sol
  dur_build = []
  dur_convert = []

  for i in range(rep):
    # create empty matrix
    if sol=="casadi":
      M = matrix_double(n,m)
    else:
      if use_lil:
        M = lil_matrix((n,m))
      else:
        M = csr_matrix((n,m))

    # fill with random elements
    t_before = time.time()
    for i in range(nel):
      row = randrange(n)
      col = randrange(m)
      res = float(row+col)
      M[row,col] = res

    # need to convert back to csr
    if sol!="casadi" and use_lil:
      M = csr_matrix(M)
      
    t_build = time.time()
    dur_build.append(t_build-t_before)

    # convert casadi matrix to csr
    if sol=="casadi":
      M = csr_matrix((M,M.col(),M.rowind()),(M.size1(),M.size2()),dtype=float)
    
    # Save the matrix
    scipyM = M

    # convert scipy to casadi
    dim = M.shape
    col = list(int(i) for i in M.indices)
    rowind = list(int(i) for i in M.indptr)
    M = matrix_double(dim[0],dim[1],col,rowind,M.data)
    
    # Save the matrix
    casadiM = M

    # convert back again
    if sol!="casadi":
      M = csr_matrix((M,M.col(),M.rowind()),(M.size1(),M.size2()),dtype=float)

    t_convert = time.time()
    dur_convert.append(t_convert-t_build)


    ## Make a couple of operations
    #if sol=="casadi":
      #M2 = trans(M)
      ##M2 = prod(M,M)
      ###M = M2 + M
    #else:
      #M2 = M.transpose()
      ##M2 = N.dot(M,M)
      ##M = M2 + M       
    ##M = M2
      
    

  print "durations to build:                      ", dur_build
  print "durations for converting back and forth: ", dur_convert
  print "number of non-zeros for the last matrix: ",
  if sol=="casadi":
    print M.size()
  else:
    print M.getnnz()

