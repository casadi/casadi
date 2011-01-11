from scipy.sparse import *
from casadi import *
from random import *
import time

n = 1000
m = 1000

# ~10 % non-zeros
nel = n*m/100

# tests
cases = ["casadi", "scipy"]

# number of repeats
rep = 4

for c in cases:
  print "testing ", c
  dur = []

  for i in range(rep):
    # create empty matrix
    if c=="casadi":
      M = matrix_double(n,m)
    else:
      M = csr_matrix((n,m))

    # fill with random elements
    t_before = time.time()
    for i in range(nel):
      row = randrange(n)
      col = randrange(m)
      res = float(row+col)
      M[row,col] = res
    t_after = time.time()
    dur.append(t_after-t_before)

  print "durations: ", dur
  print "number of non-zeros for the last matrix: ",
  if c=="casadi":
    print M.size()
  else:
    print M.getnnz()

