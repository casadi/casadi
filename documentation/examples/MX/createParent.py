#! createParent
#!======================
from casadi import *
from numpy import *

A = MX("A",2,1)         # Here a matrix
B = MX("B",1,2)         # There a matrix
C = MX("C")             # And an other little matrix
D = MX("D",sp_tril(4))  # Triangular matrix


L = [A,B,C,D]

for m in L:
  print m, " = ", array(m.sparsity())

V, (A,B,C,D) = createParent(L)

# The code below has been commented out (it relies on the mapping() function which no longer exists
#print V[A.mapping()]

#for m,n in zip(L,[A,B,C,D]):
  #print m, " => ", n, " = ", array(n.sparsity())

