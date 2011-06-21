#! createParent
#!======================
from casadi import *
from numpy import *

A = MX("A",2)           # Here a matrix
B = MX("B",2,1)         # There a matrix
C = MX("C")             # And an other little matrix
D = MX("D",sp_tril(4))  # Triangular matrix


L = veccat([A,B,C,D])
print L

#! The following statement does the same:
L = map(vecNZ,[A,B,C,D])
