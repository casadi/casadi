#! CasADi tutorial 1
#! ==================
#! This tutorial file explains the use of CasADi's Matrix<T> in a python context.
#! Matrix<T> is a general class for sparse matrices. We inspect it with the help of Matrix<double>
#! Let's start with the import statements to load CasADi.
from casadi import *
from numpy import *
#! Contructors & printing
#! --------------------------------------
#! The python name for Matrix<double> is DMatrix
a = DMatrix(3,4)
print a
#! The string representation shows only the structural non-zero entries. In this case there are none.
#! Let's make a DMatrix with some structural non-zero entries.
w = DMatrix(3,4,[1,2,1],[0,2,2,3],[3,2.3,8])
print w
#! Conversion
#! --------------
#! DMatrix can easily be converted into other data formats
print list(w)
print tuple(w)
print w.toArray()
print w.toCsr_matrix()
