#! createParent
#!======================
from casadi import *
from numpy import *
import casadi as c

# Commented out since it relies on the now removed mapping method

##! We create a bunch of matrices that depend on V
#V = MX("V",13)
#A = c.reshape(V[:2],sp_dense(2,1))  # 2x1 matrix
#B = c.reshape(V[2:3],sp_dense(1,1)) # 1x1 matrix
#C = c.reshape(V[3:],sp_tril(4))     # Triangular matrix

#print "A = ", A
#print "A.mapping() = ", A.mapping()
#print "B = ", B
#print "B.mapping() = ", B.mapping()
#print "C = ", C
#print "C.mapping() = ", C.mapping()

##! We create a bunch of matrices that depend on V
#V_= DMatrix(13,1,True)   # create a dense DMatrix

##! We can use mapping() to index into this DMatrix:
#V_[A.mapping()] = 1
#V_[B.mapping()] = 10
#V_[C.mapping()] = 100
  
#print V_
