from __future__ import division
from casadi import *
from numpy import *

# Matrix multiplication
# Given the SXMatrices A=[reA,imA] and B=[reB, imB], immul(A,B) returns an
# Mat C=[reC,imC] where C = A*B (matrix product).

def immul(A,B):
    return [mul(A[0],B[0]) - mul(A[1],B[1]),mul(A[0],B[1]) + mul(A[1],B[0])]

# Transpose of A=[reA,imA]:
def imtrans(A):
    return [A[0].T,A[1].T]

# Complex inverse of A=[reA,imA]
''' complex inverse
From paper we have
A = R + jS
K2 is diagonal e.g. I
X = -R*(S+K2*R)^-1
Y = I + X*K2
Minv = (X*R - Y*S)^-1 *(X +jY)
'''
def iminv(A):
    K2 = Mat.eye(len(A[0]))
    X = -mul(A[0],linalg.inv(A[1]+mul(K2,A[0])))
    Y = Mat.eye(len(A[0])) + mul(X,K2)
    tmp = linalg.inv(mul(X,A[0])-mul(Y,A[1]))
    return [mul(tmp,X),mul(tmp,Y)]

# Construct an imaginary Mat ([re,im]) from D = E*exp(j*angle)
def imexp(abs,angle):
    re = abs*numpy.cos(angle)
    im = abs*numpy.sin(angle)
    return [re,im]

# Conjugate for imaginary Mat
def imconj(A):
    return [A[0], -A[1]]

# Elementwise multiplication of imaginary matrices
def imelemmul(A,B):
    Are = A[0]
    Aim = A[1]
    Bre = B[0]
    Bim = B[1]
    [rAre,cAre] = shape(Are)
    [rBre,cBre] = shape(Bre)
    [rAim,cAim] = shape(Aim)
    [rBim,cBim] = shape(Bim)

    return [np.multiply(A[0],B[0]) - np.multiply(A[1],B[1]), np.multiply(A[0],B[1]) + np.multiply(A[1],B[0])]












