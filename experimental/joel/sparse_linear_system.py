# -*- coding: utf-8 -*-
from casadi import *
from numpy import *
import matplotlib.pyplot as plt

A = array( [[ 0.0,  0.0, 6.0, 7.0, 8.0, 0.0 ],
            [ 0.0,  0.0, 0.0, 9.0,10.0,11.0 ],
            [ 0.0,  1.0, 2.0, 0.0, 0.0, 0.0 ],
            [ 14.0, 0.0, 0.0, 0.0,12.0,13.0 ],
            [  0.0, 3.0, 4.0, 5.0, 0.0, 0.0 ],
            [ 16.0, 0.0, 0.0, 0.0, 0.0,15.0 ]])

b = array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])

x = linalg.solve(A,b)
print "x  = ", x

cA = DMatrix(A)
makeSparse(cA)

cb = DMatrix(b)
cx = casadi.solve(cA,cb)
print "cx = ", array(trans(cx))[0]

A = array( [[ 0.0,  1.0 ],
            [ 1.0,  0.0 ]] )
b = array([1.0, 2.0])

x = linalg.solve(A,b)
print "x  = ", x

cA = DMatrix(A)
makeSparse(cA)

cb = DMatrix(b)
cx = casadi.solve(cA,cb)
print "cx = ", array(trans(cx))[0]
