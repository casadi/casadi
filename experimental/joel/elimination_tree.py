# -*- coding: utf-8 -*-
from casadi import *
from numpy import *
import numpy as NP
import matplotlib.pyplot as plt

# Example matrix from Davis p. 39
A = DMatrix(11,11)
A[2,1] = 1
A[5,0] = 1
A[5,3] = 1
A[6,0] = 1
A[7,1] = 1
A[7,4] = 1
A[8,5] = 1
A[9,2] = 1
A[9,3] = 1
A[9,5] = 1
A[9,7] = 1
A[10,2] = 1
A[10,4] = 1
A[10,6] = 1
A[10,7] = 1
A[10,9] = 1

# Make symmetric
A = A+A.T

# Get the sparsity
s = A.sparsity()

# Get the elimination tree
etr = s.eliminationTree(False)
