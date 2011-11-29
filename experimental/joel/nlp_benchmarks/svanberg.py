# -*- coding: utf-8 -*-
from casadi import *
import numpy as NP
import matplotlib.pyplot as plt

# This is a CasADi version of svanberg.mod from the cute test collection, original by Hande Y. Benson
#   
#   Svanberg K.,
#   "Method of moving asymptots - a new method for structural optimization",
#   Int.J. Num. Meth. Eng, 24, pp. 359--373, 1987

n = 5000


#param b{i in 1..n} := i*5/n + 10;

#param a{i in 1..n} := if ( (i mod 2) == 1 ) then (i*2/n + 1) else (5 - i*3/n);

#var x{1..n} >= -0.8, <= 0.8, := 0.0;

#minimize f:
        #sum {i in 1..n-1 by 2} (a[i]/(1+x[i])) + sum {i in 2..n by 2} (a[i]/(1-x[i]));
#subject to cons1{i in 6..n-4 by 2}:
        #1/(1-x[i-4]) + 1/(1+x[i-3]) + 1/(1+x[i-2]) + 1/(1-x[i-1]) + 1/(1+x[i]) + 1/(1+x[i+1]) + 1/(1-x[i+2]) + 1/(1+x[i+3]) + 1/(1-x[i+4]) <= b[i];
#subject to cons2:
        #1/(1-x[1]) + 1/(1-x[2]) + 1/(1+x[3]) + 1/(1-x[4]) + 1/(1+x[5]) + 1/(1+x[n-3]) + 1/(1-x[n-2]) + 1/(1-x[n-1]) + 1/(1+x[n]) <= b[1];
#subject to cons3:
        #1/(1-x[1]) + 1/(1+x[2]) + 1/(1+x[3]) + 1/(1-x[4]) + 1/(1+x[5]) + 1/(1-x[6]) + 1/(1-x[n-2]) + 1/(1+x[n-1]) + 1/(1+x[n]) <= b[2];
#subject to cons4:
        #1/(1-x[1]) + 1/(1+x[2]) + 1/(1-x[3]) + 1/(1-x[4]) + 1/(1+x[5]) + 1/(1-x[6]) + 1/(1+x[7]) + 1/(1+x[n-1]) + 1/(1-x[n]) <= b[3];
#subject to cons5:
        #1/(1+x[1]) + 1/(1+x[2]) + 1/(1-x[3]) + 1/(1+x[4]) + 1/(1+x[5]) + 1/(1-x[6]) + 1/(1+x[7]) + 1/(1-x[8]) + 1/(1-x[n]) <= b[4];
#subject to cons6:
        #1/(1+x[1]) + 1/(1+x[n-7]) + 1/(1-x[n-6]) + 1/(1-x[n-5]) + 1/(1+x[n-4]) + 1/(1-x[n-3]) + 1/(1-x[n-2]) + 1/(1+x[n-1]) + 1/(1-x[n]) <= b[n-3];
#subject to cons7:
        #1/(1+x[1]) + 1/(1-x[2]) +  1/(1-x[n-6]) + 1/(1+x[n-5]) + 1/(1+x[n-4]) + 1/(1-x[n-3]) + 1/(1+x[n-2]) + 1/(1+x[n-1]) + 1/(1-x[n]) <= b[n-2];
#subject to cons8:
        #1/(1+x[1]) + 1/(1-x[2]) +  1/(1+x[3]) + 1/(1+x[n-5]) + 1/(1-x[n-4]) + 1/(1-x[n-3]) + 1/(1+x[n-2]) + 1/(1-x[n-1]) + 1/(1-x[n]) <= b[n-1];
#subject to cons9:
        #1/(1+x[1]) + 1/(1-x[2]) +  1/(1+x[3]) + 1/(1-x[4]) + 1/(1-x[n-4]) + 1/(1+x[n-3]) + 1/(1+x[n-2]) + 1/(1-x[n-1]) + 1/(1+x[n]) <= b[n];
#subject to cons10{i in 5..n-5 by 2}:
        #1/(1+x[i-4]) + 1/(1-x[i+3]) + 1/(1-x[i-2]) + 1/(1+x[i-1]) + 1/(1-x[i]) + 1/(1-x[i+1]) + 1/(1+x[i+2]) + 1/(1-x[i+3]) + 1/(1+x[i+4]) <= b[i];

