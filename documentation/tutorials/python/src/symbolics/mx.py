#! CasADi tutorial 3
#! ==================
#! This tutorial file explains the use of CasADi's MX in a python context.
from numpy import *
from casadi import *
#! Let's start with some algebra
x = MX("x",2,3)
y = MX("y",3,2)
print x
for i in range(6):
	print x[i]

for i in range(2):
	for j in range(3):
		print "x[%d,%d] = %s" % (i,j,str(x[i,j]))
		
print x[1,1]
print x[3] # Note that index is flattened. x[0,0] is illegal.
print norm_2(x)
z= prod(x,y)
print z
#! Note how the operations on MXes are lazy on the matrix level.
#! Any elementwise logic is postponed until evaluation demands it.
#! Just like, SXFunction, MXFunction can be single or multi input/output.
#! The only allowed input/output primitive is MX.
f = MXFunction([x,y],[z])


#! Evaluation
#! ----------------
f.init()
f.setInput([1,2,3,4,5,6],0);
f.setInput([1,3,0,6,0,9],1);
f.evaluate()

print f.output()
#! Note how this result is related to a numpy approach:
a=matrix(f.input(0)).reshape(3,2)
b=matrix(f.input(1)).reshape(2,3)
print a.T*b.T
#! Jacobian
#! -------------
#! BUG error
#!print f.jacobian(0)
#!type(f.hessian())

#! Numerical matrices
#! ------------------------------
X = MX(DMatrix(array([[1,2,3],[4,5,6]])))
print X
print outer_prod(X,X)
print MX(trans(DMatrix([1,2,3])))
print MX([1,2,3])
#! As before, evaluation is lazy on the matrix level
Y = MX("Y")
f = MXFunction([Y],[X])
f.init()
f.input().set([2])
f.evaluate()
print f.output()

#! Element assignement
#! -------------------
X = MX("x",2,2)
X[0,0]=MX(5)
print X
 
