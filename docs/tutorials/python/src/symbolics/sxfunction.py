#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             K.U. Leuven. All rights reserved.
#     Copyright (C) 2011-2014 Greg Horn
#
#     CasADi is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 3 of the License, or (at your option) any later version.
#
#     CasADi is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public
#     License along with CasADi; if not, write to the Free Software
#     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
#
#! CasADi tutorial 2
#! ==================
#! This tutorial file explains the use of CasADi's SXFunction in a python context.
#! We assume you have read trough the SXElement tutorial.
#! 
#! Introduction
#! --------------------
#! Let's start with creating a simple expression tree z.
from casadi import *
from numpy import *
x = SX.sym("x")
y = x**2
z = sin(y) + y
print z
#! The printout value of z may give you the false impression that the evaluation of z will involve two multiplications of x.
#! This is not the case. This is what's going on under the hood:
#!
#! The expression tree of z does not contain two subexpressions x*x, rather it contains two pointers to a signle subexpression x*x.
#! In fact, in the C++ implementation, an SXElement object is really not more than a collection of pointers. It are SXNode objects which really contain the data associated with subexpressions.
#!
#! CasADi generates SXnodes at a very fine-grained level. Even 'sin(y)' is an SXNode, even though we have not ourselves declared a variable to point to it.
#! 
#! When evaluating z for a particular numerical value of x, the product is computed only once.
#!
#! Functions
#! --------------------------------
#! CasADi's SXFunction has powerful input/output behaviour.
#! The following input/output primitives are supported:
#$ \begin{description}
#$ \item[scalar] e.g. 'SX.sym("x")'
#$ \item[vector] e.g. 'SX.sym("x",3)'
#$ \item[matrix] e.g. 'SX.sym("x",3,2)'
#$ \end{description}
#! A function that uses one primitive as input/output is said to be 'single input'/'single output'.
#!
#! In general, an SXfunction can map from-and-to list/tuples of these primitives.
#!
#! Functions with scalar valued input
#! ----------------------------------
#! The following code creates and evaluates a single input (scalar valued), single output (scalar valued) function.
#$ f : $\mathbb{R} \mapsto \mathbb{R}$
f = SXFunction([x], [z]) # z = f(x)
print "%d -> %d" % (f.getNumInputs(),f.getNumOutputs())
print f.inputExpr(), type(f.inputExpr())
print f.outputExpr(), type(f.outputExpr())
f.init()
f.setInput(2)
f.evaluate()
print f.getOutput()
print type(f.getOutput())
#! Reevaluation does not require the init call.
f.setInput(3)
f.evaluate()
print f.getOutput()
#! We can evaluate symbolically, too:
print f([y])
#! Since numbers get cast to SXConstant object, you can also write the following non-efficient code:
print f([2])
#! We can do symbolic derivatives: f' = dz/dx . 
#$ The result is $2 x \cos(x^2)+2 x$:
print f.grad()
#! The following code creates and evaluates a multi input (scalar valued), multi output (scalar valued) function.
#$ The mathematical notation could be $ f_{i,j} : $\mathbb{R} \mapsto \mathbb{R} \quad  i,j \in [0,1]$
x = SX.sym("x") # 1 by 1 matrix serves as scalar
y = SX.sym("y") # 1 by 1 matrix serves as scalar
f = SXFunction([x , y ], [x*y, x+y])
print "%d -> %d" % (f.getNumInputs(),f.getNumOutputs())
f.init()
f.setInput(2,0)
f.setInput(3,1)
f.evaluate()

print [f.getOutput(i) for i in range(2)]
print [[f.grad(i,j) for i in range(2)] for j in range(2)]

#! Symbolic function manipulation
#! ------------------------------
#$ Consider the function $f(x;a;b) = a*x + b$
x=SX.sym("x")
a=SX.sym("a")
b=SX.sym("b")
f = SXFunction([x,vertcat((a,b))],[a*x + b]) 
f.init()

print f([x,vertcat([a,b])])
print f([SX(1.0),vertcat((a,b))])
print f([x,vertcat((SX.sym("c"),SX.sym("d")))])
print f([SX(),vertcat([SX.sym("c"),SX.sym("d")])])

#$ We can make an accompanying $g(x) = f(x;a;b)$ by making a and b implicity:

k = SX(a)
print f([x,vertcat((k[0],b))])
print f([x,vertcat((SX.sym("c"),SX.sym("d")))])

#! Functions with vector valued input
#! ----------------------------------
#! The following code creates and evaluates a single input (vector valued), single output (vector valued) function.
#$ f : $\mathbb{R}^2 \mapsto \mathbb{R}^2$

x = SX.sym("x")
y = SX.sym("y")
f = SXFunction([vertcat((x , y ))], [vertcat((x*y, x+y))])
print "%d -> %d" % (f.getNumInputs(),f.getNumOutputs())
f.init()
f.setInput([2,3])
f.evaluate()
print f.getOutput()
G=f.jac().T
print G

#$ Notice how G is a 2-nd order tensor $ {\buildrel\leftrightarrow\over{G}} = \vec{\nabla}{\vec{f}} = \frac{\partial [x*y, x+y]}{\partial [x , y]}$
#$ Let's define $ \vec{v} = {\buildrel\leftrightarrow\over{G}} . \vec{p} $
#! The evaluation of v can be efficiently achieved by automatic differentiation as follows:
f.init()
df = f.derivative(1,0)
df.setInput([2,3],0)
df.setInput([7,6],1)
df.evaluate()
print df.output(1) # v

#! Functions with matrix valued input
#! ----------------------------------
x = SX.sym("x",2,2)
y = SX.sym("y",2,2)
print x*y # Not a dot product
f = SXFunction([x,y], [x*y])
f.init()
print "%d -> %d" % (f.getNumInputs(),f.getNumOutputs())
print f([x,y])
f.init()
f.setInput([1,2,3,4],0); # instead of f.setInput([[1,2],[3,4]],0);
f.setInput([4,5,6,7],1);
f.evaluate()
print f.getOutput()
print f.jac(0).T
print f.jac(1).T

print 12

f = SXFunction([x,y], [x*y,x+y])
f.init()
print type(x)
print f([x,y])
print type(f([x,y]))
print type(f([x,y])[0])
print type(f([x,y])[0][0,0])


f = SXFunction([x], [x+y])
f.init()
print type(x)
print f([x])
print type(f([x]))
print type(f([x])[0])
print type(f([x])[0][0,0])

#! A current limitation is that matrix valued input/ouput is handled through flattened vectors
#! Note the peculiar form of the gradient.
#! 
#! Conclusion
#! ----------
#! This tutorial showed how SXFunction allows for symbolic or numeric evaluation and differentiation.
