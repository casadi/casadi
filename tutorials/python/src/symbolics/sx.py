#! CasADi tutorial 1
#! ==================
#! This tutorial file explains the use of CasADi's SX in a python context.
#! Let's start with the import statements to load CasADi.
from casadi import *
from numpy import *
#! Contructors & printing
#! --------------------------------------
#! Always provide a name to your symbols.
#! This name is not used for identification, only for printing.
a = SX("z")
print type(a)
print a
#! You can explicitely create constant SX objects as follows
c = SX(5)
print c
#! Scalar algebra
#! --------------------------
#! Any operations on SX objects return SX objects.
x = SX("x")
y = SX("y")
c = x+y
print type(c)
print c
#! While you construct ever complex expressions, 
#! a graph of SX objects is created.
d = c*2 + x
print d
#! Note that, by itself, CasADi does very little expansions or simplifications of expressions.
#! Only simplifications that do not to introduce new nodes to the graph are allowed.
print d-x
print simplify(d-x)
print SX(5) + SX(7)
print 0*x + 0*y
print 1*x
#! SX objects are immutable entities.
#! The assignment and update operators are really creating new object instead of modifying the old ones.
print "object address before: %d" % id(d)
d = d - x
print d
print "object address after: %d" % id(d)
#! Consequently, updates and assignements don't have side effects for other SX objects
f = x + y
x *= 2
print x
print f
print x+y
#! Expression substitution
#! ------------------------------------
x=symbolic("x")

y=x*x
print y
substitute(y,x,symbolic("w"))
print y
#! More operators
#! ------------------------------
#! Some familiar mathematical operations are supported that mimic the standard numpy functions:
#! sqrt sin cos tan arctan arcsin arccos exp log pow fabs floor ceil erf fmin fmax.
#! Using these functions require numpy to be imported.
y = sin(x**x)
print type(y>0)
t = if_else(y>0,-10,10)
print t

#! Conclusion
#! -------------------
#! We have seen how SX objects behave like symbolic objects.
#! They can be used to contruct expression trees.
#! 
#! To see how we can efficiently evaluate and differentiate these objects, jump on to the sxfunction tutorial...
