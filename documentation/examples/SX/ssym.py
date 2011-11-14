#! ssym
#!======================
from casadi import *

#! Construct using a single name
#! =====================================
#! The names of the entries of the SXMatrix will be derived from the name provided as argument to ssym.

#! Without shape arguments, a 1-by-1 matrix is constructed:

x = ssym("x")
print type(x), x

#! Create a column matrix 
print ssym("x",2,1)

#! Create a row matrix 
print ssym("x",1,2)

#! Create a matrix 
print ssym("x",2,3)

#! Construct using multiple names
#! =====================================

#! casADi defaults to a column matrix when no shape arguments are provided
print ssym("[a,b,c]")

#! Create a row matrix 
print ssym("[a,b,c]",1,3)

#! Space can be used as a separator, too:
print ssym("[a b c]")

#! Other brackets work fine as well:
print ssym("{a b c}")
print ssym("(a b c)")

#! Or you can omit the brackets entirely:
print ssym("a b c")

#! Create a matrix 
print ssym("[a b c;d e f]",2,3) 

#! The pythonic way to create a bunch of SX'es
#! ===========================================
#!
#! Thanks to tuple unpacking, you can write:
a,b,c = ssym("[a,b,c]")

#! This assigns three variables a, b and c:
print type(a), a
print type(b), b
print type(c), c

