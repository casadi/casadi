#! CasADi
#! ==================
#! More about SX
from casadi import *
from numpy import *
#! The identity of an SX node is very persistant.
#! We demonstrate this with the help of symbolic substitution.
x=SX("x")
y=x**2
print SXFunction([[x]],[[y]]).eval([SX("w")])
#! We expect w^2.
l = x
print SXFunction([[l]],[[y]]).eval([SX("w")])
#! We expect w^2.
k=SXMatrix(x)
l=k[0]
print SXFunction([[l]],[[y]]).eval([SX("w")])
#! We expect w^2.
k=symbolic("d",2,2)
k[1] = x
l=k[1]
print SXFunction([[l]],[[y]]).eval([SX("w")])
#! We expect w^2.
#! Identity is not associated with name:
l=SX("x")
print SXFunction([[l]],[[y]]).eval([SX("w")])
