from casadi import *
a = ssym("a",2,2)
q,r = qr(a)
printCompact(q)
printCompact(r)


