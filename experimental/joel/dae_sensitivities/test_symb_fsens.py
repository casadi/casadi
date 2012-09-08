from casadi import *
x = ssym("x",3)
f = SXFunction(daeIn(x),daeOut(x))
I = CVodesIntegrator(f)
I.init()
Id = I.derivative(1,0)