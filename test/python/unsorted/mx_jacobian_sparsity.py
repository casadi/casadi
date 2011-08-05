from casadi import *

x = MX("x",5)
y = MX("y")
x2 = symbolic("x",5)
y2 = symbolic("y")

fcn = MXFunction([x,y],[4*vertcat((x[2:5],x[0:2])) + y*x])
fcn.init()
js = IMatrix(fcn.jacSparsity(),1)
js.printDense()

fcn2 = SXFunction([x2,y2],[4*vertcat((x2[2:5],x2[0:2])) + y2*x2])
fcn2.init()
js2 = IMatrix(fcn2.jacSparsity(),1)
js2.printDense()

