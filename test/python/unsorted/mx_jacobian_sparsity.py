from casadi import *

x = msym("x",5)
y = msym("y")
x2 = ssym("x",5)
y2 = ssym("y")

fcn = MXFunction([x,y],[4*vertcat((x[2:5],x[0:2])) + y*x])
fcn.init()
js = IMatrix(fcn.jacSparsityOld(),1)
js.printDense()

fcn2 = SXFunction([x2,y2],[4*vertcat((x2[2:5],x2[0:2])) + y2*x2])
fcn2.init()
js2 = IMatrix(fcn2.jacSparsityOld(),1)
js2.printDense()

fcn3 = MXFunction([x,y],fcn2.call([x,y]))
fcn3.init()
js3 = IMatrix(fcn3.jacSparsityOld(),1)
js3.printDense()
