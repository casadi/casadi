from casadi import *

x = MX.sym("x",3)


f = MXFunction("f",[x],[sumRows(sin(x))*cos(x[0])])
f.init()

h = f.hessian()

x0 = DMatrix([3.7,2,3.2])

print h([x0])[0]


print "Forward over adjoint"

fad = f.derivative(0,1)

ffwd = fad.derivative(3,0)

print horzcat(ffwd([x0,1,DMatrix([1,0,0]),0,DMatrix([0,1,0]),0,DMatrix([0,0,1]),0])[3::2]).T


print "Forward over forward"

fwd = f.derivative(3,0)

ffwd = fwd.derivative(3,0)

ins = [x0,DMatrix([1,0,0]),DMatrix([0,1,0]),DMatrix([0,0,1])]
print vertcat(ffwd(ins+[DMatrix([1,0,0]),0,0,0]+[DMatrix([0,1,0]),0,0,0]+[DMatrix([0,0,1]),0,0,0])[4:]).reshape((4,3))[1:,:]

print "Adjoint over adjoint"

fad = f.derivative(0,1)

ffad = fad.derivative(0,3)

print horzcat(ffad([x0,1,0,DMatrix([1,0,0]),0,DMatrix([0,1,0]),0,DMatrix([0,0,1])])[2::2]).T

print "Adjoint over forward"
fwd = f.derivative(3,0)
ffad = fwd.derivative(3,0)

ins = [x0,DMatrix([1,0,0]),DMatrix([0,1,0]),DMatrix([0,0,1])]
print vertcat(ffad(ins+[DMatrix([1,0,0]),0,0,0]+[DMatrix([0,1,0]),0,0,0]+[DMatrix([0,0,1]),0,0,0])[4:]).reshape((4,3))[1:,:]


