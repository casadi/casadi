from numpy import *
from casadi import *
 
# create a simple function with sparse jacobian
x = symbolic("x",10)
f = x+5
f[2] = 44
f[3] = 55
f[4:7] += inner_prod(x[1:4],x[1:4])
fcn = SXFunction([x],[f])

# Use forward or adjoint mode ad to calculate directional derivatives (uncomment one to force a mode)
#fcn.setOption("ad_mode","forward")
#fcn.setOption("ad_mode","reverse")

fcn.init()

# some point to evaluate the jacobian
x0 = linspace(-3.,12.,10)

# input index (i.e. with respect to which (matrix-valued) variable are we differentiating
iind = 0

# output index (i.e. with (matrix-valued) function output are we differentiating
oind = 0

# create the jacobian in the convensional way (i.e. using source code transformation)
J1 = fcn.jacobian(iind,oind)
J1.init()
J1.setInput(x0)
J1.evaluate()
print "J1.output().size() = ", J1.output().size()
print "J1.output().numel() = ", J1.output().numel()
print "J1(x0)", array(J1.output())

# create the jacobian using compression techniques
v2 = VectorPair_Int_Int(1,(oind,iind)) # vector with a set of Jacobian blocks that we request
J2 = Jacobian(fcn,v2) # create a "Jacobian" function instance explicitly
J2.setOption("verbose",True) # so that it prints the number of directions
J2.init()
J2.setInput(x0)
J2.evaluate()
print "J2.output().size() = ", J2.output().size()
print "J2.output().numel() = ", J2.output().numel()
print "J2(x0)", array(J2.output())

# Print difference
print "Difference: ", J2.output()-J1.output()


