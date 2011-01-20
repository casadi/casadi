# -*- coding: utf-8 -*-
from casadi import *
from numpy import *
# This example has been ported to a unittesting framework

# variables
x = SX("x")
y = SX("y")

# Create function
f = [3-sin(x*x)-y, sqrt(y)*x]
print "f = ", f
fcn = SXFunction([[x,y]],[f])

# Set some options
fcn.setOption("name","f")
fcn.setOption("ad_order",1)

# Print
print "fcn = ", fcn

# Initialize the function for numerical calculations
fcn.init()

# Pass inputs
fcn.setInput([2,3])

# Pass forward seed
fcn.setFwdSeed([0,1]);

# Pass adjoint seed
fcn.setAdjSeed([0,1])

# Evaluate numerically
fcn.evaluate(1,1)

# Get the results
res = fcn.output()
print "result = ", res

fsens = fcn.fwdSens()
print "forward derivative = ", fsens

asens = fcn.adjSens()
print "adjoint derivative = ", asens

# evaluate symbolically
fsubst = fcn.eval([[1-y,1-x]])
print "fsubst = ", fsubst

#p = FMIParser()
#p.loadXML("/home/janderss/src/jmodelica/branches/1.3.x/Python/src/jmodelica/examples/VDP_pack_VDP_Opt.xml")
#print p
#ocp = p.parseOptimica()

