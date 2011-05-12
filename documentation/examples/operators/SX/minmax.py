#! Symbolic substitution
#!======================
from casadi import *

x=SX("x")
y=SX("y")

max_ = x.fmax(y)
min_ = x.fmin(y)

print max_, min_

#! Let's construct an SXFunction with max and min as outputs
f = SXFunction([x,y],[max_,min_])
f.init()
#! We evaluate for x=4, y=6
f.input(0).set(4)
f.input(1).set(6)
f.evaluate()

#! max(4,6)=6
assert f.output(0)[0]==6
print f.output(0)[0]

#! min(4,6)=4
assert f.output(1)[0]==4
print f.output(1)[0]

#! Symbolic differentiation
#! ------------------------

#! We get a heaviside when taking the derivative of fmax
j = jacobian(max_,x)
print j


J=SXFunction([x,y],[j])
J.init()
J.input(0).set(4)
J.input(1).set(6)
J.evaluate()

#! fmax is not sensitive to the first argument (smallest)
assert J.output(0)[0]==0

#! fmax is not sensitive to the first argument (biggest)
J.input(0).set(6)
J.input(1).set(4)
J.evaluate()
assert J.output(0)[0]==1


#! We get a heaviside when taking the derivative of fmin
j = jacobian(min_,x)
print j

J=SXFunction([x,y],[j])
J.init()
J.input(0).set(4)
J.input(1).set(6)
J.evaluate()

#! fmax is not sensitive to the first argument (smallest)
assert J.output(0)[0]==1

#! fmax is not sensitive to the first argument (biggest)
J.input(0).set(6)
J.input(1).set(4)
J.evaluate()
assert J.output(0)[0]==0


#! AD forward on fmin, fmax
#! ------------------------

f.fwdSeed(0).set(1)
f.fwdSeed(1).set(0)
f.evaluate(1,0)

#! fmax is not sensitive to the first argument (smallest)
assert f.fwdSens(0)[0]==0

#! fmin is only sensitive to the first argument (smallest)
assert f.fwdSens(1)[0]==1

f.fwdSeed(0).set(0)
f.fwdSeed(1).set(1)
f.evaluate(1,0)

#! fmax is only sensitive to the second argument (largest)
assert f.fwdSens(0)[0]==1

#! fmin is not sensitive to the second argument (largest)
assert f.fwdSens(1)[0]==0

#! AD adjoint on fmin, fmax
#! ------------------------

f.adjSeed(0).set(1)
f.adjSeed(1).set(0)
f.evaluate(0,1)

#! The first argument (smallest) is not influenced by fmax
assert f.adjSens(0)[0]==0

#! The first argument (smallest) is only influenced by fmin
assert f.adjSens(1)[0]==1

f.adjSeed(0).set(0)
f.adjSeed(1).set(1)
f.evaluate(0,1)

#! The second argument (largest) is only influenced by fmax
assert f.adjSens(0)[0]==1

#! The first argument (largest) is not influenced by fmin
assert f.adjSens(1)[0]==0



#! On the borderline
#! -----------------
#! How do the sensitivities behave when both arguments are the same?
f.input(0).set(5)
f.input(1).set(5)

f.fwdSeed(0).set(1)
f.fwdSeed(1).set(0)
f.evaluate(1,0)

#! fmax sensitivity to the first argument:
print f.fwdSens(0)[0]

#! fmin sensitivity to the first argument:
print f.fwdSens(1)[0]

f.fwdSeed(0).set(0)
f.fwdSeed(1).set(1)
f.evaluate(1,0)

#! fmax sensitivity to the second argument:
print f.fwdSens(0)[0]

#! fmin sensitivity to the second argument:
print f.fwdSens(1)[0]

f.adjSeed(0).set(1)
f.adjSeed(1).set(0)
f.evaluate(0,1)

#! first argument sensitivity to fmax:
print f.adjSens(0)[0]

#! first argument sensitivity to fmin:
print f.adjSens(1)[0]

f.adjSeed(0).set(0)
f.adjSeed(1).set(1)
f.evaluate(0,1)

#! second argument sensitivity to fmax:
print f.adjSens(0)[0]

#! second argument sensitivity to fmin:
print f.adjSens(1)[0]

