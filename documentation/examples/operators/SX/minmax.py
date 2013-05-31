#
#     This file is part of CasADi.
# 
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
# 
#     CasADi is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 3 of the License, or (at your option) any later version.
# 
#     CasADi is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
# 
#     You should have received a copy of the GNU Lesser General Public
#     License along with CasADi; if not, write to the Free Software
#     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
# 
# 
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
f.setInput(4,0)
f.setInput(6,1)
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
J.setInput(4,0)
J.setInput(6,1)
J.evaluate()

#! fmax is not sensitive to the first argument (smallest)
assert J.output(0)[0]==0

#! fmax is not sensitive to the first argument (biggest)
J.setInput(6,0)
J.setInput(4,1)
J.evaluate()
assert J.output(0)[0]==1


#! We get a heaviside when taking the derivative of fmin
j = jacobian(min_,x)
print j

J=SXFunction([x,y],[j])
J.init()
J.setInput(4,0)
J.setInput(6,1)
J.evaluate()

#! fmax is not sensitive to the first argument (smallest)
assert J.output(0)[0]==1

#! fmax is not sensitive to the first argument (biggest)
J.setInput(6,0)
J.setInput(4,1)
J.evaluate()
assert J.output(0)[0]==0


#! AD forward on fmin, fmax
#! ------------------------

f.setFwdSeed(1,0)
f.setFwdSeed(0,1)
f.evaluate(1,0)

#! fmax is not sensitive to the first argument (smallest)
assert f.fwdSens(0)[0]==0

#! fmin is only sensitive to the first argument (smallest)
assert f.fwdSens(1)[0]==1

f.setFwdSeed(0,0)
f.setFwdSeed(1,1)
f.evaluate(1,0)

#! fmax is only sensitive to the second argument (largest)
assert f.fwdSens(0)[0]==1

#! fmin is not sensitive to the second argument (largest)
assert f.fwdSens(1)[0]==0

#! AD adjoint on fmin, fmax
#! ------------------------

f.setAdjSeed(1,0)
f.setAdjSeed(0,1)
f.evaluate(0,1)

#! The first argument (smallest) is not influenced by fmax
assert f.adjSens(0)[0]==0

#! The first argument (smallest) is only influenced by fmin
assert f.adjSens(1)[0]==1

f.setAdjSeed(0,0)
f.setAdjSeed(1,1)
f.evaluate(0,1)

#! The second argument (largest) is only influenced by fmax
assert f.adjSens(0)[0]==1

#! The first argument (largest) is not influenced by fmin
assert f.adjSens(1)[0]==0



#! On the borderline
#! -----------------
#! How do the sensitivities behave when both arguments are the same?
f.setInput(5,0)
f.setInput(5,1)

f.setFwdSeed(1,0)
f.setFwdSeed(0,1)
f.evaluate(1,0)

#! fmax sensitivity to the first argument:
print f.fwdSens(0)[0]

#! fmin sensitivity to the first argument:
print f.fwdSens(1)[0]

f.setFwdSeed(0,0)
f.setFwdSeed(1,1)
f.evaluate(1,0)

#! fmax sensitivity to the second argument:
print f.fwdSens(0)[0]

#! fmin sensitivity to the second argument:
print f.fwdSens(1)[0]

f.setAdjSeed(1,0)
f.setAdjSeed(0,1)
f.evaluate(0,1)

#! first argument sensitivity to fmax:
print f.adjSens(0)[0]

#! first argument sensitivity to fmin:
print f.adjSens(1)[0]

f.setAdjSeed(0,0)
f.setAdjSeed(1,1)
f.evaluate(0,1)

#! second argument sensitivity to fmax:
print f.adjSens(0)[0]

#! second argument sensitivity to fmin:
print f.adjSens(1)[0]

