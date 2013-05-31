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
from casadi import *



#! Note: the actual printing doesn't show up in this pdf, due to limitations of the reporting software.
#! Better run this example in a terminal.



a = SX("a")
b = SX("b")

c = a+b
c = c.printme(13)

d = c**2

print d

f = SXFunction([a,b],[d])
f.init()
f.setInput(4,0)
f.setInput(3,1)

#! When the graph is evaluated, a printout of c will occur (if you have set WITH_PRINTME to ON in CMakeCache.txt)
#! Printout reads '|> 13: 7'
#! 13 is an identifier of choice, 7 is the numerical value of c
f.evaluate()

J = f.jacobian(0,0)
J.init()

J.init()
J.setInput(2,0)
J.setInput(9,1)

#! The first derivative still depends on c
#! Printout reads '|> 13: 11'
J.evaluate()


J = J.jacobian(0,0)
J.init()

J.init()
J.setInput(2,0)
J.setInput(9,1)
#! second derivative doesn't, so we don't get a printout
J.evaluate()
