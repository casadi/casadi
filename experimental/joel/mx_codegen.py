#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             K.U. Leuven. All rights reserved.
#     Copyright (C) 2011-2014 Greg Horn
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
from numpy import *
from os import system
import time
import sys

x = ssym("x",3,4)
f1 = SXFunction([x],[sin(x)])
f1.init()

a = msym("A",3,4)
b = msym("b",4,1)
[z] = f1.call([a])
c = mul(z,b)
c = 2 + c
f = MXFunction([a,b],[c])
f.init()
f.generateCode("f_mx.c")
system("gcc -fPIC -shared f_mx.c  -o f_mx.so")

ef = ExternalFunction("./f_mx.so")
ef.init()

a_val = array([[1,2,3,3],[2,3,4,5],[3,4,5,6]])
b_val = array([7,6,5,4])

f.setInput(a_val,0);
f.setInput(b_val,1);
f.evaluate()

print f.getOutput()

ef.setInput(a_val,0);
ef.setInput(b_val,1);
ef.evaluate()

print f.getOutput()

