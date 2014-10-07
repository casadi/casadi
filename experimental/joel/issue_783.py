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

x = msym("x",2,1)

fun = MXFunction([x],[diag(x[[1,0]])])
fun.init()

for f,sym,Function,X in [(fun,msym,MXFunction,MX),]:
  f.init()
  print Function
  
  a = sym("a",sp_diag(2))

  _,_,[[f1]] = f.eval([X.ones(2)],[],[[a]])

  vf = Function([a],[f1])
  vf.init()
  print vf

  vf.setInput(1.0)
  vf.setAdjSeed([0,1])
  vf.evaluate(0,1)
  vf.getAdjSens().printDense()


  a2 = sym("a2",2)
  _,_,[[f2]] = vf.eval([X.ones(sp_diag(2))],[],[[a2]])

  vf2 = Function([a2],[f2])
  vf2.init()

  vf2.setInput([0,1])
  vf2.evaluate()
  vf2.getOutput().printDense()

