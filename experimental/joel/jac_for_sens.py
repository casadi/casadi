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
a = ssym("a",5)
b = ssym("b",5)

for jac_for_sens in (True,False):
  print "jac_for_sens = ", jac_for_sens, ":"
  f = SXFunction([a,b],[sqrt(b-sin(a)),inner_prod(a,b),outer_prod(a,b)])
  f.setOption("jac_for_sens",jac_for_sens)
  f.init()
  #print f.jacSparsity(0,0)

  f.setInput([1,2,3,4,5],0)
  f.setInput([10,20,30,40,50],1)

  f.setFwdSeed([1,0,0,0,1],0)
  f.setFwdSeed([0,0,0,0,0],1)

  f.setAdjSeed([0,0,0,0,0],0)
  f.setAdjSeed(1,1)
  f.setAdjSeed(DMatrix(5,5,0),2)

  f.evaluate(1,1)

  print "f.getFwdSens(0) = ", f.fwdSens(0).data()
  print "f.getFwdSens(1) = ", f.fwdSens(1).data()
  print "f.getFwdSens(2) = ", f.fwdSens(2).data()
  print "f.getAdjSens(0) = ", f.adjSens(0).data()
  print "f.getAdjSens(1) = ", f.adjSens(1).data()

  print "--"