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

def color(A):
  print "="*80
  print "Original:"
  A.spy()
  print "Colored: "
  A.unidirectionalColoring().spy()

A = sp_diag(5)
color(A)
print "One direction needed to capture all"
color(sp_dense(5,10))
print "We need 5 directions."
print "The colored response reads: each row corresponds to a direction;"
print " each column correspond to a row of the original matrix."

color(A+sp_triplet(5,5,[0],[4]))
print "First 4 rows can be taken together, the fifth row is taken seperately"
color(A+sp_triplet(5,5,[4],[0]))
print "First 4 rows can be taken together, the fifth row is taken seperately"

color(A+sp_triplet(5,5,[0]*5,range(5)))
print "The first row is taken seperately."
print "The remainding rows are lumped together in one direction."

color(A+sp_triplet(5,5,range(5),[0]*5))
print "We need 5 directions."
