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
from numpy import Inf

def reportBounds(value,lowerbound,upperbound,labels=None,tol=1e-8,showNonViolating=True):
  if hasattr(labels,"labels"):
    labels = labels.labels()

  v = list(value.nonzeros())
  lb = list(lowerbound.nonzeros())
  ub = list(upperbound.nonzeros())

  if labels is None:
    labels = [""] * len(v) 
  
  if not(len(v)==len(lb) and len(lb)==len(ub)):
    raise Exception("value, lowerbound and upperbound must all be the same size, but got %d, %d and %d. " % (len(v),len(lb),len(ub)))
  
  if len(labels)!=len(v):
    raise Exception("Labels (%d) must be same size as values (%d)" % (len(labels),len(v)))
  
  if ( all(value <= upperbound + tol) and all(value >= lowerbound - tol) ):
    print("All %d bounds are met: " % value.size())
  else:
    print("Problem with bounds : ")
  
  print("-"*60)

  # The length of the numeric fields
  fieldlength = 10
  # The length of the constraint visualizer strip
  indicator_length = 15
  

      
  # Loop over the elements of value
  for i in range(value.size()):
    violated =  (v[i] > (ub[i] + tol) or v[i] < (lb[i] - tol))
    nonregular = not is_regular([v[i]])

    if not(showNonViolating) and not(violated) and not(nonregular):
      continue
      
    if labels is None:
      identifier = "%d." % i
    else:
      identifier = labels[i]
         
   
      
    if ( abs(lb[i] - ub[i])<=tol):
      midfield = "%*s == %*s " % (fieldlength, "%.7e" % lb[i], fieldlength, "%.7e" % v[i])
      indicator = ""
    else:
      indicator = "-" * indicator_length
      if lb[i]==Inf:
        indicator = "8" +  indicator
      elif (abs(v[i]-lb[i])<=tol):
        indicator = "X" +  indicator
      else:
        indicator = "o" +  indicator
      if ub[i]==Inf:
        indicator += "8"
      elif (abs(v[i]-ub[i])<=tol):
        indicator += "X"
      else:
        indicator += "o"

      if (v[i] <= (ub[i] + tol) and v[i] >= (lb[i] - tol)):
        index = (v[i]-lb[i])/(ub[i]-lb[i])*(indicator_length-1)
        index = min(max(0,index),indicator_length-1)
        index = int(index)
        indicator = indicator[:1+index] + '=' + indicator[1+index:]
      
      midfield = "%*s <= %*s <= %*s" % (fieldlength, "%.7e" % lb[i], fieldlength, "%.7e" % v[i], fieldlength, "%.7e" % ub[i])

    if (v[i] > (ub[i] + tol) or v[i] < (lb[i] - tol)):
      indicator = " VIOLATED "
    if nonregular:
      indicator = " !REGULAR "
  
    print("%15s | %*s | %*s" % (identifier, (fieldlength + 6) * 3 , midfield, indicator_length+3, indicator))
  
  print("-" * 60)
  
