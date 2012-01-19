from casadi import *
from casadi.tools.variables import *
from numpy import Inf

def reportBounds(value,lowerbound,upperbound,labels=None,tol=1e-8,showNonViolating=True):
  if isinstance(labels,Variables):
    labels = labels.getLabels()
  
  v = list(value.data())
  lb = list(lowerbound.data())
  ub = list(upperbound.data())
  
  if not(len(v)==len(lb) and len(lb)==len(ub)):
    raise Exception("value, lowerbound and upperbound must all be the same size, but got %d, %d and %d. " % (len(v),len(lb),len(ub)))
  
  if (labels is None):
    if len(labels)!=len(v):
      raise Exception("Labels (%d) must be same size as values (%d)" % (len(labels),len(v)))
  
  if ( all(value <= upperbound + tol) and all(value >= lowerbound - tol) ):
    print "All %d bounds are met: " % value.size()
  else:
    print "Problem with bounds : "
  
  print "-"*60

  # The length of the numeric fields
  fieldlength = 10
  # The length of the constraint visualizer strip
  indicator_length = 15
  

      
  # Loop over the elements of value
  for i in range(value.size()):
    violated =  (v[i] > (ub[i] + tol) or v[i] < (lb[i] - tol))
    if not(showNonViolating) and not(violated):
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
  
    print "%15s | %*s | %*s" % (identifier, (fieldlength + 6) * 3 , midfield, indicator_length+3, indicator)
  
  print "-" * 60
  
