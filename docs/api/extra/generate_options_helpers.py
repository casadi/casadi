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

def optionDocumented(name,cl,metadata):
  acl = []
  if 'hierarchy' in metadata[cl]:
    acl += metadata[cl]['hierarchy']
  for c in acl:
    if 'options' in metadata[c] and name in metadata[c]['options']: return True
  return False
  
def extra(metadata,i,iname):
  print "Adding to ", metadata
  for name in i.optionNames():
    print "found option", name
    if optionDocumented(name,"casadi::%s" % iname,metadata):
      continue
    print "Adding it."
    meta = metadata["casadi::%s" % iname]["options"][name] = dict()
    meta['name'] = name
    meta['type'] = i.optionTypeName(name)
    meta['used'] = "casadi::%s" % iname
    meta['inherit'] = False
    meta['description'] = i.optionDescription(name)
    try:
      meta['default'] = i.optionDefault(name)
    except:
      meta['default'] = ''
      pass #too bad
    print meta

def addExtra(metadata):
  # This is not possible anymore
  return
  
  print "Adding extra"

  x=SX.sym("x")
  f = {'x':x, 'f':x**2}
  
  i = nlpsol("mysolver", "ipopt", f)
  extra(metadata,i,"IpoptInterface")

  x=SX.sym("x")
  i = nlpsol("mysolver", "worhp", f)
  extra(metadata,i,"WorhpInterface")

  try:
    x=SX.sym("x")
    i = nlpsol("mysolver", "snopt", f)
    extra(metadata,i,"SnoptInterface")
  except:
    pass
  
  i = qpsol("mysolver", "qpoases", {"h": Sparsity.dense(3,3),"a":Sparsity.dense(1,3)})
  extra(metadata,i,"QpoasesInterface")
 
  G = sparsify(DM([[1,0],[0,1]])).T
  E = sparsify(DM([0,0]))

  A = DM(0,2)

