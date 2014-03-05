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

def optionDocumented(name,cl,metadata):
  acl = []
  if 'hierarchy' in metadata[cl]:
    acl += metadata[cl]['hierarchy']
  for c in acl:
    if 'options' in metadata[c] and name in metadata[c]['options']: return True
  return False

def addExtra(metadata):

  x=SX.sym("x")
  f = SXFunction(nlpIn(x=x),nlpOut(f=x**2))
  f.init()
  i = IpoptSolver(f)
  
  for name in i.getOptionNames():
    if optionDocumented(name,"CasADi::IpoptInternal",metadata):
      continue
    meta = metadata["CasADi::IpoptInternal"]["options"][name] = dict()
    meta['name'] = name
    meta['type'] = i.getOptionTypeName(name)
    meta['used'] = 'CasADi::IpoptInternal'
    meta['inherit'] = False
    meta['description'] = i.getOptionDescription(name)
    try:
      meta['default'] = i.getOptionDefault(name)
    except:
      meta['default'] = ''
      pass #too bad
    #if (len(i.getOptionAllowed(name))>1):
    #  meta['description'] += "(" + "|".join(i.getOptionAllowed(name))  + ")"

  
  x=SX.sym("x")
  f = SXFunction(nlpIn(x=x),nlpOut(f=x**2))
  f.init()
  try:
    i = WorhpSolver(f)
  except:
    return
    
  for name in i.getOptionNames():
    if optionDocumented(name,"CasADi::WorhpInternal",metadata): continue
    meta = metadata["CasADi::WorhpInternal"]["options"][name] = dict()
    meta['name'] = name
    meta['type'] = i.getOptionTypeName(name)
    meta['used'] = 'CasADi::WorhpInternal'
    meta['inherit'] = False
    meta['description'] = i.getOptionDescription(name)
    try:
      meta['default'] = i.getOptionDefault(name)
    except:
      meta['default'] = ''
      pass #too bad
    #if (len(i.getOptionAllowed(name))>1):
    #  meta['description'] += "(" + "|".join(i.getOptionAllowed(name))  + ")"

  try:
    i = QPOasesSolver(qpStruct(h=sp_dense(3,3),a=sp_dense(1,3)))
  except:
    return
    
  for name in i.getOptionNames():
    if optionDocumented(name,"CasADi::QPOasesInternal",metadata): continue
    meta = metadata["CasADi::QPOasesInternal"]["options"][name] = dict()
    meta['name'] = name
    meta['type'] = i.getOptionTypeName(name)
    meta['used'] = 'CasADi::QPOasesInternal'
    meta['inherit'] = False
    meta['description'] = i.getOptionDescription(name)
    try:
      meta['default'] = i.getOptionDefault(name)
    except:
      meta['default'] = ''
      pass #too bad
    #if (len(i.getOptionAllowed(name))>1):
    #  meta['description'] += "(" + "|".join(i.getOptionAllowed(name))  + ")"
