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

def addExtra(metadata):

  x=SX("x")
  f = SXFunction([x],[x**2])
  f.init()
  i = IpoptSolver(f)

  for name in i.getOptionNames():
    if name in metadata["CasADi::IpoptInternal"]["options"]:
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

  
  x=SX("x")
  f = SXFunction([x],[x**2])
  f.init()
  try:
    i = WorhpSolver(f)
  except:
    return
    
  for name in i.getOptionNames():
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
