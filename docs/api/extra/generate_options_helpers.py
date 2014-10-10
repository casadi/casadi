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
  for name in i.getOptionNames():
    if optionDocumented(name,"casadi::%s" % iname,metadata):
      continue
    meta = metadata["casadi::%s" % iname]["options"][name] = dict()
    meta['name'] = name
    meta['type'] = i.getOptionTypeName(name)
    meta['used'] = "casadi::%s" % iname
    meta['inherit'] = False
    meta['description'] = i.getOptionDescription(name)
    try:
      meta['default'] = i.getOptionDefault(name)
    except:
      meta['default'] = ''
      pass #too bad

def addExtra(metadata):

  x=SX.sym("x")
  f = SXFunction(nlpIn(x=x),nlpOut(f=x**2))
  f.init()
  try:
    NlpSolver.loadPlugin("ipopt")
    i = NlpSolver("ipopt", f)
    extra(metadata,i,"IpoptInterface")
  except Exception as e:
    print e 
  
  x=SX.sym("x")
  f = SXFunction(nlpIn(x=x),nlpOut(f=x**2))
  f.init()
  try:
    NlpSolver.loadPlugin("worhp")
    i = NlpSolver("worhp", f)
    extra(metadata,i,"WorhpInterface")
  except Exception as e:
    print e
    
  x=SX.sym("x")
  f = SXFunction(nlpIn(x=x),nlpOut(f=x**2))
  f.init()
  try:
    NlpSolver.loadPlugin("snopt")
    i = NlpSolver("snopt", f)
    extra(metadata,i,"SnoptInterface")
  except Exception as e:
    print e
 
  try:
    i = QpSolver("qpoases", qpStruct(h=Sparsity.dense(3,3),a=Sparsity.dense(1,3)))
    extra(metadata,i,"QpoasesInterface")
  except Exception as e:
    print e
    

