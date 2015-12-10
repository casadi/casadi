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
# -*- coding: utf-8 -*-

from graph import *
from bounds import *
from structure import repeated, entry, struct_symSX, struct_symMX, struct_SX, struct_MX, struct_MX_mutable, nesteddict, index, indexf, struct, struct_load
from io import nice_stdout, capture_stdout

def print_subclasses(myclass, depth=0):
  print ("  " * depth) + " - " + myclass.__name__
  for s in myclass.__subclasses__():
    print_subclasses(s,depth=depth+1)
    
def loadAllCompiledPlugins():
  for k in CasadiMeta.getPlugins().split(";"):
    cls, name = k.split("::")
    print "Testing: ", cls, name
    if cls=='Integrator':
      casadi.load_integrator(name)
    elif cls=='Nlpsol':
      casadi.load_nlpsol(name)
    elif cls=='Qpsol':
      casadi.load_qpsol(name)
    elif cls=='Rootfinder':
      casadi.load_rootfinder(name)
    elif cls=='Linsol':
      casadi.load_linsol(name)
    else:
      getattr(casadi,cls).loadPlugin(name)
