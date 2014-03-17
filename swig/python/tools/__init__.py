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
# -*- coding: utf-8 -*-

from graph import *
from bounds import *
from structure import repeated, entry, struct_symSX, struct_symMX, struct_SX, struct_MX, struct_MX_mutable, nesteddict, index, indexf, struct, struct_load
from io import nice_stdout, capture_stdout

def print_subclasses(myclass, depth=0):
  print ("  " * depth) + " - " + myclass.__name__
  for s in myclass.__subclasses__():
    print_subclasses(s,depth=depth+1)
    
def deprecatedfun(func,message):

    def depfun(*args, **kwargs):
        if not func._fired:
          warnings.warn(message,SyntaxWarning,stacklevel=2)
          func._fired = True
        return func(*args, **kwargs)
    func._fired = False
    depfun.__dict__.update(func.__dict__)
    depfun.__name__ = func.__name__
    depfun.__doc__ = func.__doc__
    return depfun

struct_ssym = deprecatedfun(struct_symSX,"struct_ssym will soon be replaced by struct_symSX")
struct_msym = deprecatedfun(struct_symMX,"struct_msym will soon be replaced by struct_symMX")

