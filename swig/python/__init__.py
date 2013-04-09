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

from casadi import *
import casadi

import types
  
def wrapper(f, warning,error=False):
    def new(*args, **kwargs):
      print "*" * 40
      print "Deprecation Warning"
      print "-" * 40
      print warning
      print "*" * 40
      if error:
        raise Exception("Deprecation error: " + warning)
      return f(*args, **kwargs)
    return new


class Deprecate(object):
    def __new__(self, o, warning,error=False):
        class temp(o.__class__): pass
        temp.__name__ = "Deprecated_%s" % o.__class__.__name__
        output = temp.__new__(temp, o)

        output.warned = True
        wrappable_types = (type(int.__add__), type(zip), types.FunctionType)
        unwrappable_names = ("__str__", "__unicode__", "__repr__", "__getattribute__", "__setattr__")

        for method_name in dir(temp):
            if not type(getattr(temp, method_name)) in wrappable_types: continue
            if method_name in unwrappable_names: continue

            setattr(temp, method_name, wrapper(getattr(temp, method_name), warning,error=error))

        output.warned = False
        return output

deprecated = {}
for oldenum, newenum, string in [
  ("NLP_X_INIT","NLP_SOLVER_X0","x0"),
  ("NLP_LBX","NLP_SOLVER_LBX","lbx"),
  ("NLP_UBX","NLP_SOLVER_UBX","ubx"),
  ("NLP_LBG","NLP_SOLVER_LBG","lbg"),
  ("NLP_UBG","NLP_SOLVER_UBG","ubg"),
  ("NLP_LAMBDA_INIT","NLP_SOLVER_LAMB_G0","lam_g0"),
  ("NLP_P","NLP_SOLVER_P","p"),
  ("NLP_X_OPT","NLP_SOLVER_X","x"),
  ("NLP_COST","NLP_SOLVER_F","f"),
  ("NLP_LAMBDA_G","NLP_SOLVER_LAM_G","lam_g"),
  ("NLP_LAMBDA_X","NLP_SOLVER_LAM_X","lam_x"),
  ("NLP_G","NLP_SOLVER_G","g"),
  ]:
    deprecated[oldenum] = Deprecate(0,"""%s is dropped. Use %s instead.\nYou can in fact avoid using either of them alltogether:\n  solver.input(%s)\n     becomes\n  solver.input("%s")\nreference: https://github.com/casadi/casadi/issues/566""" % (oldenum,newenum,oldenum,string),error=True)
          
for oldenum, newenum in [
    ("NLP_NUM_IN","NLP_SOLVER_NUM_IN"),
    ("NLP_NUM_OUT","NLP_SOLVER_NUM_OUT"),
  ]:
    deprecated[oldenum] = Deprecate(0,"""%s is dropped. Use %s instead.\nreference: https://github.com/casadi/casadi/issues/566""" % (oldenum,newenum),error=True)

NLP_X_INIT = deprecated["NLP_X_INIT"]
NLP_LBX = deprecated["NLP_LBX"]
NLP_UBX = deprecated["NLP_UBX"]
NLP_LBG = deprecated["NLP_LBG"]
NLP_UBG = deprecated["NLP_UBG"]
NLP_LAMBDA_INIT = deprecated["NLP_LAMBDA_INIT"]
NLP_P = deprecated["NLP_P"]
NLP_X_OPT = deprecated["NLP_X_OPT"]
NLP_COST = deprecated["NLP_COST"]
NLP_LAMBDA_G = deprecated["NLP_LAMBDA_G"]
NLP_LAMBDA_X = deprecated["NLP_LAMBDA_X"]
NLP_G = deprecated["NLP_G"]
NLP_NUM_IN = deprecated["NLP_NUM_IN"]
NLP_NUM_OUT = deprecated["NLP_NUM_OUT"]

__version__ = CasadiMeta.version
