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

import warnings
import os

# workaround for issue #1012
# Since ipopt uses "dlopen" internally, we have to make sure that all the
# libraries we linked against are visible to the libraries dlopen loads.
# Specifically, hsl.so needs blas and lapack.
try:
    import numpy as np
except:
    pass

import sys
import ctypes

if hasattr(sys,"getdlopenflags"):
    flags0 = sys.getdlopenflags() # get the original flags
    sys.setdlopenflags( flags0 | ctypes.RTLD_GLOBAL ) # set our workaround flags

# add to PATH to make dlopen find the libraries
if "PATH" in os.environ:
  os.environ["PATH"] = os.path.dirname(__file__) + os.pathsep +  os.environ["PATH"]
else:
  os.environ["PATH"] = os.path.dirname(__file__)

if "LD_LIBRARY_PATH" in os.environ:
  os.environ["LD_LIBRARY_PATH"] = os.path.dirname(__file__) + os.pathsep + os.environ["LD_LIBRARY_PATH"]
else:
  os.environ["LD_LIBRARY_PATH"] = os.path.dirname(__file__)

from casadi import *    # import everything
import casadi as casadi # import everything
  
if 'casadi_core' in failed_modules:
    raise Exception("Error while loading casadi: %s" % str(failed_modules["casadi_core"]))

if hasattr(sys,"getdlopenflags"):
    sys.setdlopenflags( flags0 ) # set the old flags back
  
import os
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

import warnings
warnings.filterwarnings("default",".*This CasADi function.*",DeprecationWarning)

import contextlib

@contextlib.contextmanager
def internalAPI():
    backup = CasadiOptions.getAllowedInternalAPI()
    CasadiOptions.setAllowedInternalAPI(True)
    yield
    CasadiOptions.setAllowedInternalAPI(backup)

__version__ = CasadiMeta.getVersion()
if '+' in __version__ and CasadiMeta.getGitDescribe()!='':
    __version__  = CasadiMeta.getGitDescribe()
  
