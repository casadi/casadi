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

import warnings
import os
import numpy as np

import sys
if sys.version_info >= (3, 0):
  from casadi.casadi import *
else:
  from casadi import *
  import casadi

# For plugin loading
GlobalOptions.setCasadiPath(os.path.dirname(__file__))

import types

def wrapper(f, warning,error=False):
    def new(*args, **kwargs):
        print(("*" * 40))
        print("Deprecation Warning")
        print(("-" * 40))
        print(warning)
        print(("*" * 40))
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


import inspect
import re


__version__ = CasadiMeta.version()
if '+' in __version__ and CasadiMeta.git_describe()!='':
    __version__  = CasadiMeta.git_describe()
