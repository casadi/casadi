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

import ctypes

# add to PATH to make dlopen find the libraries
if "PATH" in os.environ:
  os.environ["PATH"] = os.path.dirname(__file__) + os.pathsep +  os.environ["PATH"]
else:
  os.environ["PATH"] = os.path.dirname(__file__)

if "LD_LIBRARY_PATH" in os.environ:
  os.environ["LD_LIBRARY_PATH"] = os.path.dirname(__file__) + os.pathsep + os.environ["LD_LIBRARY_PATH"]
else:
  os.environ["LD_LIBRARY_PATH"] = os.path.dirname(__file__)

from casadi_loader import *    # import everything
import casadi_loader as casadi # import everything
  
if 'casadi_core' in failed_modules:
    raise Exception("Error while loading casadi: %s" % str(failed_modules["casadi_core"]))

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


import inspect
import re


def pythonify(s):

  s = s.replace("C/C++ prototypes","Python usages")
  s = s.replace("casadi::","")
  s = s.replace("std::string","str")
  s = s.replace(" const &","")
  s = re.sub("(const )?Matrix< SXElement >( &)?",r"SX",s)
  s = re.sub("(const )?GenericMatrix< ?(\w+) *>( ?&)?",r"\2 ",s)
  s = re.sub("(const )?Matrix< ?(\w+) *>( ?&)?",r"array(\2) ",s)
  s = re.sub("(const )?GenericMatrix< ?([\w\(\)]+) *>( ?&)?",r"\2 ",s)
  s = re.sub(r"const (\w+) &",r"\1 ",s)
  s = re.sub(r"< [\w\(\)]+ +>\(",r"(",s)
  s = re.sub(r"\b(\w+)(< \w+ >)?::\1",r"\1",s)
  s = re.sub(r"(const )? ?std::pair< ?([\w\(\) ]+?) ?, ?([\w\(\) ]+?) ?> ?&?",r"(\2,\3) ",s)
  for i in range(5):
    s = re.sub(r"(const )? ?std::vector< ?([\w\(\)\[\] ]+) ?(, ?std::allocator< ?\2 ?>)? ?> ?&?",r"[\2] ",s)
  s = re.sub(r"StructIOSchemeVector(< ?([\w\[\] ]+) ?>)?",r"Structure",s)
  s = re.sub(r"IOSchemeVector< ?(\w+) ?>",r"scheme(\1)",s)
  s = re.sub(r"\b(\w+)(< \w+ >)?::\1",r"\1",s)
  s = s.replace("casadi::","")
  s = s.replace("IOInterface< Function >","Function")
  s = s.replace("::",".")
  s = s.replace(".operator ()","")
  s = re.sub(r"([A-Z]\w+)Vector",r"[\1]",s)
  return s
  
def type_descr(a):
  if isinstance(a,list):
    if len(a)>0:
      return "[%s]" % ",".join(set([type_descr(i) for i in a]))
    else:
      return "[]"
  elif isinstance(a,tuple):
    return "(%s)" % ",".join([type_descr(i) for i in a])
  elif isinstance(a,np.ndarray):
    return "np.array(%s)" % ",".join(set([type_descr(i) for i in np.array(a).flatten().tolist()])) 
  if type(a).__name__.startswith("IOSchemeVector"):
    return "scheme(%s)" % type(a).__name__[len("IOSchemeVector"):]
  else:
    return type(a).__name__

def monkeypatch(v,cl=True):
  if hasattr(v,"__monkeypatched__"):
    return v
  def foo(*args,**kwargs):
    try:
      return v(*args,**kwargs)
    except NotImplementedError as e:
      import sys
      exc_info = sys.exc_info()
      if e.message.startswith("Wrong number or type of arguments for overloaded function"):

        s = e.args[0]
        s = s.replace("'new_","'")
        s = re.sub(r"overloaded function '(\w+?)_(\w+)'",r"overloaded function '\1.\2'",s)
        m = re.search("overloaded function '([\w\.]+)'",s)
        if m:
          name = m.group(1)
          name = name.replace(".__call__","")
        else:
          name = "method"
        ne = NotImplementedError(pythonify(s)+"You have: %s(%s)\n" % (name,", ".join(map(type_descr,args[1:] if cl else args)+ ["%s=%s" % (k,type_descr(vv)) for k,vv in kwargs.items()])))
        raise ne.__class__, ne, exc_info[2].tb_next
      else:
        raise exc_info[1], None, exc_info[2].tb_next
    except TypeError as e:
      import sys
      exc_info = sys.exc_info()
      
      methodname = "method"
      try:
        methodname = exc_info[2].tb_next.tb_frame.f_code.co_name
      except:
        pass

      if e.message.startswith("in method '"):
        s = e.args[0]
        s = re.sub(r"method '(\w+?)_(\w+)'",r"method '\1.\2'",s)
        m = re.search("method '([\w\.]+)'",s)
        if m:
          name = m.group(1)
          name = name.replace(".__call__","")
        else:
          name = "method"
        ne = TypeError(pythonify(s)+" expected.\nYou have: %s(%s)\n" % (name,", ".join(map(type_descr,args[1:] if cl else args))))
        raise ne.__class__, ne, exc_info[2].tb_next
      elif e.message.startswith("Expecting one of"):
        s = e.args[0]
        conversion = {"mul": "*", "div": "/", "add": "+", "sub": "-","le":"<=","ge":">=","lt":"<","gt":">","eq":"==","pow":"**"}
        if methodname.startswith("__") and methodname[2:-2] in conversion:
          ne = TypeError(pythonify(s)+"\nYou try to do: %s %s %s.\n" % (  type_descr(args[0]),conversion[methodname[2:-2]] ,type_descr(args[1]) ))
        elif methodname.startswith("__r") and methodname[3:-2] in conversion:
          ne = TypeError(pythonify(s)+"\nYou try to do: %s %s %s.\n" % ( type_descr(args[1]),  conversion[methodname[3:-2]], type_descr(args[0]) ))
        else:
          ne = TypeError(pythonify(s)+"\nYou have: (%s)\n" % (", ".join(map(type_descr,args[1:] if cl else args))))
        raise ne.__class__, ne, exc_info[2].tb_next
      else:
        s = e.args[0]
        ne = TypeError(s+"\nYou have: (%s)\n" % (", ".join(map(type_descr,args[1:] if cl else args) + ["%s=%s" % (k,type_descr(vv)) for k,vv in kwargs.items()]  )))
        raise ne.__class__, ne, exc_info[2].tb_next
    except Exception as e:
      import sys
      exc_info = sys.exc_info()
      raise exc_info[1], None, exc_info[2].tb_next
      
  if v.__doc__ is not None:
    foo.__doc__ = pythonify(v.__doc__)
  foo.__name__ = v.__name__
  foo.__monkeypatched__ = True
  return foo

def improvedcall(v):
  def newcall(self,*args,**kwargs):
    if len(args)>0 and len(kwargs)>0:
      raise Exception("You cannot mix positional and keyword arguments in __call__")
    if len(kwargs)>0:
      scheme = self.getInputScheme()
      if scheme.known():
        return v(self,scheme(**kwargs))
      else:
        raise Exception("The CasADi Function that you call has no known input scheme. You cannot use keyword arguments, just possitional arguments.")
    else:
      return v(self,*args)

  newcall.__name__ = v.__name__
  newcall.__doc__ = v.__doc__ + "\nYou can also call with keyword arguments if the Function has a known scheme\nExample: nlp(x=x)\n"
  return newcall
  
for name,cl in inspect.getmembers(casadi, inspect.isclass):
  for k,v in inspect.getmembers(cl, inspect.ismethod):
    if k=="__init__" and  "SchemeVectorSX" in name: continue
    if k == "__del__" or v.__name__ == "<lambda>": continue
    vv = v
    if k=="__call__" and issubclass(cl,Function):
      vv = improvedcall(v)
    setattr(cl,k,monkeypatch(vv))
  for k,v in inspect.getmembers(cl, inspect.isfunction):
    setattr(cl,k,staticmethod(monkeypatch(v,cl=False)))
  
for name,v in inspect.getmembers(casadi, inspect.isfunction):
  p = monkeypatch(v,cl=False)
  setattr(casadi,name,p)
  import sys
  setattr(sys.modules[__name__], name, p)
  
  
class IOSchemeVectorExtractor(object):
  def __init__(self,schemevector):
    self._list = []
    self._schemevector = schemevector
    
  def __getattr__(self, name):
    self._list.append(name)
    return self
  
  def __iter__(self):
    for k in self._list:
      try:
        yield self._schemevector[k]
      except:
        import sys
        exc_info = sys.exc_info()
        raise exc_info[1], None, exc_info[2].tb_next
    
def extract(self):
  """
    This is a convenience function to extract multiple outputs from a scheme.
    
    Use::
      
      [xf,qf] = integrator(x0=x).get.xf.qf
      
    If you need just a single output, you could also use the slicing approach::
    
      xf = integrator(x0=x)["xf"]
    
  """
  return IOSchemeVectorExtractor(self)

for name,cl in inspect.getmembers(casadi, inspect.isclass):
  if "IOSchemeVector" in name:
    setattr(cl,"get",property(extract))
    getattr(cl,"__swig_getmethods__")["get"] = extract
    
@contextlib.contextmanager
def internalAPI():
    backup = CasadiOptions.getAllowedInternalAPI()
    CasadiOptions.setAllowedInternalAPI(True)
    try:
      yield
    finally:
      CasadiOptions.setAllowedInternalAPI(backup)

__version__ = CasadiMeta.getVersion()
if '+' in __version__ and CasadiMeta.getGitDescribe()!='':
    __version__  = CasadiMeta.getGitDescribe()
  
