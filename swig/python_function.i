/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


namespace CasADi{
%extend CFunction {
  void __setUserData__(PyObject * obj) {
    $self->setOption("user_data",static_cast<void*>(obj));
  }
}
}

%pythoncode %{


class PyFunction(CFunction):
   def __init__(self,fun, inputscheme=None, outputscheme=None):
     if inputscheme == None:
        inputscheme = []
     if outputscheme == None:
        outputscheme = []
     CFunction.__init__(self,py_c_wrapper,inputscheme,outputscheme)
     self.__setUserData__(self)
     self.__fun__ = fun
     self.userdata = {}
   def __call__(self,a,b):
     return self.__fun__(self,a,b,self.userdata)
     
   def setUserData(self,data):
     self.userdata = data
%}

%{
#include "swig/python_function.h"
%}
%include "swig/python_function.h"
