
namespace CasADi{
%extend CFunction {
  void __setUserData__(PyObject * obj) {
    $self->setUserData(obj);
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
