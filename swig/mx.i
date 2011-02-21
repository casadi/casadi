%{
#include "casadi/mx/mx.hpp"
#include "casadi/mx/mx_tools.hpp"
%}

%include "casadi/mx/mx.hpp"
%include "casadi/mx/mx_tools.hpp"

namespace CasADi{
  %extend MX{
    %pythoncode %{
    def __getitem__(self,s):
      if isinstance(s,int):
        if s < 0:
          s = s + self.size()
      elif isinstance(s,slice):
        s = (s,[0])
      if isinstance(s,tuple):
        if len(s)!=2:
          raise Exception("get/setitem can only do 1D or 2D indexing")
        s = list(s)
        if isinstance(s[0],int) and isinstance(s[1],int):
          for k in range(2):
            if s[k]<0:
              s[k]=s[k]+self.shape[k]
        else:
          for k in range(2):
            if isinstance(s[k],slice):
              J = s[k].indices(self.shape[k])
              s[k] = range(J[0],J[1],J[2])
            elif isinstance(s[k],int):
              if s[k]<0:
                s[k]=s[k]+self.shape[k]
              s[k] = [s[k]]
      return self.getitem(s)
    %}

    %pythoncode %{
    def __setitem__(self,s,val):
      if isinstance(s,int):
        if s < 0:
          s = s + self.size()
      elif isinstance(s,slice):
        s = (s,[0])
      if isinstance(s,tuple):
        if len(s)!=2:
          raise Exception("get/setitem can only do 1D or 2D indexing")
        s = list(s)
        if isinstance(s[0],int) and isinstance(s[1],int):
          for k in range(2):
            if s[k]<0:
              s[k]=s[k]+self.shape[k]
        else:
          for k in range(2):
            if isinstance(s[k],slice):
              J = s[k].indices(self.shape[k])
              s[k] = range(J[0],J[1],J[2])
            elif isinstance(s[k],int):
              if s[k]<0:
                s[k]=s[k]+self.shape[k]
              s[k] = [s[k]]
      self.setitem(s,val)
    %}


    %pythoncode %{
        @property
        def shape(self):
            return (self.size1(),self.size2())
    %}
    
    %pythoncode %{
        @property
        def T(self):
            return trans(self)
    %}

  };
} // namespace CasADi
