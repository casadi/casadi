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
      if isinstance(s,slice):
        s = (s,[0])
      if isinstance(s,tuple) and len(s)==2 and (isinstance(s[1],slice) or isinstance(s[0],slice)):
        s = list(s)
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
      if isinstance(s,slice):
        s = (s,[0])
      if isinstance(s,tuple) and len(s)==2 and (isinstance(s[1],slice) or isinstance(s[0],slice)):
        s = list(s)
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

  };
} // namespace CasADi
