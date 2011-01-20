%{
#include "casadi/matrix/crs_sparsity.hpp"
#include "casadi/matrix/matrix.hpp"
#include "casadi/sx/sx.hpp"
#include "casadi/sx/sx_tools.hpp"
%}

%include "typemaps.i"
%include "casadi/matrix/crs_sparsity.hpp"
%include "casadi/matrix/matrix.hpp"
%include "casadi/sx/sx.hpp"

%extend std::vector<CasADi::SX>{
  std::string __repr__(){ return CasADi::getRepresentation(*$self); }
  std::string __str__(){ return CasADi::getDescription(*$self); }
};

%extend std::vector<CasADi::Matrix< CasADi::SX> >{
  std::string __repr__(){ return CasADi::getRepresentation(*$self); }
  std::string __str__(){ return CasADi::getDescription(*$self); }
};

%extend std::vector<std::vector< CasADi::SX> >{
  std::string __repr__(){ return CasADi::getRepresentation(*$self); }
  std::string __str__(){ return CasADi::getDescription(*$self); }
};

namespace CasADi {
  %extend Matrix<double> {
    %pythoncode %{
    def toArray(self):
      import numpy as n
      r = n.array((),dtype=float)
      r.resize(self.size1(),self.size2())
      for i in range(self.size1()):  # loop over rows
        for el in range(self.rowind(i),self.rowind(i+1)): # loop over the non-zero elements
          j=self.col(el)  # column
          r[i,j] = self[el] # add the non-zero element

      return r
    %}
        
    %pythoncode %{
    def toMatrix(self):
      import numpy as n
      return n.matrix(self.toArray())
    %}
    
    %pythoncode %{
    def __getitem__(self,s):
      if isinstance(s,tuple) and len(s)==2:
        s = list(s)
        for k in range(2):
          if isinstance(s[k],slice):
            J = s[k].indices(self.size1())
            s[k] = range(J[0],J[1],J[2])
          elif isinstance(s[k],int):
            s[k] = [s[k]]
      return self.getitem(s)
    %}
    
    %pythoncode %{
    def __setitem__(self,s,val):
      if isinstance(s,tuple) and len(s)==2:
        s = list(s)
        for k in range(2):
          if isinstance(s[k],slice):
            J = s[k].indices(self.size1())
            s[k] = range(J[0],J[1],J[2])
          elif isinstance(s[k],int):
            s[k] = [s[k]]
      self.setitem(s,val)
    %}

  };

%extend Matrix<SX>{
    // The constructor has to be added since SX::operator Matrix<SX does not work
    // Matrix<SX>(const SX&){ *$self
    
    // These methods must be added since the implicit type cast does not work
    Matrix<SX> __pow__ (double b) const{ return $self->__pow__(CasADi::Matrix<CasADi::SX>(b));}
    Matrix<SX> __rpow__(double b) const{ return CasADi::Matrix<CasADi::SX>(b).__pow__(*$self);}
    Matrix<SX> __add__ (double b) const{ return *$self + CasADi::Matrix<CasADi::SX>(b);}
    Matrix<SX> __radd__(double b) const{ return CasADi::Matrix<CasADi::SX>(b) + *$self;}
    Matrix<SX> __sub__ (double b) const{ return *$self - CasADi::Matrix<CasADi::SX>(b);}
    Matrix<SX> __rsub__(double b) const{ return CasADi::Matrix<CasADi::SX>(b) - *$self;}
    Matrix<SX> __mul__ (double b) const{ return *$self * CasADi::Matrix<CasADi::SX>(b);}
    Matrix<SX> __rmul__(double b) const{ return CasADi::Matrix<CasADi::SX>(b) * *$self;}
    Matrix<SX> __div__ (double b) const{ return *$self / CasADi::Matrix<CasADi::SX>(b);}
    Matrix<SX> __rdiv__(double b) const{ return CasADi::Matrix<CasADi::SX>(b) / *$self;}

    %pythoncode %{
    def toArray(self):
      import numpy as n
      r = n.array((),dtype=object)
      r.resize(self.size1(),self.size2())
      for i in range(self.size1()):  # loop over rows
        for el in range(self.rowind(i),self.rowind(i+1)): # loop over the non-zero elements
          j=self.col(el)  # column
          r[i,j] = self[el] # add the non-zero element

      return r
    %}
        
    %pythoncode %{
    def toMatrix(self):
      import numpy as n
      return n.matrix(self.toArray())
    %}
    
    %pythoncode %{
    def __getitem__(self,s):
      if isinstance(s,tuple) and len(s)==2:
        s = list(s)
        for k in range(2):
          if isinstance(s[k],slice):
            J = s[k].indices(self.size1())
            s[k] = range(J[0],J[1],J[2])
          elif isinstance(s[k],int):
            s[k] = [s[k]]
      return self.getitem(s)
    %}

    %pythoncode %{
    def __setitem__(self,s,val):
      if isinstance(s,tuple) and len(s)==2:
        s = list(s)
        for k in range(2):
          if isinstance(s[k],slice):
            J = s[k].indices(self.size1())
            s[k] = range(J[0],J[1],J[2])
          elif isinstance(s[k],int):
            s[k] = [s[k]]
      self.setitem(s,val)
    %}
};
} // namespace CasADi


#ifdef WITH_NUMPY
#include <numpy/arrayobject.h>
#endif // WITH_NUMPY

// Template instantiations
%template(vector_PyObject)    std::vector<PyObject*>;
