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


%extend CasADi::Matrix<CasADi::SX>{
  // The constructor has to be added since SX::operator Matrix<SX does not work
  // CasADi::Matrix<CasADi::SX>(const CasADi::SX&){ *$self
  
  // These methods must be added since the implicit type cast does not work
  CasADi::Matrix<CasADi::SX> __pow__ (double b) const{ return $self->__pow__(CasADi::Matrix<CasADi::SX>(b));}
  CasADi::Matrix<CasADi::SX> __rpow__(double b) const{ return CasADi::Matrix<CasADi::SX>(b).__pow__(*$self);}
  CasADi::Matrix<CasADi::SX> __add__ (double b) const{ return *$self + CasADi::Matrix<CasADi::SX>(b);}
  CasADi::Matrix<CasADi::SX> __radd__(double b) const{ return CasADi::Matrix<CasADi::SX>(b) + *$self;}
  CasADi::Matrix<CasADi::SX> __sub__ (double b) const{ return *$self - CasADi::Matrix<CasADi::SX>(b);}
  CasADi::Matrix<CasADi::SX> __rsub__(double b) const{ return CasADi::Matrix<CasADi::SX>(b) - *$self;}
  CasADi::Matrix<CasADi::SX> __mul__ (double b) const{ return *$self * CasADi::Matrix<CasADi::SX>(b);}
  CasADi::Matrix<CasADi::SX> __rmul__(double b) const{ return CasADi::Matrix<CasADi::SX>(b) * *$self;}
  CasADi::Matrix<CasADi::SX> __div__ (double b) const{ return *$self / CasADi::Matrix<CasADi::SX>(b);}
  CasADi::Matrix<CasADi::SX> __rdiv__(double b) const{ return CasADi::Matrix<CasADi::SX>(b) / *$self;}
};





#ifdef WITH_NUMPY
#include <numpy/arrayobject.h>
#endif // WITH_NUMPY

// Template instantiations
%template(vector_PyObject)    std::vector<PyObject*>;
