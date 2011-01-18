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
  std::string __repr__(){ return repr(*$self); }
  std::string __str__(){ return print(*$self); }
};

%extend std::vector<CasADi::Matrix< CasADi::SX> >{
  std::string __repr__(){ return repr(*$self); }
  std::string __str__(){ return print(*$self); }
};

%extend std::vector<std::vector< CasADi::SX> >{
  std::string __repr__(){ return repr(*$self); }
  std::string __str__(){ return print(*$self); }
};

#ifdef WITH_NUMPY
#include <numpy/arrayobject.h>
#endif // WITH_NUMPY

// Template instantiations
%template(vector_PyObject)    std::vector<PyObject*>;
