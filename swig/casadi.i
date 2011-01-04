%module casadi

#ifdef WITH_DOXDOC
%include "../doc/doc.i"
#endif 

#ifndef WITH_DOXDOC
#warning "Not using doxygen. Run /doc/doc2swig.py to include it. A make rebuild-cache may be necessary."
#endif

%feature("autodoc", "1");

// STL
%include "std_string.i"
%include "std_vector.i"

// Template instantiations
namespace std {
%template(vector_int) vector<int>;
%template(vector_double) vector<double>;
} // namespace std;

#ifndef WITH_NUMPY
#warning "Not using numpy. option(WITH_NUMPY = OFF)"
#endif

#ifdef WITH_NUMPY
#warning "Using numpy. option(WITH_NUMPY = ON)"
%include "numpy.i"
#endif

// Auxilliary casadi functions
%include "casadi_aux.i"

// SX class
%include "sx.i"

// MX class
%include "mx.i"

// FX
%include "fx.i"

// Expression tools
%include "casadi/expression_tools.hpp"

// IPOPT
#ifdef WITH_IPOPT
%include "ipopt_interface.i"
#endif

// Modelica
#ifdef WITH_MODELICA
%include "modelica.i"
#endif

// ACADO
#ifdef WITH_ACADO
%include "acado_interface.i"
#endif

// Sundials
#ifdef WITH_SUNDIALS
%include "sundials_interface.i"
#endif

// SuperLU
#ifdef WITH_SUPERLU
%include "superlu.i"
#endif

// LAPACK
#ifdef WITH_LAPACK
%include "lapack_interface.i"
#endif


