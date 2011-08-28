#ifdef WITH_SWIG_SPLIT
%module casadi_primitive

%include "common.i"

%import "casadi_main.i"
#endif //WITH_SWIG_SPLIT


#ifdef SWIGPYTHON
#ifdef WITH_SWIG_SPLIT
%pythoncode %{
import _casadi_primitive_tools as _casadi_global
%}
#endif // WITH_SWIG_SPLIT
#ifndef WITH_SWIG_SPLIT
%pythoncode %{
_casadi_global = _casadi
%}
#endif // WITH_SWIG_SPLIT
#endif // SWIGPYTHON

// Matrix typemaps class
%include "matrix.i"

// SX class
%include "sx.i"

// MX class
%include "mx.i"
