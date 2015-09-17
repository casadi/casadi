%module(package="extending_casadi") extending_casadi

// WORKAROUNDS BEGINS: Due to Python-related issues in casadi.i
#define CASADI_NOT_IN_DERIVED
#ifdef SWIGPYTHON
%{
#define SWIG_FILE_WITH_INIT
#include "numpy.hpp"
#define SWIG_PYTHON_CAST_MODE 1
%}

%init %{
import_array();
%}
#endif // SWIGPYTHON
// WORKAROUNDS END

%import "../casadi.i"

%{
#include <swig/extending_casadi/extending_casadi.hpp>
%}

%include <swig/extending_casadi/extending_casadi.hpp>

