%module(package="extending_casadi") extending_casadi

// WORKAROUNDS BEGINS: Due to Python-related issues in casadi.i
#ifdef SWIGPYTHON
%{
#define SWIG_FILE_WITH_INIT
#include "casadi_numpy.hpp"
%}

%init %{
import_array();
%}
#endif // SWIGPYTHON
// WORKAROUNDS END

%include "exception.i"
%import "../casadi.i"
%exception {
  try {
    $action
   } catch(const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

%{
#include <swig/extending_casadi/extending_casadi.hpp>
%}

%include <swig/extending_casadi/extending_casadi.hpp>

