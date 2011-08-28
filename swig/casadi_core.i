#ifdef WITH_SWIG_SPLIT
%module casadi_core

%include "common.i"

%import "casadi_main.i"
#endif //WITH_SWIG_SPLIT

// SX, Matrix, MX
%include "casadi_primitive.i"

// tools for SX, matrix, MX
%include "casadi_primitive_tools.i"

