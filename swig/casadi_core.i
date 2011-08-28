%module casadi_core

//  The things needed to make each casadi_*.i  compilable by itself: typemaps
%include "common.i"

//  init hooks
%include "casadi_runtime.i"

// Auxilliary casadi functions:  printing for std::vector, printable_object, shared_object, casadi_types, generic_type, options_functionality
%include "casadi_aux.i"

// SX, Matrix, MX
%include "casadi_primitive.i"

// tools for SX, matrix, MX
%include "casadi_primitive_tools.i"

