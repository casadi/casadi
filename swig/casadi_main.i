#ifdef WITH_SWIG_SPLIT
%module casadi_main

%include "common.i"

#endif //WITH_SWIG_SPLIT

//  init hooks
%include "casadi_runtime.i"

// Auxilliary casadi functions:  printing for std::vector, printable_object, shared_object, casadi_types, generic_type, options_functionality
%include "casadi_aux.i"
