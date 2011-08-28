#ifdef WITH_SWIG_SPLIT
%module casadi_noncore

%include "common.i"

%import "casadi_core.i"
#endif //WITH_SWIG_SPLIT

// FX
%include "casadi_fx.i"

// optimal_control
%include "casadi_optimal_control.i"

// interfaces
%include "casadi_interfaces.i"


