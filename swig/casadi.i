%module casadi

// STL
%include "std_string.i"
%include "std_vector.i"

// Template instantiations
namespace std {
%template(vector_int) vector<int>;
%template(vector_double) vector<double>;
} // namespace std;

// Auxilliary casadi functions
%include "casadi_aux.i"

// SX class
%include "sx.i"

// MX class
%include "mx.i"

// FX
%include "fx.i"

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


