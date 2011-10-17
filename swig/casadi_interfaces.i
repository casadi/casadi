
#ifdef WITH_OOQP
%{
#include "interfaces/ooqp/ooqp_solver.hpp"
%}
%include "interfaces/ooqp/ooqp_solver.hpp"
#endif

// IPOPT
#ifdef WITH_IPOPT
%include "ipopt_interface.i"
#endif

// ACADO
#ifdef WITH_ACADO
%include "acado_interface.i"
%{
#include "interfaces/qpoases/qpoases_solver.hpp"
%}
%include "interfaces/qpoases/qpoases_solver.hpp"
#endif

// Sundials
#ifdef WITH_SUNDIALS
%include "sundials_interface.i"
#endif

// SuperLU
#ifdef WITH_SUPERLU
%include "superlu.i"
#endif

#ifdef WITH_LAPACK
%include "lapack_interface.i"
#endif

#ifdef WITH_KNITRO
%{ 
  #include "interfaces/knitro/knitro_solver.hpp"
%}
%include "interfaces/knitro/knitro_solver.hpp"
#endif

#ifdef WITH_CPLEX
%{ 
  #include "interfaces/cplex/cplex_solver.hpp"
%}
#include "interfaces/cplex/cplex_solver.hpp"
#endif

#ifdef WITH_LIFTOPT
%{
  #include "interfaces/liftopt/liftopt_solver.hpp"
%}
%include "interfaces/liftopt/liftopt_solver.hpp"
#endif

#ifdef WITH_CSPARSE
%{
#include "interfaces/csparse/csparse.hpp"
%}
%include "interfaces/csparse/csparse.hpp"
#endif

#ifdef WITH_GSL
%{
#include "interfaces/gsl/gsl_integrator.hpp"
%}
%include "interfaces/gsl/gsl_integrator.hpp"
#endif
