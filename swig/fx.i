%{
#include "casadi/fx/function_io.hpp"
#include "casadi/fx/fx.hpp"
#include "casadi/fx/jacobian.hpp"
#include "casadi/fx/mx_function.hpp"
#include "casadi/fx/sx_function.hpp"
#include "casadi/fx/integrator.hpp"
#include "casadi/fx/simulator.hpp"
#include "casadi/fx/nlp_solver.hpp"
#include "casadi/fx/external_function.hpp"
#include "casadi/fx/fx_tools.hpp"
#include "casadi/fx/parallelizer.hpp"
%}

%include "casadi/fx/function_io.hpp"
%include "casadi/fx/fx.hpp"
%include "casadi/fx/jacobian.hpp"
%include "casadi/fx/sx_function.hpp"
%include "casadi/fx/mx_function.hpp"
%include "casadi/fx/linear_solver.hpp"
%include "casadi/fx/integrator_jacobian.hpp"
%include "casadi/fx/integrator.hpp"
%include "casadi/fx/simulator.hpp"
%include "casadi/fx/nlp_solver.hpp"
%include "casadi/fx/external_function.hpp"
%include "casadi/fx/fx_tools.hpp"
%include "casadi/fx/parallelizer.hpp"

%template(IntegratorVector) std::vector<CasADi::Integrator>;

