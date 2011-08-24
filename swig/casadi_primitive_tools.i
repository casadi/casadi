// Matrix tools
%include "casadi/matrix/matrix_tools.hpp"

// Instsantiate the functions
MATRIX_TOOLS_TEMPLATES(double)
MATRIX_TOOLS_TEMPLATES(CasADi::SX)

// Sparsity tools
%{
#include "casadi/matrix/sparsity_tools.hpp"
%}
%include "casadi/matrix/sparsity_tools.hpp"

// SX tools
%include "casadi/sx/sx_tools.hpp"
