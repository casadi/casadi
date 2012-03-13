%{
#include "optimal_control/variable.hpp"
#include "optimal_control/symbolic_ocp.hpp"
#include "optimal_control/flat_ocp.hpp"
#include "optimal_control/variable_tools.hpp"
#include "optimal_control/ocp_tools.hpp"
#include "optimal_control/multiple_shooting.hpp"
#include "optimal_control/collocation.hpp"
%}

%include "optimal_control/variable.hpp"
%include "optimal_control/symbolic_ocp.hpp"
%include "optimal_control/flat_ocp.hpp"
%include "optimal_control/variable_tools.hpp"
%include "optimal_control/ocp_tools.hpp"
%include "optimal_control/multiple_shooting.hpp"
%include "optimal_control/collocation.hpp"

#ifdef SWIGPYTHON
%pythoncode %{
  class VariableStruct(object):
    """Structure for browsing through a variable tree."""
    def __repr__(self):
      return repr(self.__dict__)
%}

#endif // SWIGPYTHON


