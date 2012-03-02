%{
#include "optimal_control/variable.hpp"
#include "optimal_control/flat_ocp.hpp"
#include "optimal_control/variable_tools.hpp"
#include "optimal_control/ocp_tools.hpp"
#include "optimal_control/multiple_shooting.hpp"
#include "optimal_control/collocation.hpp"
%}

%include "optimal_control/variable.hpp"
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

/*%extend CasADi::OptimalControl::OCP{
    %pythoncode %{
      def getSubTree(self,vars):
        if vars.var_.isNull():
          varstruct = VariableStruct()
          names = vars.getNames()
          for i in range(len(names)):
            varstruct.__setattr__(names[i],self.getSubTree(vars.subByIndex(i+1)))
          return varstruct
        else:
          return vars.var_
      
      def getVariables(self):
        return self.getSubTree(self.variables_)

    %}
};*/
#endif // SWIGPYTHON


