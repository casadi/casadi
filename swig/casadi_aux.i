%{
#include <sstream>
#include "casadi/stl_vector_tools.hpp"
#include "casadi/printable_object.hpp"
#include "casadi/shared_object.hpp"
#include "casadi/generic_type.hpp"
#include "casadi/options_functionality.hpp"
%}

%include "casadi/printable_object.hpp"
%include "casadi/shared_object.hpp"
%include "casadi/casadi_types.hpp"
%include "casadi/generic_type.hpp"
%include "casadi/options_functionality.hpp"

namespace CasADi {
  %extend OptionsFunctionality {
    void setOption(const std::string &name, const std::string& val){$self->setOption(name,val);} 
    void setOption(const std::string &name, const std::vector<int>& val){$self->setOption(name,val);} 
    void setOption(const std::string &name, const std::vector<double>& val){$self->setOption(name,val);} 
    void setOption(const std::string &name, double val){$self->setOption(name,val);}
    void setOption(const std::string &name, int val){$self->setOption(name,val);} 
    void setOption(const std::string &name, bool val){$self->setOption(name,val);}  
  }
} // namespace CasADi

