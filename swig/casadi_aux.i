%{
#include <sstream>
#include "casadi/stl_vector_tools.hpp"
#include "casadi/printable_object.hpp"
#include "casadi/shared_object.hpp"
#include "casadi/options_functionality.hpp"
%}

%include "casadi/printable_object.hpp"
%include "casadi/shared_object.hpp"
%include "casadi/options_functionality.hpp"

// Exceptions handling
%include "exception.i"
%exception {
try {
  $action
  } catch (const std::exception& e) {
  SWIG_exception(SWIG_RuntimeError, e.what());
  } catch (const char* e) { // depreciated!!
    SWIG_exception(SWIG_RuntimeError, e);
  }
}

namespace CasADi {
  %extend OptionsFunctionality {
    void setOption(const std::string &name, const std::string& val){$self->setOption(name,val);} 
    void setOption(const std::string &name, const std::vector<double>& val){$self->setOption(name,val);} 
    void setOption(const std::string &name, double val){$self->setOption(name,val);} 
  }
} // namespace CasADi

