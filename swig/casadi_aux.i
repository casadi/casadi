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
%include "casadi/generic_type.hpp"

namespace CasADi {
%typemap(out) GenericType {
if ($1.isBool()) {
  $result=PyBool_FromLong($1.toBool());
} else if ($1.isInt()) {
  $result=PyInt_FromLong($1.toInt());
} else if ($1.isDouble()) {
  $result=PyFloat_FromDouble($1.toDouble());
} else if ($1.isString()) {
  $result=PyString_FromString($1.toString().c_str());
} else {
  SWIG_exception_fail(SWIG_TypeError,"GenericType not yet implemented");
}
}
}

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

