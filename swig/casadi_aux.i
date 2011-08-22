%{
#include <sstream>
#include "casadi/stl_vector_tools.hpp"
#include "casadi/printable_object.hpp"
#include "casadi/shared_object.hpp"
#include "casadi/generic_type.hpp"
#include "casadi/options_functionality.hpp"
%}

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

%include "casadi/printable_object.hpp"
%include "casadi/shared_object.hpp"


#ifdef SWIGPYTHON

namespace CasADi {
%typemap(out) GenericType {
bool ret=meta<  CasADi::GenericType >::toPython($1,$result);
if (!ret) {
  SWIG_exception_fail(SWIG_TypeError,"GenericType not yet implemented");
}
}

%typemap(out) const GenericType::Dictionary&  {
bool ret=meta<  CasADi::GenericType::Dictionary >::toPython(*$1,$result);
if (!ret) {
  SWIG_exception_fail(SWIG_TypeError,"GenericType not yet implemented");
}
}

} // namespace CasADi
#endif // SWIGPYTHON

%my_generic_const_typemap(PRECEDENCE_GENERICTYPE,CasADi::GenericType)
#ifdef SWIGPYTHON
%my_generic_const_typemap(PRECEDENCE_DICTIONARY ,CasADi::GenericType::Dictionary)
#endif

#ifdef SWIGOCTAVE
%my_generic_const_typemap(PRECEDENCE_DVector,std::vector<double>);
%my_generic_const_typemap(PRECEDENCE_IVector,std::vector<int>);
#endif //SWIGOCTAVE

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

