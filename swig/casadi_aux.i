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

#define memberbinopsr(Type,uname) \
Type __r##uname##__(const Type& b) const{ return b.__##uname##__(*$self);}

#define memberbinopsr_un(Type,uname) \
Type __r##uname##__(const Type& b) const{ return b.##uname##(*$self);}

#define memberbinopsr_nn(Type,uname) \
Type r##uname##(const Type& b) const{ return b.##uname##(*$self);}

#define binopsrFull(Type) \
memberbinopsr(Type,pow) \
memberbinopsr(Type,add) \
memberbinopsr(Type,sub) \
memberbinopsr(Type,mul) \
memberbinopsr(Type,div) \
memberbinopsr(Type,mldivide) \
memberbinopsr(Type,mrdivide) \
memberbinopsr(Type,mpower) \
memberbinopsr(Type,constpow) \
memberbinopsr_un(Type,fmin) \
memberbinopsr_un(Type,fmax) \
memberbinopsr_nn(Type,prod)

#define memberbinops(uname,argtype,argCast,selfCast,returntype) \
returntype __##uname##__ (argtype) const{ return selfCast(*$self).__##uname##__(argCast(b));} \
returntype __r##uname##__(argtype) const{ return argCast(b).__##uname##__(selfCast(*$self));} \

// These methods must be added since the implicit type cast does not work.
// Consider a+b  with a DMatrix and b SXMatrix
// In C++, operator+(SXMatrix,SXMatrix) will be called (implicit cast)
// In octave, a.__add__(b) will be called   (no implicit cast)
// In python, __array_priority__ will be checked and b.__radd__(a) will be called (effectively implicit casting)

// This is a list of all operators:
#define binopsFull(argtype,argCast,selfCast,returntype) \
memberbinops(pow,argtype,argCast,selfCast,returntype) \
memberbinops(add,argtype,argCast,selfCast,returntype) \
memberbinops(sub,argtype,argCast,selfCast,returntype) \
memberbinops(mul,argtype,argCast,selfCast,returntype) \
memberbinops(div,argtype,argCast,selfCast,returntype) \
returntype prod (argtype) const{ return prod(selfCast(*$self) , argCast(b));} \
returntype rprod (argtype) const{ return prod(argCast(b) , selfCast(*$self));} \
memberbinops(mldivide,argtype,argCast,selfCast,returntype) \
memberbinops(mrdivide,argtype,argCast,selfCast,returntype) \
memberbinops(mpower,argtype,argCast,selfCast,returntype) 

// This is a list of operators that do not check __array_priority__ in python
#define binopsNoPriority(argtype,argCast,selfCast,returntype) \
memberbinops(pow,argtype,argCast,selfCast,returntype) \
