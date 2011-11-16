// Turn off the warnings that certain methods are effectively ignored, this seams to be a false warning, 
// for example vertcat(SXMatrixVector), vertcat(DMatrixVector) and vertcat(MXVector) appears to work fine
#pragma SWIG nowarn=509,303

%include "doc.i"

%feature("autodoc", "1");

// Make sure that a copy constructor is created
%copyctor;

// STL
%include "std_string.i"
%include "std_vector.i"
%include "std_pair.i"



#ifdef SWIG_MAIN_MODULE
%template(StringVector) std::vector<std::string>;

%template(BVector)             std::vector<bool> ;
%template(BVectorVector)       std::vector<std::vector<bool> > ;
%template(BVectorVectorVector) std::vector< std::vector<std::vector<bool> > > ;

%template(IVector)             std::vector<int>;
%template(IVectorVector)       std::vector<std::vector<int> > ;
%template(IVectorVectorVector) std::vector< std::vector<std::vector<int> > > ;

%template(DVector)             std::vector<double>;
%template(DVectorVector)       std::vector<std::vector<double> > ;
%template(DVectorVectorVector) std::vector< std::vector<std::vector<double> > > ;

%template(Pair_Int_Int) std::pair<int,int>;
%template(VectorPair_Int_Int) std::vector< std::pair<int,int> >;
#endif //SWIG_MAIN_MODULE
#ifndef SWIG_MAIN_MODULE
%template() std::vector<std::string>;

%template() std::vector<std::vector<bool> > ;
%template() std::vector< std::vector<std::vector<bool> > > ;

%template() std::vector<int>;
%template() std::vector<std::vector<int> > ;
%template() std::vector< std::vector<std::vector<int> > > ;

%template() std::vector<double>;
%template() std::vector<std::vector<double> > ;
%template() std::vector< std::vector<std::vector<double> > > ;

%template() std::pair<int,int>;
%template() std::vector< std::pair<int,int> >;
#endif //SWIG_MAIN_MODULE


// The following is a work-around since it appears not possible to use the standard print functions from stl_vector tools,
// nor the std::stringstream class, since these are included _after_ std::vector in the C++ generated wrapper code
%extend std::vector<double>{  
  std::string __repr__(){
    char buffer[32];
    std::string ret;
    ret += "[";
    for(int i=0; i<$self->size(); ++i){
      if(i!=0) ret += ",";

      // Print to buffer and append
      snprintf(buffer, 32, "%g", $self->at(i));
      ret += buffer;
    }
    ret += "]";
    return ret;
  }
  std::string __str__(){
    char buffer[32];
    std::string ret;
    
    // Add dimension
    snprintf(buffer, 32, "[%lu]", (unsigned long) $self->size());
    ret += buffer; 
    ret += "(";
    for(int i=0; i<$self->size(); ++i){
      if(i!=0) ret += ",";

      // Print to buffer and append
      snprintf(buffer, 32, "%g", $self->at(i));
      ret += buffer;
    }
    ret += ")";
    return ret;
  }
};

// The same workaround for vector<double>
%extend std::vector<int>{  
  std::string __repr__(){
    char buffer[32];
    std::string ret;
    ret += "[";
    for(int i=0; i<$self->size(); ++i){
      if(i!=0) ret += ",";

      // Print to buffer and append
      snprintf(buffer, 32, "%d", $self->at(i));
      ret += buffer;
    }
    ret += "]";
    return ret;
  }
  std::string __str__(){
    char buffer[32];
    std::string ret;
    
    // Add dimension
    snprintf(buffer, 32, "[%lu]", (unsigned long) $self->size());
    ret += buffer; 
    ret += "(";
    for(int i=0; i<$self->size(); ++i){
      if(i!=0) ret += ",";

      // Print to buffer and append
      snprintf(buffer, 32, "%d", $self->at(i));
      ret += buffer;
    }
    ret += ")";
    return ret;
  }
};

#ifdef SWIGOCTAVE
%inline %{
#include <octave/dim-vector.h>
%}
#endif

#ifdef SWIGPYTHON
#ifndef WITH_NUMPY
//#warning "Not using numpy. option(WITH_NUMPY = OFF)"
#endif


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

#ifdef WITH_NUMPY
//#warning "Using numpy. option(WITH_NUMPY = ON)"
%{
#define SWIG_FILE_WITH_INIT
%}
// Get the NumPy typemaps
%{
//#define NO_IMPORT_ARRAY
#include "numpy.hpp"
%}
%{
#define SWIG_PYTHON_CAST_MODE 1
%}

%init %{
import_array();
%}
#endif // WITH_NUMPY
#endif // SWIGPYTHON

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
memberbinops(mldivide,argtype,argCast,selfCast,returntype) \
memberbinops(mrdivide,argtype,argCast,selfCast,returntype) \
memberbinops(mpower,argtype,argCast,selfCast,returntype) 

// This is a list of operators that do not check __array_priority__ in python
#define binopsNoPriority(argtype,argCast,selfCast,returntype) \
memberbinops(pow,argtype,argCast,selfCast,returntype) \

// typemaphelpers
%include "typemaphelpers.i"

// typemap meta implementations
%include "meta.i"

// common typemaps
%include "commontypemaps.i"

%{
#include "casadi/matrix/matrix.hpp" 
#include "casadi/matrix/matrix_tools.hpp" 
#include "casadi/matrix/sparsity_tools.hpp" 
	 
// Scalar expressions 
#include "casadi/sx/sx.hpp" 
#include "casadi/sx/sx_tools.hpp" 
#include "casadi/fx/sx_function.hpp" 
	 
// Matrix expressions 
#include "casadi/mx/mx.hpp" 
#include "casadi/mx/mx_tools.hpp" 

#include "casadi/fx/mx_function.hpp" 
 	
#include "casadi/fx/mx_function.hpp"
#include "casadi/fx/c_function.hpp"
#include "casadi/fx/jacobian.hpp"
#include "casadi/fx/ocp_solver.hpp"
#include "casadi/fx/simulator.hpp"
#include "casadi/fx/parallelizer.hpp"
#include "casadi/fx/external_function.hpp"


#include "optimal_control/ocp_tools.hpp"
#include "optimal_control/multiple_shooting.hpp"
#include "optimal_control/flat_ocp.hpp"

%}

%{
namespace std {
void dummy(CasADi::SX foo,
	std::vector< std::vector<double> > foo1,
	std::vector<double> &foo2,
	std::vector<CasADi::MX> &foo3,
	CasADi::MX foo4,
	CasADi::Matrix<double> foo5,
	CasADi::CRSSparsity foo6,
	std::vector<CasADi::SX> foo7,
	std::vector< std::vector<CasADi::SX> > foo8,
	CasADi::Matrix<CasADi::SX> foo9,
	CasADi::GenericType foo10,
	int &bar,
	double &baz) {}
};
%}

#ifdef SWIGPYTHON
#ifdef WITH_SWIG_SPLIT
%pythoncode %{
import _casadi_main as _casadi_main_module
%}
#endif // WITH_SWIG_SPLIT
#ifndef WITH_SWIG_SPLIT
%pythoncode %{
_casadi_main_module = _casadi
%}
#endif // WITH_SWIG_SPLIT
#endif // SWIGPYTHON

namespace std {
void dummy(CasADi::SX foo,
	std::vector< std::vector<double> > foo1,
	std::vector<double> &foo2,
	std::vector<CasADi::MX> &foo3,
	CasADi::MX foo4,
	CasADi::Matrix<double> foo5,
	CasADi::CRSSparsity foo6,
	std::vector<CasADi::SX> foo7,
	std::vector< std::vector<CasADi::SX> > foo8,
	CasADi::Matrix<CasADi::SX> foo9,
  CasADi::GenericType foo10,
	int &bar,
	double &baz);
};

