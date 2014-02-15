/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

// Turn off the warnings that certain methods are effectively ignored, this seams to be a false warning, 
// for example vertcat(SXMatrixVector), vertcat(DMatrixVector) and vertcat(MXVector) appears to work fine
#pragma SWIG nowarn=509,303,302

// Incude cmath early on, see #622
%begin %{
#include <cmath>
#ifdef _XOPEN_SOURCE
#undef _XOPEN_SOURCE
#endif
#ifdef _POSIX_C_SOURCE
#undef _POSIX_C_SOURCE
#endif
%}

#ifdef SWIGPYTHON
%pythoncode %{

class _copyableObject(_object):
  def __copy__(self):
    return self.__class__(self)

  def __deepcopy__(self,dummy=None):
    shallow = self.__class__(self)
    if hasattr(self,'makeUnique'):
      shallow.makeUnique()
    return shallow

_object = _copyableObject
%}
#endif // WITH_SWIGPYTHON

%include "doc.i"

%feature("autodoc", "1");

// Make sure that a copy constructor is created
%copyctor;

// STL
#ifdef SWIGXML
//%include <typemaps/std_string.swg>
//%include <std/std_vector.i>
//%include <std/std_pair.i>
namespace std {

template<class T>
class vector {};

template<class A,class B>
class pair {};

}
#else
%include "std_string.i"
%include "std_vector.i"
%include "std_pair.i"
#endif


#ifdef SWIG_MAIN_MODULE
%template(StringVector) std::vector<std::string>;

%template(BVector)             std::vector<bool> ;
%template(BVectorVector)       std::vector<std::vector<bool> > ;
%template(BVectorVectorVector) std::vector< std::vector<std::vector<bool> > > ;

%template(UCharVector)         std::vector<unsigned char>;

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

%template()  std::vector<unsigned char>;

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

%{
#include "symbolic/casadi_options.hpp" 
#include "symbolic/casadi_meta.hpp" 
%}
%include "symbolic/casadi_options.hpp"
%include "symbolic/casadi_meta.hpp"

#ifdef CASADI_MODULE
%{
#define START \
  if (CasADi::CasadiOptions::catch_errors_python){ \
  try {
  
#define STOP \
  } catch (const std::exception& e) { \
  SWIG_exception(SWIG_RuntimeError, e.what()); \
  } catch (const char* e) { \
    SWIG_exception(SWIG_RuntimeError, e); \
  } \
} else
%}

// Exceptions handling
%include "exception.i"
%exception {
  START $action STOP { $action }
}

// Python sometimes takes an approach to not check, but just try.
// It expects a python error to be thrown.
%exception __int__ {
 try {
    $action
  } catch (const std::exception& e) { \
  SWIG_exception(SWIG_RuntimeError, e.what()); \
  } catch (const char* e) { \
    SWIG_exception(SWIG_RuntimeError, e); \
  }
}

// See https://github.com/casadi/casadi/issues/701
// Recent numpys will only catch TypeError or ValueError in printing logic
%exception __nonzero__ {
 try {
    $action
    // foobar
  } catch (const std::exception& e) { \
  SWIG_exception(SWIG_TypeError, e.what()); \
  } catch (const char* e) { \
    SWIG_exception(SWIG_TypeError, e); \
  }
}
#endif // CASADI_MODULE

#ifdef SWIGPYTHON
#ifndef WITH_NUMPY
//#warning "Not using numpy. option(WITH_NUMPY = OFF)"
#endif

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

#define memberbinopsr_custom(Type,uname,custom) \
Type __r##uname##__(const Type& b) const{ return b ##custom## (*$self);}


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
memberbinopsr(Type,truediv) \
memberbinopsr(Type,mldivide) \
memberbinopsr(Type,mrdivide) \
memberbinopsr(Type,mpower) \
memberbinopsr(Type,constpow) \
memberbinopsr_custom(Type,ge,>=) \
memberbinopsr_custom(Type,gt,>) \
memberbinopsr_custom(Type,le,<=) \
memberbinopsr_custom(Type,lt,<) \
memberbinopsr_custom(Type,eq,==) \
memberbinopsr_custom(Type,ne,!=) \
memberbinopsr_un(Type,fmin) \
memberbinopsr_un(Type,fmax) \
memberbinopsr_nn(Type,mul) \
memberbinopsr_un(Type,arctan2) \
memberbinopsr(Type,copysign)

#define memberbinops(uname,argtype,argCast,selfCast,returntype) \
returntype __##uname##__ (argtype) const{ return selfCast(*$self).__##uname##__(argCast(b));} \
returntype __r##uname##__(argtype) const{ return argCast(b).__##uname##__(selfCast(*$self));} \

#define memberbinops_custom(uname,custom,argtype,argCast,selfCast,returntype) \
returntype __##uname##__ (argtype) const{ return selfCast(*$self) ##custom## argCast(b);} \
returntype __r##uname##__(argtype) const{ return argCast(b) ##custom## selfCast(*$self);} \

#define memberbinops_un(uname,argtype,argCast,selfCast,returntype) \
returntype __r##uname##__(argtype) const{ return argCast(b).##uname##(selfCast(*$self));}

// These methods must be added since the implicit type cast does not work.
// Consider a+b  with a DMatrix and b SXMatrix
// In C++, operator+(SXMatrix,SXMatrix) will be called (implicit cast)
// In octave, a.__add__(b) will be called   (no implicit cast)
// In python, __array_priority__ will be checked and b.__radd__(a) will be called (effectively implicit casting)

// This is a list of all operators:
#define binopsFull(argtype,argCast,selfCast,returntype) \
memberbinops_un(fmin,argtype,argCast,selfCast,returntype) \
memberbinops_un(fmax,argtype,argCast,selfCast,returntype) \
memberbinops(constpow,argtype,argCast,selfCast,returntype) \
memberbinops_un(arctan2,argtype,argCast,selfCast,returntype) \
memberbinops(copysign,argtype,argCast,selfCast,returntype) \
memberbinops(pow,argtype,argCast,selfCast,returntype) \
memberbinops(add,argtype,argCast,selfCast,returntype) \
memberbinops(sub,argtype,argCast,selfCast,returntype) \
memberbinops(mul,argtype,argCast,selfCast,returntype) \
memberbinops_custom(ge,>=,argtype,argCast,selfCast,returntype) \
memberbinops_custom(le,<=,argtype,argCast,selfCast,returntype) \
memberbinops_custom(gt,>,argtype,argCast,selfCast,returntype) \
memberbinops_custom(lt,<,argtype,argCast,selfCast,returntype) \
memberbinops_custom(eq,==,argtype,argCast,selfCast,returntype) \
memberbinops_custom(ne,!=,argtype,argCast,selfCast,returntype) \
returntype mul (argtype) const{ return mul(selfCast(*$self) , argCast(b));} \
returntype rmul (argtype) const{ return mul(argCast(b) , selfCast(*$self));} \
memberbinops(div,argtype,argCast,selfCast,returntype) \
memberbinops(truediv,argtype,argCast,selfCast,returntype) \
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

%include "symbolic/fx/schemes_metadata.hpp"

// common typemaps
%include "commontypemaps.i"

%{
#include "symbolic/matrix/matrix.hpp" 
#include "symbolic/matrix/matrix_tools.hpp" 
#include "symbolic/matrix/sparsity_tools.hpp" 
	 
// Scalar expressions 
#include "symbolic/sx/sx.hpp" 
#include "symbolic/sx/sx_tools.hpp" 
#include "symbolic/fx/sx_function.hpp" 
	 
// Matrix expressions 
#include "symbolic/mx/mx.hpp" 
#include "symbolic/mx/mx_tools.hpp" 

#include "symbolic/fx/mx_function.hpp" 
 	
#include "symbolic/fx/mx_function.hpp"
#include "symbolic/fx/custom_function.hpp"
#include "symbolic/fx/ocp_solver.hpp"
#include "symbolic/fx/simulator.hpp"
#include "symbolic/fx/parallelizer.hpp"
#include "symbolic/fx/external_function.hpp"


#include "optimal_control/ocp_tools.hpp"
#include "optimal_control/direct_multiple_shooting.hpp"
#include "optimal_control/symbolic_ocp.hpp"

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
	std::vector < CasADi::Matrix<double> > foo11,
  std::vector < std::vector < CasADi::Matrix<double> > > foo12,
  std::vector < CasADi::Matrix<int> > foo13,
  std::vector < std::vector < CasADi::Matrix<int> > > foo14,
  std::vector < CasADi::Matrix<CasADi::SX> > foo15,
  std::vector < std::vector < CasADi::Matrix<CasADi::SX> > > foo16,
  std::vector < std::vector < CasADi::MX > > foo17,
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
  std::vector < CasADi::Matrix<double> > foo11,
  std::vector < std::vector < CasADi::Matrix<double> > > foo12,
  std::vector < CasADi::Matrix<int> > foo13,
  std::vector < std::vector < CasADi::Matrix<int> > > foo14,
  std::vector < CasADi::Matrix<CasADi::SX> > foo15,
  std::vector < std::vector < CasADi::Matrix<CasADi::SX> > > foo16,
  std::vector < std::vector < CasADi::MX > > foo17,
	int &bar,
	double &baz);
};
