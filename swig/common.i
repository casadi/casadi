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
// for example vertcat(SXVector), vertcat(DMatrixVector) and vertcat(MXVector) appears to work fine
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

%begin %{
#define SWIG_PYTHON_OUTPUT_TUPLE
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

%naturalvar;

// Make sure that a copy constructor is created
%copyctor;

#ifndef SWIGXML
%feature("compactdefaultargs","1");
%feature("compactdefaultargs","0") CasADi::taylor; // taylor function has a default argument for which the namespace is not recognised by SWIG
#endif //SWIGXML

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

template<class A,class B>
class map {};

}
#else

#ifdef SWIGPYTHON
%fragment("StdSequenceTraits","header",
	  fragment="StdTraits",
	  fragment="SwigPySequence_Cont")
{
namespace swig {
  template <class SwigPySeq, class Seq>
  inline void
  assign(const SwigPySeq& swigpyseq, Seq* seq) {
    // seq->assign(swigpyseq.begin(), swigpyseq.end()); // not used as not always implemented
    typedef typename SwigPySeq::value_type value_type;
    typename SwigPySeq::const_iterator it = swigpyseq.begin();
    for (;it != swigpyseq.end(); ++it) {
      seq->insert(seq->end(),(value_type)(*it));
    }
  }

  template <class Seq, class T = typename Seq::value_type >
  struct traits_asptr_stdseq {
    typedef Seq sequence;
    typedef T value_type;

    static int asptr(PyObject *obj, sequence **seq) {
      if (obj == Py_None || SWIG_Python_GetSwigThis(obj)) {
	sequence *p;
	if (::SWIG_ConvertPtr(obj,(void**)&p,
			      swig::type_info<sequence>(),0) == SWIG_OK) {
	  if (seq) *seq = p;
	  return SWIG_OLDOBJ;
	}
      } else if (PySequence_Check(obj)) {
	try {
	  SwigPySequence_Cont<value_type> swigpyseq(obj);
	  if (seq) {
	    sequence *pseq = new sequence();
	    assign(swigpyseq, pseq);
	    *seq = pseq;
	    return SWIG_NEWOBJ;
	  } else {
	    return swigpyseq.check() ? SWIG_OK : SWIG_ERROR;
	  }
	} catch (std::exception& e) {
	  if (seq) {
	    if (!PyErr_Occurred()) {
	      PyErr_SetString(PyExc_TypeError, e.what());
	    }
	  }
	  return SWIG_ERROR;
	}
      }
      return SWIG_ERROR;
    }
  };

  template <class Seq, class T = typename Seq::value_type >
  struct traits_from_stdseq {
    typedef Seq sequence;
    typedef T value_type;
    typedef typename Seq::size_type size_type;
    typedef typename sequence::const_iterator const_iterator;

    static PyObject *from(const sequence& seq) {
%#ifdef SWIG_PYTHON_EXTRA_NATIVE_CONTAINERS
      swig_type_info *desc = swig::type_info<sequence>();
      if (desc && desc->clientdata) {
	return SWIG_NewPointerObj(new sequence(seq), desc, SWIG_POINTER_OWN);
      }
%#endif
      size_type size = seq.size();
      if (size <= (size_type)INT_MAX) {
	PyObject *obj = PyTuple_New((int)size);
	int i = 0;
	for (const_iterator it = seq.begin();
	     it != seq.end(); ++it, ++i) {
	  PyTuple_SetItem(obj,i,swig::from<value_type>(*it));
	}
	return obj;
      } else {
	PyErr_SetString(PyExc_OverflowError,"sequence size not valid in python");
	return NULL;
      }
    }
  };
  
  template <class T >
  struct traits_from_stdseq< std::vector<T> , T > {
    typedef std::vector<T> sequence;
    typedef T value_type;
    typedef typename sequence::size_type size_type;
    typedef typename sequence::const_iterator const_iterator;

    static PyObject *from(const sequence& seq) {
      size_type size = seq.size();
      if (size <= (size_type)INT_MAX) {
	PyObject *obj = PyList_New((int)size);
	int i = 0;
	for (const_iterator it = seq.begin();
	     it != seq.end(); ++it, ++i) {
	  PyList_SetItem(obj,i,swig::from<value_type>(*it));
	}
	return obj;
      } else {
	PyErr_SetString(PyExc_OverflowError,"sequence size not valid in python");
	return NULL;
      }
    }
  };
}
}
%fragment("StdVectorTraits","header",fragment="StdSequenceTraits")
%{
  namespace swig {
    template <class T>
    struct traits_asptr<std::vector<T> >  {
      static int asptr(PyObject *obj, std::vector<T> **vec) {
	return traits_asptr_stdseq<std::vector<T> >::asptr(obj, vec);
      }
    };
    
    template <class T>
    struct traits_from<std::vector<T> > {
      static PyObject *from(const std::vector<T>& vec) {
	return traits_from_stdseq<std::vector<T> >::from(vec);
      }
    };
  }
%}
#endif


%include "std_string.i"

#ifdef SWIGPYTHON

%{
#include <vector>
%}

%include <std/std_common.i>
%include <std_container.i>

namespace std {

  template<class _Tp, class _Alloc = allocator< _Tp > >
  class vector {
  public:
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    typedef _Tp value_type;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef _Tp& reference;
    typedef const _Tp& const_reference;
    typedef _Alloc allocator_type;

    %traits_swigtype(_Tp);

    %fragment(SWIG_Traits_frag(std::vector<_Tp, _Alloc >), "header",
	      fragment=SWIG_Traits_frag(_Tp),
	      fragment="StdVectorTraits") {
      namespace swig {
	template <>  struct traits<std::vector<_Tp, _Alloc > > {
	  typedef pointer_category category;
	  static const char* type_name() {
	    return "std::vector<" #_Tp "," #_Alloc " >";
	  }
	};
      }
    }

    %typemap_traits_ptr(SWIG_TYPECHECK_VECTOR, std::vector<_Tp, _Alloc >);
  
  };

}
#else
%include "std_vector.i"
#endif
%include "std_pair.i"
%include "std_map.i"
#endif

//#ifdef SWIG_MAIN_MODULE
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
//#endif //SWIG_MAIN_MODULE
//#ifndef SWIG_MAIN_MODULE
/**%template() std::vector<std::string>;

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
%template() std::vector< std::pair<int,int> >;*/
//#endif //SWIG_MAIN_MODULE


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

#ifdef SWIGPYTHON
%pythoncode %{
  CasadiOptions.setDeprecatedWarningAsException(True)
%}
#endif

#ifdef CASADI_MODULE

#ifdef SWIGPYTHON
%{
#define START \
  if (CasADi::CasadiOptions::catch_errors_python){ \
  try {
  
#define STOP \
  } catch (const CasADi::CasadiDeprecationException& e) { \
    PyErr_WarnEx(PyExc_SyntaxWarning,e.what(),2); \
  } catch (const std::exception& e) { \
  SWIG_exception(SWIG_RuntimeError, e.what()); \
  } catch (const char* e) { \
    SWIG_exception(SWIG_RuntimeError, e); \
  } \
} else
%}
#endif

#ifdef SWIGOCTAVE
%{
#define START \
  if (CasADi::CasadiOptions::catch_errors_python){ \
  try {
  
#define STOP \
  } catch (const CasADi::CasadiDeprecationException& e) { \
    warning(e.what()); \
  } catch (const std::exception& e) { \
  SWIG_exception(SWIG_RuntimeError, e.what()); \
  } catch (const char* e) { \
    SWIG_exception(SWIG_RuntimeError, e); \
  } \
} else
%}
#endif

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
// Consider a+b  with a DMatrix and b SX
// In C++, operator+(SX,SX) will be called (implicit cast)
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

//%traits_swigtype(CasADi::GenericType);
//%traits_swigtype(std::mapCasADi::Dictionary);

//%fragment(SWIG_Traits_frag(std::map< std::string, CasADi::GenericType, std::less<std::string > , allocator<std::pair<const std::string, CasADi::GenericType > > >));

//%fragment(SWIG_Traits_frag(CasADi::Dictionary));

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
#include "symbolic/sx/sx_element.hpp" 
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


#include "optimal_control/direct_multiple_shooting.hpp"
#include "optimal_control/symbolic_ocp.hpp"

#include "symbolic/fx/io_interface.hpp"
#include "symbolic/fx/fx.hpp"
#include "symbolic/fx/sx_function.hpp"
#include "symbolic/fx/mx_function.hpp"
#include "symbolic/fx/linear_solver.hpp"
#include "symbolic/fx/symbolic_qr.hpp"
#include "symbolic/fx/implicit_function.hpp"
#include "symbolic/fx/integrator.hpp"
#include "symbolic/fx/simulator.hpp"
#include "symbolic/fx/control_simulator.hpp"
#include "symbolic/fx/nlp_solver.hpp"
#include "symbolic/fx/qp_solver.hpp"
#include "symbolic/fx/stabilized_qp_solver.hpp"
#include "symbolic/fx/lp_solver.hpp"
#include "symbolic/fx/ocp_solver.hpp"
#include "symbolic/fx/sdp_solver.hpp"
#include "symbolic/fx/socp_solver.hpp"
#include "symbolic/fx/qcqp_solver.hpp"
#include "symbolic/fx/sdqp_solver.hpp"
#include "symbolic/fx/external_function.hpp"
#include "symbolic/fx/parallelizer.hpp"
#include "symbolic/fx/custom_function.hpp"
#include "symbolic/fx/nullspace.hpp"

#include "nonlinear_programming/symbolic_nlp.hpp"
#include "nonlinear_programming/sqp_method.hpp"
#include "nonlinear_programming/stabilized_sqp_method.hpp"
#include "nonlinear_programming/scpgen.hpp"
#include "nonlinear_programming/nlp_qp_solver.hpp"
#include "nonlinear_programming/nlp_implicit_solver.hpp"
#include "nonlinear_programming/newton_implicit_solver.hpp"

#include "integration/fixed_step_integrator.hpp"
#include "integration/implicit_fixed_step_integrator.hpp"
#include "integration/rk_integrator.hpp"
#include "integration/collocation_integrator.hpp"
#include "integration/old_collocation_integrator.hpp"
#include "integration/integration_tools.hpp"


#include "optimal_control/variable.hpp"
#include "optimal_control/symbolic_ocp.hpp"
#include "optimal_control/direct_single_shooting.hpp"
#include "optimal_control/direct_multiple_shooting.hpp"
#include "optimal_control/direct_collocation.hpp"

#include "control/dple_solver.hpp"
#include "control/simple_indef_dple_solver.hpp"

%}

#ifndef SWIGXML
%traits_swigtype(CasADi::DerivativeGenerator);
%fragment(SWIG_Traits_frag(CasADi::DerivativeGenerator));
%traits_swigtype(CasADi::Callback);
%fragment(SWIG_Traits_frag(CasADi::Callback));
%traits_swigtype(CasADi::CustomEvaluate);
%fragment(SWIG_Traits_frag(CasADi::CustomEvaluate));
%traits_swigtype(CasADi::IndexList);
%fragment(SWIG_Traits_frag(CasADi::IndexList));

%template(Dictionary) std::map<std::string,CasADi::GenericType>;

%traits_swigtype(CasADi::FX);
%fragment(SWIG_Traits_frag(CasADi::FX));

#endif


// These dummy things would go away when we properly use fragments
// %traits_swigtype

%{
namespace std {
void dummy(CasADi::SXElement foo,
	std::vector< std::vector<double> > foo1,
	std::vector<double> &foo2,
	std::vector<CasADi::MX> &foo3,
	CasADi::MX foo4,
	CasADi::Matrix<double> foo5,
	CasADi::Sparsity foo6,
	std::vector<CasADi::SXElement> foo7,
	std::vector< std::vector<CasADi::SXElement> > foo8,
	CasADi::Matrix<CasADi::SXElement> foo9,
	CasADi::GenericType foo10,
	std::vector < CasADi::Matrix<double> > foo11,
  std::vector < std::vector < CasADi::Matrix<double> > > foo12,
  std::vector < CasADi::Matrix<int> > foo13,
  std::vector < std::vector < CasADi::Matrix<int> > > foo14,
  std::vector < CasADi::Matrix<CasADi::SXElement> > foo15,
  std::vector < std::vector < CasADi::Matrix<CasADi::SXElement> > > foo16,
  std::vector < std::vector < CasADi::MX > > foo17,
  std::vector < std::vector < CasADi::MX* > > foo17b,
  CasADi::Dictionary foo18,
  std::string& foo19,
  CasADi::Matrix<int> foo20,
  CasADi::CustomFunction foo24,
  CasADi::FX foo25,
	int &bar,
	double &baz) {}


#ifdef SWIGPYTHON
void dummy2(
  CasADi::DerivativeGenerator foo1,
  CasADi::Callback foo2,
  CasADi::CustomEvaluate foo3
  ) {}
#endif// SWIGPYTHON
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
void dummy(CasADi::SXElement foo,
	std::vector< std::vector<double> > foo1,
	std::vector<double> &foo2,
	std::vector<CasADi::MX> &foo3,
	CasADi::MX foo4,
	CasADi::Matrix<double> foo5,
	CasADi::Sparsity foo6,
	std::vector<CasADi::SXElement> foo7,
	std::vector< std::vector<CasADi::SXElement> > foo8,
	CasADi::Matrix<CasADi::SXElement> foo9,
  CasADi::GenericType foo10,
  std::vector < CasADi::Matrix<double> > foo11,
  std::vector < std::vector < CasADi::Matrix<double> > > foo12,
  std::vector < CasADi::Matrix<int> > foo13,
  std::vector < std::vector < CasADi::Matrix<int> > > foo14,
  std::vector < CasADi::Matrix<CasADi::SXElement> > foo15,
  std::vector < std::vector < CasADi::Matrix<CasADi::SXElement> > > foo16,
  std::vector < std::vector < CasADi::MX > > foo17,
  std::vector < std::vector < CasADi::MX* > > foo17b,
  CasADi::Dictionary foo18,
  std::string& foo19,
  CasADi::Matrix<int> foo20,
  CasADi::CustomFunction foo24,
  CasADi::FX foo25,
	int &bar,
	double &baz);


#ifdef SWIGPYTHON
void dummy2(
  CasADi::DerivativeGenerator foo1,
  CasADi::Callback foo2,
  CasADi::CustomEvaluate foo3
  ) {}
#endif

};
