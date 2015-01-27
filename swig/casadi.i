/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


%module(package="casadi") casadi

 // Include all public CasADi C++
%{
#include <casadi/casadi.hpp>
%}

// Turn off the warnings that certain methods are effectively ignored, this seams to be a false warning, 
// for example vertcat(SXVector), vertcat(DMatrixVector) and vertcat(MXVector) appears to work fine
#pragma SWIG nowarn=509,303,302

#define CASADI_EXPORT

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

%ignore *::operator->;

#ifdef SWIGMATLAB
%rename(disp) repr;
#else
%ignore print;
%ignore repr;
#endif

%begin %{
#define SWIG_PYTHON_OUTPUT_TUPLE
%}

// Print representation
#ifdef SWIGMATLAB
#define SWIG_REPR disp
#else
#define SWIG_REPR __repr__
#endif

// Print description
#ifdef SWIGMATLAB
#define SWIG_STR print
#else
#define SWIG_STR __str__
#endif


//#endif // SWIGPYTHON


#ifdef SWIGPYTHON
%pythoncode %{

import contextlib

@contextlib.contextmanager
def internalAPI():
    backup = CasadiOptions.getAllowedInternalAPI()
    CasadiOptions.setAllowedInternalAPI(True)
    try:
      yield
    finally:
      CasadiOptions.setAllowedInternalAPI(backup)

class _copyableObject(_object):
  def __copy__(self):
    return self.__class__(self)

  def __deepcopy__(self,dummy=None):
    shallow = self.__class__(self)
    if hasattr(self,'makeUnique'):
      shallow.makeUnique()
    return shallow

_object = _copyableObject

_swig_repr_default = _swig_repr
def _swig_repr(self):
  if hasattr(self,'getRepresentation'):
    return self.getRepresentation()
  else:
    return _swig_repr_default(self)

%}
#endif // WITH_SWIGPYTHON

#ifdef SWIGPYTHON
%include "doc_merged.i"
#else
%include "doc.i"
#endif


%feature("autodoc", "1");

%naturalvar;

// Make sure that a copy constructor is created
%copyctor;

#ifndef SWIGXML
%feature("compactdefaultargs","1");
//%feature("compactdefaultargs","0") casadi::taylor; // taylor function has a default argument for which the namespace is not recognised by SWIG
%feature("compactdefaultargs","0") casadi::solve; // buggy
%feature("compactdefaultargs","0") casadi::Function::generateCode; // buggy
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

#ifndef SWIGMATLAB
%template(Pair_Int_Int) std::pair<int,int>;
%template(VectorPair_Int_Int) std::vector< std::pair<int,int> >;
#endif // SWIGMATLAB


// The following is a work-around since it appears not possible to use the standard print functions from stl_vector tools,
// nor the std::stringstream class, since these are included _after_ std::vector in the C++ generated wrapper code
%extend std::vector<double>{  
  std::string SWIG_REPR(){
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
  std::string SWIG_STR(){
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
  std::string SWIG_REPR(){
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
  std::string SWIG_STR(){
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

%define DEPRECATED_MSG(MSG)
if (deprecated("$decl",MSG)) SWIG_fail;
%enddef

%define INTERNAL_MSG()
if (internal("$decl")) SWIG_fail;
%enddef

#ifdef SWIGPYTHON
%wrapper %{
int deprecated(const std::string & c,const std::string & a) {
  std::string msg = "This CasADi function (" + c + ") is deprecated. " + a;
  return PyErr_WarnEx(PyExc_DeprecationWarning,msg.c_str(),2);
}
int internal(const std::string & c) {
  if (CasadiOptions::allowed_internal_api) return 0;
  std::string msg = "This CasADi function (" + c + ") is not part of the public API. Use at your own risk.";
  return PyErr_WarnEx(PyExc_SyntaxWarning,msg.c_str(),2);
}
%}
#endif // SWIGPYTHON

#ifdef SWIGMATLAB
%wrapper %{
int deprecated(const std::string & c,const std::string & a) {
  std::string msg = "This CasADi function (" + c + ") is deprecated. " + a;
  mexWarnMsgIdAndTxt("SWIG:DeprecationWarning",msg.c_str());
  return 1;
}
int internal(const std::string & c) {
  if (CasadiOptions::allowed_internal_api) return 0;
  std::string msg = "This CasADi function (" + c + ") is not part of the public API. Use at your own risk.";
  mexWarnMsgIdAndTxt("SWIG:SyntaxWarning",msg.c_str());
  return 1;
}
%}
#endif // SWIGMATLAB

#ifndef SWIGXML
%{
// TODO(jgillis): remove after internal.i was updated
#define CATCH_OR_RETHROW \
  catch (const std::exception& e) { \
    if (casadi::CasadiOptions::catch_errors_swig) { \
      SWIG_exception(SWIG_RuntimeError, e.what()); \
    } else { \
      throw e; \
    } \
  }
#define CATCH_OR_NOT(...) \
if (casadi::CasadiOptions::catch_errors_swig) { \
  try { \
    __VA_ARGS__ \
  } catch(const std::exception& e) { \
    SWIG_exception(SWIG_RuntimeError, e.what()); \
  } \
} else { \
  __VA_ARGS__ \
}

%}
#endif

// Exceptions handling
%include "exception.i"
%exception {
  CATCH_OR_NOT($action)
}

// Python sometimes takes an approach to not check, but just try.
// It expects a python error to be thrown.
%exception __int__ {
 try {
    $action
  } catch (const std::exception& e) { \
  SWIG_exception(SWIG_RuntimeError, e.what()); \
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
  }
}
%include "internal.i"
%include "deprecated.i"

#ifdef SWIGPYTHON

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
#endif // SWIGPYTHON


#define memberbinopsr(Type,uname) \
Type __r##uname##__(const Type& b) const{ return b.__##uname##__(*$self);}

#define memberbinopsr_custom(Type,uname,custom) \
Type __r##uname##__(const Type& b) const{ return b ##custom## (*$self);}


#define memberbinopsr_un(Type,uname) \
Type __r##uname##__(const Type& b) const{ return b.##uname##(*$self);}

#define memberbinopsr_nn(Type,uname) \
Type r##uname##(const Type& b) const{ return b.##uname##(*$self);}

%define binopsrFull(Type)
Type __rpow__(const Type& b) const{ return pow(b, *$self);}
Type __radd__(const Type& b) const{ return b + *$self;}
Type __rsub__(const Type& b) const{ return b - *$self;}
Type __rmul__(const Type& b) const{ return b * *$self;}
Type __rdiv__(const Type& b) const{ return b / *$self;}
memberbinopsr(Type,truediv)
memberbinopsr(Type,mldivide)
memberbinopsr(Type,mrdivide)
Type __rmpower__(const Type& b) const{ return b.zz_mpower(*$self);}
memberbinopsr(Type,constpow)
Type __rge__(const Type& b) const{ return b >= (*$self);}
Type __rgt__(const Type& b) const{ return b > (*$self);}
Type __rle__(const Type& b) const{ return b <= (*$self);}
Type __rlt__(const Type& b) const{ return b < (*$self);}
Type __req__(const Type& b) const{ return b == (*$self);}
Type __rne__(const Type& b) const{ return b != (*$self);}
Type __rfmin__(const Type& b) const { return fmin(b, *$self);}
Type __rfmax__(const Type& b) const { return fmax(b, *$self);}
Type rmul(const Type& b) const{ return b.zz_mtimes(*$self);}
Type __rarctan2__(const Type& b) const{ return b.zz_atan2(*$self);}
memberbinopsr(Type,copysign)
%enddef

%define memberbinops(uname,argtype,argCast,selfCast,returntype)
returntype __##uname##__ (argtype) const{ return selfCast(*$self).__##uname##__(argCast(b));}
returntype __r##uname##__(argtype) const{ return argCast(b).__##uname##__(selfCast(*$self));}
%enddef

%define memberbinops_custom(uname,custom,argtype,argCast,selfCast,returntype)
returntype __##uname##__ (argtype) const{ return selfCast(*$self) ##custom## argCast(b);}
returntype __r##uname##__(argtype) const{ return argCast(b) ##custom## selfCast(*$self);}
%enddef

%define memberbinops_un(uname,argtype,argCast,selfCast,returntype)
returntype __r##uname##__(argtype) const{ return argCast(b).##uname##(selfCast(*$self));}
%enddef

// These methods must be added since the implicit type cast does not work.
// Consider a+b  with a DMatrix and b SX
// In C++, operator+(SX,SX) will be called (implicit cast)
// In python, __array_priority__ will be checked and b.__radd__(a) will be called (effectively implicit casting)

// This is a list of all operators:
%define binopsFull(argtype,argCast,selfCast,returntype)
returntype __rfmin__(argtype) const { return fmin(argCast(b), selfCast(*$self));}
returntype __rfmax__(argtype) const { return fmax(argCast(b), selfCast(*$self));}
memberbinops(constpow,argtype,argCast,selfCast,returntype)
returntype __rarctan2__(argtype) const{ return argCast(b).zz_atan2(selfCast(*$self));}
memberbinops(copysign,argtype,argCast,selfCast,returntype)
returntype __pow__ (argtype) const { return selfCast(*$self).zz_power(argCast(b));}
returntype __rpow__(argtype) const { return argCast(b).zz_power(selfCast(*$self));}
returntype __add__ (argtype) const { return selfCast(*$self).zz_plus(argCast(b));}
returntype __radd__(argtype) const { return argCast(b).zz_plus(selfCast(*$self));}
returntype __sub__ (argtype) const { return selfCast(*$self).zz_minus(argCast(b));}
returntype __rsub__(argtype) const { return argCast(b).zz_minus(selfCast(*$self));}
returntype __mul__ (argtype) const { return selfCast(*$self).zz_times(argCast(b));}
returntype __rmul__(argtype) const { return argCast(b).zz_times(selfCast(*$self));}
returntype __div__ (argtype) const { return selfCast(*$self).zz_rdivide(argCast(b));}
returntype __rdiv__(argtype) const { return argCast(b).zz_rdivide(selfCast(*$self));}
memberbinops_custom(ge,>=,argtype,argCast,selfCast,returntype)
memberbinops_custom(le,<=,argtype,argCast,selfCast,returntype)
memberbinops_custom(gt,>,argtype,argCast,selfCast,returntype)
memberbinops_custom(lt,<,argtype,argCast,selfCast,returntype)
memberbinops_custom(eq,==,argtype,argCast,selfCast,returntype)
memberbinops_custom(ne,!=,argtype,argCast,selfCast,returntype)
returntype mul (argtype) const{ return mul(selfCast(*$self) , argCast(b));}
returntype rmul (argtype) const{ return mul(argCast(b) , selfCast(*$self));}
memberbinops(truediv,argtype,argCast,selfCast,returntype)
memberbinops(mldivide,argtype,argCast,selfCast,returntype)
memberbinops(mrdivide,argtype,argCast,selfCast,returntype)
returntype __mpower__ (argtype) const{ return selfCast(*$self).zz_mpower(argCast(b));}
returntype __rmpower__(argtype) const{ return argCast(b).zz_mpower(selfCast(*$self));}
%enddef

// This is a list of operators that do not check __array_priority__ in python
%define binopsNoPriority(argtype,argCast,selfCast,returntype)
returntype __pow__ (argtype) const { return pow(selfCast(*$self), argCast(b));}
returntype __rpow__(argtype) const { return pow(argCast(b), selfCast(*$self));}
%enddef

// typemaphelpers
%include "typemaphelpers.i"

#ifdef SWIGPYTHON
%include "python/meta_python.i"
#endif

#ifdef SWIGMATLAB
%include "matlab/meta_matlab.i"
#endif

// common typemaps
%include "commontypemaps.i"

%{
using namespace casadi;
%}

#ifdef SWIGMATLAB
%typemap(out) std::pair<int, int> {
  $result = mxCreateDoubleMatrix(1, 2, mxREAL);
  double* data = static_cast<double*>(mxGetData($result));
  data[0] = $1.first;
  data[1] = $1.second;
}

// The natural way to represent an integer in MATLAB is as a double
%typemap(out) int {
  $result = mxCreateDoubleScalar($1);
}
#endif // SWIGMATLAB

#ifndef SWIGXML
%traits_swigtype(casadi::DerivativeGenerator);
%fragment(SWIG_Traits_frag(casadi::DerivativeGenerator));
%traits_swigtype(casadi::Callback);
%fragment(SWIG_Traits_frag(casadi::Callback));
%traits_swigtype(casadi::CustomEvaluate);
%fragment(SWIG_Traits_frag(casadi::CustomEvaluate));

%template(Dictionary) std::map<std::string,casadi::GenericType>;

%traits_swigtype(casadi::Function);
%fragment(SWIG_Traits_frag(casadi::Function));

#endif


// These dummy things would go away when we properly use fragments
// %traits_swigtype

%{
namespace std {
void dummy(casadi::SXElement foo,
           std::vector< std::vector<double> > foo1,
           std::vector<double> &foo2,
           std::vector<casadi::MX> &foo3,
           casadi::MX foo4,
           casadi::Matrix<double> foo5,
           casadi::Sparsity foo6,
           std::vector<casadi::SXElement> foo7,
           std::vector< std::vector<casadi::SXElement> > foo8,
           casadi::SX foo9,
           casadi::GenericType foo10,
           std::vector < casadi::Matrix<double> > foo11,
           std::vector < std::vector < casadi::Matrix<double> > > foo12,
           std::vector < casadi::Matrix<int> > foo13,
           std::vector < std::vector < casadi::Matrix<int> > > foo14,
           std::vector < casadi::SX > foo15,
           std::vector < std::vector < casadi::SX > > foo16,
           std::vector < std::vector < casadi::MX > > foo17,
           std::vector < std::vector < casadi::MX* > > foo17b,
           casadi::Dictionary foo18,
           std::string& foo19,
           casadi::Matrix<int> foo20,
           casadi::CustomFunction foo24,
           casadi::Function foo25,
           int &bar,
           double &baz) {}


#ifdef SWIGPYTHON
void dummy2(
  casadi::DerivativeGenerator foo1,
  casadi::Callback foo2,
  casadi::CustomEvaluate foo3
  ) {}
#endif// SWIGPYTHON
};
%}

#ifdef SWIGPYTHON
%pythoncode %{
import _casadi
%}
#endif // SWIGPYTHON

namespace std {
  void dummy(casadi::SXElement foo,
             std::vector< std::vector<double> > foo1,
             std::vector<double> &foo2,
             std::vector<casadi::MX> &foo3,
             casadi::MX foo4,
             casadi::Matrix<double> foo5,
             casadi::Sparsity foo6,
             std::vector<casadi::SXElement> foo7,
             std::vector< std::vector<casadi::SXElement> > foo8,
             casadi::SX foo9,
             casadi::GenericType foo10,
             std::vector < casadi::Matrix<double> > foo11,
             std::vector < std::vector < casadi::Matrix<double> > > foo12,
             std::vector < casadi::Matrix<int> > foo13,
             std::vector < std::vector < casadi::Matrix<int> > > foo14,
             std::vector < casadi::SX > foo15,
             std::vector < std::vector < casadi::SX > > foo16,
             std::vector < std::vector < casadi::MX > > foo17,
             std::vector < std::vector < casadi::MX* > > foo17b,
             casadi::Dictionary foo18,
             std::string& foo19,
             casadi::Matrix<int> foo20,
             casadi::CustomFunction foo24,
             casadi::Function foo25,
             int &bar,
             double &baz);


#ifdef SWIGPYTHON
void dummy2(
  casadi::DerivativeGenerator foo1,
  casadi::Callback foo2,
  casadi::CustomEvaluate foo3
  ) {}
#endif

};

%include <casadi/core/function/schemes_metadata.hpp>

// Init hooks
#ifdef SWIGPYTHON
#ifdef WITH_PYTHON_INTERRUPTS
%{
#include <pythonrun.h>
void SigIntHandler(int) {
  std::cerr << "Keyboard Interrupt" << std::endl;
  signal(SIGINT, SIG_DFL);
  kill(getpid(), SIGINT);
}
%}

%init %{
PyOS_setsig(SIGINT, SigIntHandler);
%}
#endif // WITH_PYTHON_INTERRUPTS

%pythoncode%{
try:
  from numpy import pi, inf
except:
  pass

try:
  from numpy import sin, cos, tan, sqrt, log, exp, floor, ceil, fmod, fmin, fmax, sinh, cosh, tanh, arcsin, arccos, arctan, arctan2, fabs, sign, arctanh, arcsinh, arccosh, copysign
except:
  sin = lambda x: x.sin()
  cos = lambda x: x.cos()
  tan = lambda x: x.tan()
  arcsin = lambda x: x.arcsin()
  arccos = lambda x: x.arccos()
  arctan = lambda x: x.arctan()
  sqrt = lambda x: x.sqrt()
  log = lambda x: x.log()
  exp = lambda x: x.exp()
  floor = lambda x: x.floor()
  ceil = lambda x: x.ceil()
  fmin = lambda x,y: x.fmin(y)
  fmax = lambda x,y: x.fmax(y)
  sinh = lambda x: x.sinh()
  cosh = lambda x: x.cosh()
  tanh = lambda x: x.tanh()
  fabs = lambda x: x.fabs()
  sign = lambda x: x.sign()
  arctan2 = lambda x,y: x.arctan2(y)
  arctanh = lambda x: x.arctanh()
  arcsinh = lambda x: x.arcsinh()
  arccosh = lambda x: x.arccosh()
  copysign = lambda x,y: x.copysign(y)
%}
#endif // SWIGPYTHON

%rename("%(regex:/zz_(?!ML)(.*)/\\1/)s") ""; // Strip leading zz_ unless followed by ML
#ifdef SWIGPYTHON
%rename(__add__) zz_plus;
%rename(__sub__) zz_minus;
%rename(__mul__) zz_times;
%rename(__div__) zz_rdivide;
%rename(arcsin) zz_asin;
%rename(arccos) zz_acos;
%rename(arctan) zz_atan;
%rename(arctan2) zz_atan2;
%rename(arcsinh) zz_asinh;
%rename(arccosh) zz_acosh;
%rename(arctanh) zz_atanh;
%rename(fabs) zz_abs;
%rename(fmod) zz_mod;
%rename(fmin) zz_min;
%rename(fmax) zz_max;
%rename(mul) zz_mtimes;
%ignore T;
%ignore shape;
%rename(__lt__) zz_lt;
%rename(__gt__) zz_gt;
%rename(__le__) zz_le;
%rename(__ge__) zz_ge;
%rename(__ne__) zz_ne;
%rename(__eq__) zz_eq;
%rename(logic_and) zz_and;
%rename(logic_or) zz_or;
%rename(logic_not) zz_not;
%rename(logic_all) zz_all;
%rename(logic_any) zz_any;
%rename(__pow__) zz_power;
#endif // SWIGPYTHON

#ifdef SWIGMATLAB
%rename(uminus) operator-;
%rename(uplus) operator+;
%rename(ldivide) __rdiv__;
%rename(mrdivide) __mrdivide__;
%rename(mldivide) __mldivide__;
%rename(transpose) T;
%ignore size;
%rename(size) shape;
%feature("varargin","1") zz_vertcat;
%feature("varargin","1") zz_horzcat;
%feature("nonstatic","1") zz_vertcat;
%feature("nonstatic","1") zz_horzcat;

// Explicit type conversion of the first argument of const member functions i.e. this/self
%feature("convertself","1");

// Workarounds, pending proper fix
%rename(truediv) __truediv__;
%rename(nonzero) __nonzero__;
%rename(constpow) __constpow__;
%rename(copysign) __copysign__;
%rename(rcopysign) __rcopysign__;
%rename(hash) __hash__;
#endif // SWIGMATLAB

#ifdef SWIGPYTHON
%pythoncode %{
def prod(self,*args):
    raise Exception("'prod' is not supported anymore in CasADi. Use 'mul' to do matrix multiplication.")
def dot(self,*args):
    raise Exception("'dot' is not supported anymore in CasADi. Use 'mul' to do matrix multiplication.")
    
class NZproxy:
  def __init__(self,matrix):
    self.matrix = matrix
    
  def __getitem__(self,s):
    return self.matrix.getNZ(False, s)

  def __setitem__(self,s,val):
    return self.matrix.setNZ(val, False, s)

  def __len__(self):
    return self.matrix.nnz()
    
  def __iter__(self):
    for k in range(len(self)):
      yield self[k]
%}

%define %matrix_convertors
%pythoncode %{
        
    def toMatrix(self):
        import numpy as n
        return n.matrix(self.toArray())

    def __iter__(self):
      for k in self.nz:
        yield k
        
%}
%enddef 
%define %matrix_helpers(Type)
%pythoncode %{
    @property
    def shape(self):
        return (self.size1(),self.size2())
        
    def reshape(self,arg):
        return _casadi.reshape(self,arg)
        
    @property
    def T(self):
        return _casadi.transpose(self)
        
    def __getitem__(self, s):
        with internalAPI():
          if isinstance(s, tuple) and len(s)==2:
            return self.getSub(False, s[0], s[1])
          return self.getSub(False, s)

    def __setitem__(self,s,val):
        with internalAPI():
          if isinstance(s,tuple) and len(s)==2:
            return self.setSub(val, False, s[0], s[1])  
          return self.setSub(val, False, s)
        
    @property
    def nz(self):
      return NZproxy(self)
        
    def prod(self,*args):
        raise Exception("'prod' is not supported anymore in CasADi. Use 'mul' to do matrix multiplication.")
     
%}
%enddef 
    
%define %python_array_wrappers(arraypriority)
%pythoncode %{

  __array_priority__ = arraypriority

  def __array_wrap__(self,out_arr,context=None):
    if context is None:
      return out_arr
    name = context[0].__name__
    args = list(context[1])
    
    if len(context[1])==3:
      raise Exception("Error with %s. Looks like you are using an assignment operator, such as 'a+=b' where 'a' is a numpy type. This is not supported, and cannot be supported without changing numpy." % name)

    if "vectorized" in name:
        name = name[:-len(" (vectorized)")]
    
    conversion = {"multiply": "mul", "divide": "div", "true_divide": "div", "subtract":"sub","power":"pow","greater_equal":"ge","less_equal": "le", "less": "lt", "greater": "gt"}
    if name in conversion:
      name = conversion[name]
    if len(context[1])==2 and context[1][1] is self and not(context[1][0] is self):
      name = 'r' + name
      args.reverse()
    if not(hasattr(self,name)) or ('mul' in name):
      name = '__' + name + '__'
    fun=getattr(self, name)
    return fun(*args[1:])
     
     
  def __array__(self,*args,**kwargs):
    import numpy as n
    if len(args) > 1 and isinstance(args[1],tuple) and isinstance(args[1][0],n.ufunc) and isinstance(args[1][0],n.ufunc) and len(args[1])>1 and args[1][0].nin==len(args[1][1]):
      if len(args[1][1])==3:
        raise Exception("Error with %s. Looks like you are using an assignment operator, such as 'a+=b'. This is not supported when 'a' is a numpy type, and cannot be supported without changing numpy itself. Either upgrade a to a CasADi type first, or use 'a = a + b'. " % args[1][0].__name__)
      return n.array([n.nan])
    else:
      if hasattr(self,'__array_custom__'):
        return self.__array_custom__(*args,**kwargs)
      else:
        return self.toArray()
      
%}
%enddef
#endif // SWIGPYTHON

#ifdef SWIGXML
%define %matrix_helpers(Type)
%enddef
#endif

#ifdef SWIGMATLAB
%define %matrix_helpers(Type)
    // Get a submatrix (index-1)
    const Type getitem(const Slice& rr) const { return $self->getSub(true, rr);}
    const Type getitem(const Matrix<int>& rr) const { return $self->getSub(true, rr);}
    const Type getitem(const Sparsity& sp) const { return $self->getSub(true, sp);}
    const Type getitem(const Slice& rr, const Slice& cc) const { return $self->getSub(true, rr, cc);}
    const Type getitem(const Slice& rr, const Matrix<int>& cc) const { return $self->getSub(true, rr, cc);}
    const Type getitem(const Matrix<int>& rr, const Slice& cc) const { return $self->getSub(true, rr, cc);}
    const Type getitem(const Matrix<int>& rr, const Matrix<int>& cc) const { return $self->getSub(true, rr, cc);}

    // Set a submatrix (index-1)
    void setitem(const Type& m, const Slice& rr) { $self->setSub(m, true, rr);}
    void setitem(const Type& m, const Matrix<int>& rr) { $self->setSub(m, true, rr);}
    void setitem(const Type& m, const Sparsity& sp) { $self->setSub(m, true, sp);}
    void setitem(const Type& m, const Slice& rr, const Slice& cc) { $self->setSub(m, true, rr, cc);}
    void setitem(const Type& m, const Slice& rr, const Matrix<int>& cc) { $self->setSub(m, true, rr, cc);}
    void setitem(const Type& m, const Matrix<int>& rr, const Slice& cc) { $self->setSub(m, true, rr, cc);}
    void setitem(const Type& m, const Matrix<int>& rr, const Matrix<int>& cc) { $self->setSub(m, true, rr, cc);}

    // Get nonzeros (index-1)
    const Type getitemcurl(const Slice& rr) const { return $self->getNZ(true, rr);}
    const Type getitemcurl(const Matrix<int>& rr) const { return $self->getNZ(true, rr);}

    // Set nonzeros (index-1)
    void setitemcurl(const Type& m, const Slice& rr) { $self->setNZ(m, true, rr);}
    void setitemcurl(const Type& m, const Matrix<int>& rr) { $self->setNZ(m, true, rr);}

    // 'end' function (needed for end syntax in MATLAB)
    inline int end(int i, int n) const {
      return n==1 ? $self->numel() : i==1 ? $self->size1() : $self->size2();
    }

    // Transpose using the A' syntax in addition to A.'
    Type ctranspose() const { return $self->T();}

%enddef
#endif

#ifndef SWIGPYTHON
%define %matrix_convertors
%enddef
#endif

#ifndef SWIGXML
%include "typemaps.i"
#endif

%include "std_vector_tools.i"
%include "weak_ref.i"
%include "options_functionality.i"
%include "casadi_calculus.i"
%include "sparsity.i"
%include "slice.i"
%include "generic_expression.i"
%include "generic_matrix.i"
%include "matrix.i"
%include "sx_element.i"
%include "mx.i"
%include "matrix_tools.i"
%include "generic_expression_tools.i"
%include "sx_tools.i"
%include "mx_tools.i"

#ifdef SWIGPYTHON
%pythoncode %{
def attach_return_type(f,t):
  if not(hasattr(f,'func_annotations')):
    f.func_annotations = {}
  if not(isinstance(getattr(f,'func_annotations'),dict)):
    raise Exception("Cannot annotate this python Method to be a sparsitygenerator. Method has func_annotations attribute with unknown type.")
  f.func_annotations["return"] = t
  return f

def pyderivativegenerator(f):
  return attach_return_type(f,Function)

def pyevaluate(f):
  return attach_return_type(f,None)
  
def pycallback(f):
  return attach_return_type(f,int)
  

def pyfunction(inputs,outputs):
  def wrap(f):
    
    @pyevaluate
    def fcustom(f2):
      res = f([f2.getInput(i) for i in range(f2.getNumInputs())])
      if not isinstance(res,list):
        res = [res]
      for i in range(f2.getNumOutputs()):
        f2.setOutput(res[i],i)
    Fun = CustomFunction(fcustom,inputs,outputs)
    Fun.setOption("name","CustomFunction")
    return Fun
  return wrap
  
def PyFunction(obj,inputs,outputs):
    @pyevaluate
    def fcustom(f):
      obj.evaluate([f.input(i) for i in range(f.getNumInputs())],[f.output(i) for i in range(f.getNumOutputs())])
      
    Fun = CustomFunction(fcustom,inputs,outputs)
    Fun.setOption("name","CustomFunction")
    if hasattr(obj,'getDerivative'):
      @pyderivativegenerator
      def derivativewrap(f,nfwd,nadj):
        return obj.getDerivative(f,nfwd,nadj)
      Fun.setOption("derivative_generator",derivativewrap)
      
    elif hasattr(obj,'fwd') or hasattr(obj,'adj'):
      @pyderivativegenerator
      def derivativewrap(f,nfwd,nadj):
        num_in = f.getNumInputs()
        num_out = f.getNumOutputs()
        
        @pyevaluate
        def der(f2):
          all_inputs = [f2.input(i) for i in range(f2.getNumInputs())]
          all_outputs = [f2.output(i) for i in range(f2.getNumOutputs())]
          inputs=all_inputs[:num_in]
          outputs=all_outputs[:num_out]
          fwd_seeds=zip(*[iter(all_inputs[num_in:num_in*(nfwd+1)])]*num_in)
          fwd_sens=zip(*[iter(all_outputs[num_out:num_out*(nfwd+1)])]*num_out)
          adj_seeds=zip(*[iter(all_inputs[num_in*(nfwd+1):])]*num_out)
          adj_sens=zip(*[iter(all_outputs[num_out*(nfwd+1):])]*num_in)
          if hasattr(obj,'fwd') and nfwd>0:
            obj.fwd(inputs,outputs,fwd_seeds,fwd_sens)
          if hasattr(obj,'adj') and nadj>0:
            obj.adj(inputs,outputs,adj_seeds,adj_sens)
          
        DerFun = CustomFunction(der,inputs+nfwd*inputs+nadj*outputs,outputs+nfwd*outputs+nadj*inputs)
        DerFun.setOption("name","CustomFunction_derivative")
        DerFun.init()
        return DerFun
 
      Fun.setOption("derivative_generator",derivativewrap)
    
    if not(hasattr(obj,'getDerivative')) and hasattr(obj,'fwd') and not hasattr(obj,'adj'):
      Fun.setOption("ad_mode","forward")
    if not(hasattr(obj,'getDerivative')) and not hasattr(obj,'fwd') and hasattr(obj,'adj'):
      Fun.setOption("ad_mode","reverse")
    return Fun
  
%}
#endif

%include "io_interface.i"
%include "io_scheme.i"
%include "io_scheme_vector.i"

%include <casadi/core/function/function.hpp>
%template(Pair_Function_Function) std::pair<casadi::Function,casadi::Function>;
%include <casadi/core/function/sx_function.hpp>
%include <casadi/core/function/mx_function.hpp>
%include <casadi/core/function/linear_solver.hpp>
%include <casadi/core/function/implicit_function.hpp>
%include <casadi/core/function/integrator.hpp>
%template(IntegratorVector) std::vector<casadi::Integrator>;
%include <casadi/core/function/simulator.hpp>
%include <casadi/core/function/control_simulator.hpp>
%include <casadi/core/function/nlp_solver.hpp>
%include <casadi/core/function/homotopy_nlp_solver.hpp>
%include <casadi/core/function/qp_solver.hpp>
%include <casadi/core/function/stabilized_qp_solver.hpp>
%include <casadi/core/function/lp_solver.hpp>
%include <casadi/core/function/sdp_solver.hpp>
%include <casadi/core/function/socp_solver.hpp>
%include <casadi/core/function/qcqp_solver.hpp>
%include <casadi/core/function/sdqp_solver.hpp>
%include <casadi/core/function/external_function.hpp>
%include <casadi/core/function/parallelizer.hpp>
%include <casadi/core/function/custom_function.hpp>
%include <casadi/core/functor.hpp>
%include <casadi/core/function/nullspace.hpp>
%include <casadi/core/function/dple_solver.hpp>
%include <casadi/core/function/dle_solver.hpp>
%include <casadi/core/function/lr_dple_solver.hpp>
%include <casadi/core/function/lr_dle_solver.hpp>
%include <casadi/core/function/cle_solver.hpp>

%include "autogenerated.i"

%include <casadi/core/casadi_options.hpp>
%include <casadi/core/casadi_meta.hpp>
%include <casadi/core/misc/integration_tools.hpp>
%include <casadi/core/misc/symbolic_nlp.hpp>
%include <casadi/core/misc/variable.hpp>
%include <casadi/core/misc/symbolic_ocp.hpp>
%include <casadi/core/misc/xml_file.hpp>
