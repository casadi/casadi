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

// Order in typemap matching: Lower value means will be checked first
%define PRECEDENCE_IVector 92 %enddef
%define PRECEDENCE_IVectorVector 92 %enddef
%define PRECEDENCE_PAIR_SLICE_SLICE 93 %enddef
%define PRECEDENCE_SLICE 94 %enddef
%define PRECEDENCE_PAIR_IVector_IVector 96 %enddef
%define PRECEDENCE_IMatrix 97 %enddef
%define PRECEDENCE_IMatrixVector 98 %enddef
%define PRECEDENCE_IMatrixVectorVector 98 %enddef
%define PRECEDENCE_DVector 99 %enddef
%define PRECEDENCE_DMatrix 100 %enddef
%define PRECEDENCE_DMatrixVector 101 %enddef
%define PRECEDENCE_DMatrixVectorVector 101 %enddef
%define PRECEDENCE_SX 103 %enddef
%define PRECEDENCE_SXVector 103 %enddef
%define PRECEDENCE_SXVectorVector 103 %enddef
%define PRECEDENCE_MX 104 %enddef
%define PRECEDENCE_MXVector 105 %enddef
%define PRECEDENCE_MXVectorVector 106 %enddef
%define PRECEDENCE_CREATOR 150 %enddef
%define PRECEDENCE_DERIVATIVEGENERATOR 21 %enddef
%define PRECEDENCE_CUSTOMEVALUATE 21 %enddef
%define PRECEDENCE_CALLBACK 21 %enddef
%define PRECEDENCE_GENERICTYPE 22 %enddef
%define PRECEDENCE_DICT 21 %enddef

  /// Data structure in the target language holding data
#ifdef SWIGPYTHON
%{
#define GUESTOBJECT PyObject
%}
#elif defined(SWIGMATLAB)
%{
#define GUESTOBJECT mxArray
%}
#else
%{
#define GUESTOBJECT void
%}
#endif

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
namespace std {
  template<class T> class vector {};
  template<class A, class B> class pair {};
  template<class A, class B> class map {};
}
#elif !defined(SWIGPYTHON)
%include "std_string.i"
%include "std_vector.i"
%include "std_pair.i"
%include "std_map.i"
#else // SWIGPYTHON
 /* TODO(@jaeandersson): Not maintainable, refactor */

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


%include "std_string.i"

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

%include "std_pair.i"
%include "std_map.i"
#endif // SWIGPYTHON


%template() std::vector<std::string>;
%template() std::vector<bool> ;
%template() std::vector<std::vector<bool> > ;
%template() std::vector< std::vector<std::vector<bool> > > ;
%template() std::vector<unsigned char>;
%template() std::vector<int>;
%template() std::vector<std::vector<int> > ;
%template() std::vector< std::vector<std::vector<int> > > ;
%template() std::vector<double>;
%template() std::vector<std::vector<double> > ;
%template() std::vector< std::vector<std::vector<double> > > ;

#ifndef SWIGMATLAB
%template() std::pair<int,int>;
%template() std::vector< std::pair<int,int> >;
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

#ifndef SWIGXML
%include "typemaps.i"
#endif

/// Generic typemap structure


/// Check if Python object is of type T
%fragment("is_a", "header") {
  bool is_null(GUESTOBJECT *p) {
#ifdef SWIGPYTHON
    if (p == Py_None) return true;
#endif
#ifdef SWIGMATLAB
    if (p == 0) return true;
#endif
    return false;
  }

  bool is_a(GUESTOBJECT *p, swig_type_info *type) {
    void *dummy = 0;
    return !is_null(p) && SWIG_ConvertPtr(p, &dummy, type, 0) >= 0;
  }
}

%fragment("vector_size", "header", fragment="is_a") {
  int vector_size(GUESTOBJECT * p) {
#ifdef SWIGPYTHON
    if (PySequence_Check(p)
        && !PyString_Check(p)
        && !is_a(p, $descriptor(casadi::SX*))
        && !is_a(p, $descriptor(casadi::MX*))
        && !is_a(p, $descriptor(casadi::Matrix<int>*))
        && !is_a(p, $descriptor(casadi::Matrix<double>*))
        && !PyObject_HasAttrString(p,"__DMatrix__")
        && !PyObject_HasAttrString(p,"__SX__")
        && !PyObject_HasAttrString(p,"__MX__")) {
      return PySequence_Size(p);
    } else {
      return -1;
    }
#endif // SWIGPYTHON
#ifdef SWIGMATLAB
    if (mxGetClassID(p)==mxCELL_CLASS && mxGetM(p)==1) {
      return mxGetN(p);
    } else {
      return -1;
    }
#endif // SWIGMATLAB
    return -1;
  }
 }

%fragment("to_vector", "header") {
  int to_vector(GUESTOBJECT * p, void *mv, int (*f)(GUESTOBJECT *p, void *mv, int offs)) {
#ifdef SWIGPYTHON
    PyObject *it = PyObject_GetIter(p);
    if (!it) {
      PyErr_Clear();
      return false;
    }
    PyObject *pe;
    int size = PySequence_Size(p);
    if (size==-1) {
      PyErr_Clear();
      return false;
    }
    int i=0;
    while ((pe = PyIter_Next(it))) {
      // Iterate over the sequence inside the sequence
      if (!f(pe, mv, i++)) {
        Py_DECREF(pe);
        Py_DECREF(it);
        return false;
      }
      Py_DECREF(pe);
    }
    Py_DECREF(it);
    return true;
#endif // SWIGPYTHON
#ifdef SWIGMATLAB
    if (mxGetClassID(p)==mxCELL_CLASS && mxGetM(p)==1) {
      int sz = mxGetN(p);
      for (int i=0; i<sz; ++i) {
        GUESTOBJECT *pi = mxGetCell(p, i);
        if (pi==0 || !f(pi, mv, i)) return false;
      }
      return true;
    }
#endif // SWIGMATLAB
    return false;
  }
 }

%fragment("make_vector", "header", fragment="vector_size,to_vector") {
  template<typename T>
  bool make_vector(GUESTOBJECT * p, std::vector<T>* m, int (*f)(GUESTOBJECT *p, void *mv, int offs)) {
    int sz = vector_size(p);
    if (sz<0) return false;
    if (m) m->resize(sz);
    if (sz>0 && !to_vector(p, m ? &m->front() : 0, f)) return false;
    return true;
  }
 }

%fragment("make_vector2", "header", fragment="make_vector") {
  template<typename T>
    bool make_vector2(GUESTOBJECT * p, std::vector<std::vector<T> >* m, int (*f)(GUESTOBJECT *p, void *mv, int offs)) {
    int sz = vector_size(p);
    if (sz<0) return false;
    if (m) m->resize(sz);
#ifdef SWIGPYTHON
    PyObject *it = PyObject_GetIter(p);
    if (!it) { PyErr_Clear();  return false;}
    PyObject *pe;
    int i=0;
    while ((pe = PyIter_Next(it))) {
      if (!make_vector(pe, m ? &m->at(i++) : 0, f)) {
        Py_DECREF(pe);
        Py_DECREF(it);
        return false;
      }
      Py_DECREF(pe);
    }
    Py_DECREF(it);
    return true;
#endif // SWIGPYTHON
#ifdef SWIGMATLAB
    for (int i=0; i<sz; ++i) {
      GUESTOBJECT *pi = mxGetCell(p, i);
      if (pi==0 || !make_vector(pi, m ? &m->at(i) : 0, f)) return false;
    }
    return true;
#endif // SWIGMATLAB
    return false;
  }
}

%define %casadi_in_typemap(xName, xType...)
%typemap(in, noblock=1, fragment="to"{xName}) xType (xType m) {
  if (!to_##xName($input, &m)) SWIG_exception_fail(SWIG_TypeError,"Failed to convert input to xName.");
  $1 = m;
}
%enddef

%fragment("conv_constref", "header", fragment="is_a") {
  template<typename T>
  bool conv_constref(GUESTOBJECT *p, T* &ptr, T &m, swig_type_info *type, int (*f)(GUESTOBJECT *p, void *mv, int offs)) {
    if (is_null(p)) {
      m = T();
      ptr = &m;
      return true;
    }
    if (!is_null(p) && SWIG_ConvertPtr(p, (void **)&ptr, type, 0) == -1) {
      if (!f(p, &m, 0)) return false;
      ptr = &m;
    }
    return true;
  }
 }

%define %casadi_in_typemap_constref2(xName, xType...)
%typemap(in, noblock=1, fragment="to"{xName}) const xType & (xType m) {
  if (!to_##xName($input, &m)) SWIG_exception_fail(SWIG_TypeError,"Failed to convert input to xName.");
  $1 = &m;
}
%enddef

%define %casadi_in_typemap_constref(xName, xType...)
%typemap(in, noblock=1, fragment="to"{xName}) const xType & (xType m) {
  if (!conv_constref($input, $1, m, $1_descriptor, to_##xName)) SWIG_exception_fail(SWIG_TypeError,"Failed to convert input to xName.");
 }
%enddef

%fragment("conv_vector", "header", fragment="make_vector") {
  template<typename T>
  bool conv_vector(GUESTOBJECT *p, std::vector<T>* &ptr, std::vector<T> &m, swig_type_info *type, int (*f)(GUESTOBJECT *p, void *mv, int offs)) {
    if (SWIG_ConvertPtr(p, (void **) &ptr, type, 0) == -1) {
      if (!make_vector(p, &m, f)) return false;
      ptr = &m;
    }
    return true;
  }
 }

%define %casadi_in_typemap_vector(xName,xType...)
%typemap(in, noblock=1, fragment="to"{xName}) const std::vector< xType > & (std::vector< xType > m) {
  if (!conv_vector($input, $1, m, $1_descriptor, to_##xName)) SWIG_exception_fail(SWIG_TypeError,"Cannot convert input to std::vector<xName>.");
 }
%enddef

%fragment("conv_vector2", "header", fragment="make_vector2") {
template<typename T>
bool conv_vector2(GUESTOBJECT *p, std::vector< std::vector<T> >* &ptr, std::vector< std::vector<T> > &m,
                 swig_type_info *type, int (*f)(GUESTOBJECT *p, void *mv, int offs)) {
  if (SWIG_ConvertPtr(p, (void **) &ptr, type, 0) == -1) {
    if (!make_vector2(p, &m, f)) return false;
    ptr = &m;
  }
  return true;
 }
}

%define %casadi_in_typemap_vector2(xName,xType...)
%typemap(in, noblock=1, fragment="conv_vector2") const std::vector< std::vector< xType > > & (std::vector< std::vector< xType > > m) {
  if (!conv_vector2($input, $1, m, $1_descriptor, to_##xName)) SWIG_exception_fail(SWIG_TypeError,"Cannot convert input to std::vector< std::vector<xName> >.");
 }
%enddef

%define %casadi_freearg_typemap(xType...)
%typemap(freearg, noblock=1) xType {}
%enddef

%define %casadi_typecheck_typemap(xName, xPrec, xType...)
%typemap(typecheck, noblock=1, fragment="to"{xName}, precedence=xPrec) xType {
  $1 = to_##xName($input, static_cast< xType *>(0));
 }
%enddef

%define %casadi_typecheck_typemap_constref(xName, xPrec, xType...)
%typemap(typecheck, noblock=1, fragment="to"{xName}, precedence=xPrec) const xType& {
  $1 = to_##xName($input, static_cast< xType *>(0));
 }
%enddef

%define %casadi_typecheck_typemap_vector(xName, xPrec, xType...)
%typemap(typecheck, noblock=1, precedence=xPrec, fragment="is_a,to_vector") const std::vector< xType > & {
  $1 = is_a($input,$1_descriptor) || make_vector< xType >($input, 0, to_##xName);
 }
%enddef

%define %casadi_typecheck_typemap_vector2(xName, xPrec, xType...)
%typemap(typecheck, noblock=1, precedence=xPrec, fragment="is_a,make_vector2") const std::vector< std::vector< xType > > & {
  $1 = is_a($input, $1_descriptor) || make_vector2< xType >($input, 0, to_##xName);
 }
%enddef

%fragment("conv_genericmatrix", "header", fragment="is_a") {
template<typename T>
bool conv_genericmatrix(GUESTOBJECT *p, casadi::GenericMatrix< T >* &ptr, T &m,
                        swig_type_info *type, int (*f)(GUESTOBJECT *p, void *mv, int offs)) {
  if (is_null(p) || SWIG_ConvertPtr(p, (void **) &ptr, type, 0) == -1) {
    if (!f(p, &m, 0)) return false;
    ptr = &m;
  }
  return true;
 }
}

%define %casadi_in_typemap_genericmatrix(xName, xType...) 
%typemap(in, noblock=1, fragment="conv_genericmatrix") const casadi::GenericMatrix< xType > & (xType m) {
  if (!conv_genericmatrix($input, $1, m, $descriptor(xType *), to_##xName)) SWIG_exception_fail(SWIG_TypeError, "Input type conversion failure ($1_type)");
 }
%enddef

%define %casadi_typecheck_typemap_genericmatrix(xName, xPrec, xType...)
%typemap(typecheck, noblock=1, precedence=xPrec, fragment="is_a") const casadi::GenericMatrix< xType > & {
  $1 = is_a($input, $descriptor(xType *)) || to_##xName($input, 0);
 }
%enddef

%define %casadi_typemaps(xName, xPrec, xType...)
%casadi_in_typemap(xName, xType)
%casadi_freearg_typemap(xType)
%casadi_typecheck_typemap(xName, xPrec, xType)
%enddef

%define %casadi_typemaps_constref(xName, xPrec, xType...)
%casadi_in_typemap_constref(xName, xType)
%casadi_freearg_typemap(const xType&)
%casadi_typecheck_typemap_constref(xName, xPrec, xType)
%enddef

%define %casadi_typemaps_vector(xName, xPrec, xType...)
%casadi_in_typemap_vector(xName, xType)
%casadi_freearg_typemap(const std::vector< xType >&)
%casadi_typecheck_typemap_vector(xName, xPrec, xType)
%enddef

%define %casadi_typemaps_vector2(xName, xPrec, xType...)
%casadi_in_typemap_vector2(xName, xType)
%casadi_freearg_typemap(const std::vector< std::vector< xType > >&)
%casadi_typecheck_typemap_vector2(xName, xPrec, xType)
%enddef

%define %casadi_typemaps_genericmatrix(xName, xPrec, xType...)
%casadi_in_typemap_genericmatrix(xName, xType)
%casadi_freearg_typemap(const casadi::GenericMatrix< xType >  &)
%casadi_typecheck_typemap_genericmatrix(xName, xPrec, xType)
%enddef

#ifdef SWIGPYTHON

%define %my_creator_typemap(Precedence,Type...)

%typemap(in) Type {
  int res = SWIG_ConvertFunctionPtr($input, (void**)(&$1), $descriptor);
  if (!SWIG_IsOK(res)) {
    if (PyType_Check($input) && PyObject_HasAttrString($input,"creator")) {
      PyObject *c = PyObject_GetAttrString($input,"creator");
      res = SWIG_ConvertFunctionPtr(c, (void**)(&$1), $descriptor);
      Py_DECREF(c);
    }
    if (!SWIG_IsOK(res)) {
      %argument_fail(res,"$type",$symname, $argnum); 
    }
  }
}

%typemap(typecheck,precedence=Precedence) Type { 
  void *ptr = 0;
  int res = SWIG_ConvertFunctionPtr($input, &ptr, $descriptor);
  $1 = SWIG_CheckState(res);
  if (!$1 && PyType_Check($input) && PyObject_HasAttrString($input,"creator")) {
    PyObject *c = PyObject_GetAttrString($input,"creator");
    res = SWIG_ConvertFunctionPtr(c, &ptr, $descriptor);
    $1 = SWIG_CheckState(res);
    Py_DECREF(c);
  };
}
%enddef

#else // SWIGPYTHON
%define %my_creator_typemap(Precedence,Type...)
%enddef
#endif // SWIGPYTHON



// Create an output typemap for a const ref such that a copy is made
%define %outputConstRefCopy(Type)
%typemap(out, noblock=1) const Type & {
  $result = SWIG_NewPointerObj((new Type(*$1)), $descriptor(Type*), SWIG_POINTER_OWN |  0 );
}
%enddef

#ifdef SWIGPYTHON
%inline%{
// Indicates that self is derived from parent.
// by setting a special attribute __swigref_parent__
void PySetParent(PyObject* self, PyObject* parent) {
  Py_INCREF(parent);
  PyObject_SetAttrString(self,"__swigref_parent__",parent);
}

// Look for the __swigref_parent__ attribute and DECREF the parent if there is one
void PyDECREFParent(PyObject* self) {
  if (!PyObject_HasAttrString(self,"__swigref_parent__")) return;
  
  PyObject* parent = PyObject_GetAttrString(self,"__swigref_parent__");
  if (!parent) return;
  Py_DECREF(parent); // Once for PyObject_GetAttrString
  if (!parent) return;
  if (parent!=Py_None) {
    Py_DECREF(parent); // Once for the actual DECREF
  }
}

%}
#endif //SWIGPYTHON

// Create an output typemap for a ref such that ownership is implied
// The following text is obsolete:
// We make use of SWIG_POINTER_OWN to ensure that the "delete_*" routine is called, where we have put another hook via %unrefobject.
// We do not really imply that this SWIG objects owns the pointer
// We are actually abusing the term SWIG_POINTER_OWN: a non-const ref is usually created with SWIG_NewPointerObj(..., 0 |  0 )
%define %outputRefOwn(Type)
%typemap(out, noblock=1) Type & {
  $result = SWIG_NewPointerObj($1, $descriptor(Type*), 0 |  0 );
  PySetParent($result, obj0);
}
%typemap(out, noblock=1) const Type & {
   $result = swig::from(static_cast< Type * >($1));
   PySetParent($result, obj0);
}
%extend Type {
%pythoncode%{
    def __del__(self):
      if not(_casadi is None):
         _casadi.PyDECREFParent(self)

%}

}
%enddef

// Convert reference output to a new data structure
%define %outputRefNew(Type)
%typemap(out, noblock=1) Type & {
   $result = swig::from(static_cast< Type >(*$1));
}
%typemap(out, noblock=1) const Type & {
   $result = swig::from(static_cast< Type >(*$1));
}
%extend Type {
%pythoncode%{
    def __del__(self):
      if not(_casadi is None):
         _casadi.PyDECREFParent(self)

%}

}
%enddef

#ifdef SWIGPYTHON
%inline%{
/** Check PyObjects by class name */
bool PyObjectHasClassName(PyObject* p, const char * name) {
  PyObject * classo = PyObject_GetAttrString( p, "__class__");
  PyObject * classname = PyObject_GetAttrString( classo, "__name__");
  
  bool ret = strcmp(PyString_AsString(classname),name)==0;
  Py_DECREF(classo);Py_DECREF(classname);
	return ret;
}

%}
#endif // SWIGPYTHON

%fragment("try_copy", "header", fragment="is_a") {
  template<typename FromType>
  struct CopyTraits {
    template<typename ToType>
    static bool try_copy(GUESTOBJECT *p, swig_type_info *type, ToType* m) {
      FromType *mp = 0;
      if (!is_null(p) && SWIG_ConvertPtr(p, (void **) &mp, type, 0) != -1) {
        if (m) *m=*mp;
        return true;
      }
      return false;
    }
  };
}
#define TRY_COPY(fromPtr, fromClass, fromType, toPtr) CopyTraits<fromClass>::try_copy(fromPtr, fromType, toPtr)

#ifdef SWIGMATLAB
// Get sparsity pattern
%fragment("get_sparsity", "header") {
  Sparsity getSparsity(const mxArray* p) {
    // Get sparsity pattern
    size_t nrow = mxGetM(p);
    size_t ncol = mxGetN(p);

    if (mxIsSparse(p)) {
      // Sparse storage in MATLAB
      mwIndex *Jc = mxGetJc(p);
      mwIndex *Ir = mxGetIr(p);

      // Store in vectors
      std::vector<int> colind(ncol+1);
      std::copy(Jc, Jc+colind.size(), colind.begin());
      std::vector<int> row(colind.back());
      std::copy(Ir, Ir+row.size(), row.begin());

      // Create pattern and return
      return casadi::Sparsity(nrow, ncol, colind, row);
    } else {
      return casadi::Sparsity::dense(nrow, ncol);
    }
  }
 }

// Number of nonzeros
%fragment("get_nnz", "header") {
  size_t getNNZ(const mxArray* p) {
    // Dimensions
    size_t nrow = mxGetM(p);
    size_t ncol = mxGetN(p);
    if (mxIsSparse(p)) {
      // Sparse storage in MATLAB
      mwIndex *Jc = mxGetJc(p);
      return Jc[ncol];
    } else {
      return nrow*ncol;
    }
  }
}
#endif // SWIGMATLAB

%{
#define SWIG_Error_return(code, msg)  { std::cerr << "Error occured in CasADi SWIG interface code:" << std::endl << "  "<< msg << std::endl;SWIG_Error(code, msg); return 0; }
%}

// Forward declarations
%fragment("fwd", "header",
          fragment="vector_size,to_vector,make_vector,make_vector2,conv_constref,conv_vector,conv_vector2,conv_genericmatrix,try_copy"
#ifdef SWIGMATLAB
          ,fragment="get_sparsity,get_nnz"
#endif // SWIGMATLAB
          ) {
  namespace casadi {
    // Traits for CasADi typemaps
    template <typename Type> struct traits_casptr {};

    /* Convert a pointer in interfaced language to C++
     * Input: GUESTOBJECT pointer p
     * Output: Pointer to pointer: At input, pointer to pointer to temporary
     * The routine will either:
     *   - Do nothing, if 0
     *   - Change the pointer
     *   - Change the temporary object
     * Returns true upon success, else false
     */
    template <typename Type> bool casptr(GUESTOBJECT *p, Type** m) {
      return traits_casptr<Type>::casptr(p, m);
    }
  } // namespace casadi

  int to_int(GUESTOBJECT *p, void *mv, int offs=0);
  int to_double(GUESTOBJECT *p, void *mv, int offs=0);
  int to_Dict(GUESTOBJECT *p, void *mv, int offs=0);
  int to_GenericType(GUESTOBJECT *p, void *mv, int offs=0);
  int to_DerivativeGenerator(GUESTOBJECT *p, void *mv, int offs=0);
  int to_CustomEvaluate(GUESTOBJECT *p, void *mv, int offs=0);
  int to_Callback(GUESTOBJECT *p, void *mv, int offs=0);
  int to_DVector(GUESTOBJECT *p, void *mv, int offs=0);
  int to_IVector(GUESTOBJECT *p, void *mv, int offs=0);
  int to_SX(GUESTOBJECT *p, void *mv, int offs=0);
  int to_MX(GUESTOBJECT *p, void *mv, int offs=0);
  int to_DMatrix(GUESTOBJECT *p, void *mv, int offs=0);
  int to_IMatrix(GUESTOBJECT *p, void *mv, int offs=0);
  int to_Slice(GUESTOBJECT *p, void *mv, int offs=0);
  int to_string(GUESTOBJECT *p, void *mv, int offs=0);
  int to_Function(GUESTOBJECT *p, void *mv, int offs=0);

  inline int to_Casadi(GUESTOBJECT *p, casadi::MX* m) { return to_MX(p, m);}
  inline int to_Casadi(GUESTOBJECT *p, casadi::Matrix<double>* m) { return to_DMatrix(p, m);}
  inline int to_Casadi(GUESTOBJECT *p, casadi::Matrix<casadi::SXElement>* m) { return to_SX(p, m);}

  template<typename M> int to_Map(GUESTOBJECT *p, std::map<std::string, M> *m);
  template<typename M> int to_PairMap(GUESTOBJECT *p, std::pair<std::map<std::string, M>, std::vector<std::string> > *m);

  GUESTOBJECT * from_GenericType(const casadi::GenericType &a);
  GUESTOBJECT * from_Dict(const casadi::GenericType::Dict &a);

  swig_type_info * type_SX() { return $descriptor(casadi::Matrix<casadi::SXElement> *); }
  swig_type_info * type_DMatrix() { return $descriptor(casadi::Matrix<double> *); }
  swig_type_info * type_IMatrix() { return $descriptor(casadi::Matrix<int> *); }
  swig_type_info * type_MX() { return $descriptor(casadi::MX *); }
}

#ifdef SWIGPYTHON
%{
  // Returns a new reference
  PyObject * getReturnType(PyObject* p) {
    if (!p) return 0;
    if (!PyCallable_Check(p)) return 0;
    if (!PyObject_HasAttrString( p, "func_annotations")) return 0;
    PyObject * func_annotations = PyObject_GetAttrString(p,"func_annotations");
    if (!PyDict_Check(func_annotations)) {
      Py_DECREF(func_annotations);
      return 0;
    }
    PyObject * return_type = PyDict_GetItemString(func_annotations, "return"); // Borrowed
    Py_INCREF(return_type); // Make a new reference
    Py_DECREF(func_annotations);
    if (return_type==0) {
      return 0;
    }
    if (!PyType_Check(return_type) && return_type!=Py_None) {
      Py_DECREF(return_type);
      return 0;
    }
    return return_type;
  }

  // Returns a new reference
  PyObject * getCasadiObject(const std::string &s) {
    PyObject* pPyObjectModuleName = PyString_FromString("casadi");
    if (pPyObjectModuleName) {
      PyObject* pObjectModule = PyImport_Import(pPyObjectModuleName);
      Py_DECREF(pPyObjectModuleName);
      if (pObjectModule) {
        PyObject* pObjectDict = PyModule_GetDict(pObjectModule); // Borrowed
        Py_DECREF(pObjectModule);
        if (pObjectDict) {
          PyObject* ret = PyDict_GetItemString(pObjectDict,  s.c_str()); // Borrowed
          if (ret) {
            // New reference
            Py_INCREF(ret);
            return ret;
          }
        }
      }
    }

    // Unsuccessful
    PyErr_Clear();
    return 0;
  }

#include <casadi/core/functor_internal.hpp>

  namespace casadi {
    class FunctorPythonInternal {
    public:
    FunctorPythonInternal(PyObject *p) : p_(p) { Py_INCREF(p_); }
      ~FunctorPythonInternal() { Py_DECREF(p_); }
    protected:
      PyObject *p_;
    };
  
    class DerivativeGeneratorPythonInternal : public DerivativeGeneratorInternal, FunctorPythonInternal {
      friend class DerivativeGeneratorPython;
    
    DerivativeGeneratorPythonInternal(PyObject *p) : FunctorPythonInternal(p) {}
      virtual Function call(Function& fcn, int ndir, void* user_data);
      virtual DerivativeGeneratorPythonInternal* clone() const { return new DerivativeGeneratorPythonInternal(p_); }
    };
   
    class CustomEvaluatePythonInternal : public CustomEvaluateInternal, FunctorPythonInternal {
      friend class CustomEvaluatePython;
    
    CustomEvaluatePythonInternal(PyObject *p) : FunctorPythonInternal(p) {}
      virtual void call(CustomFunction& fcn, void* user_data);
      virtual CustomEvaluatePythonInternal* clone() const { return new CustomEvaluatePythonInternal(p_); }
    };
  
    class DerivativeGeneratorPython : public DerivativeGenerator {
    public:
      DerivativeGeneratorPython(PyObject *p) { assignNode(new DerivativeGeneratorPythonInternal(p)); }
    };
  
    class CallbackPythonInternal : public CallbackInternal, FunctorPythonInternal {
      friend class CallbackPython;
    
    CallbackPythonInternal(PyObject *p) : FunctorPythonInternal(p) {}
      virtual int call(Function& fcn, void* user_data);
      virtual CallbackPythonInternal* clone() const { return new CallbackPythonInternal(p_); }
    };
  
    class CustomEvaluatePython : public CustomEvaluate {
    public:
      CustomEvaluatePython(PyObject *p) { assignNode(new CustomEvaluatePythonInternal(p)); }
    };

    class CallbackPython : public Callback {
    public:
      CallbackPython(PyObject *p) { assignNode(new CallbackPythonInternal(p)); }
    };
  }
  %}
#endif

#ifdef SWIGMATLAB
%{
#include <streambuf>
#include <ostream>
  namespace casadi {
    /** Stream buffer to allow printing to mex conveniently
        Adapted from answer to stackoverflow question 243696.
    */
    class mex_buf : public std::streambuf {
    public:
      mex_buf() {}
    protected:
      virtual int_type overflow(int_type ch) {
        if(ch != traits_type::eof()) {
          mexPrintf("%c", static_cast<char>(ch));
        }
        return ch;
      }
      /* virtual int sync() { // needed?
         mexEvalString("drawnow;");
         return 0;
      } */
      virtual std::streamsize xsputn(const char* s, std::streamsize num) {
        mexPrintf("%.*s", static_cast<int>(num), s);
        return num;
      }
    };

    // Corresponding output stream
    class mexostream : public std::ostream {
    protected:
      mex_buf buf;
    public:
      mexostream() : std::ostream(&buf) {}
    };

    // Instantiation (cf. cout)
    static mexostream mexout;
  } // namespace casadi
%}
#endif

%template() std::vector< casadi::Sparsity > ;
%template() std::vector< std::vector< casadi::Sparsity> > ;
%template() std::vector<casadi::Matrix<casadi::SXElement> > ;
%template() std::vector< std::vector<casadi::Matrix<casadi::SXElement> > > ;
%template() std::vector<casadi::MX>;
%template() std::vector< std::vector<casadi::MX> >;
%template() std::vector<casadi::Matrix<int> > ;
%template() std::vector<casadi::Matrix<double> > ;
%template() std::vector< std::vector<casadi::Matrix<double> > > ;
%template() std::vector< std::vector<casadi::Matrix<int> > > ;

%fragment("to"{int}, "header", fragment="fwd") {
  // Traits specialization for int
  namespace casadi {
    template <> struct traits_casptr<int> {
      static bool casptr(GUESTOBJECT *p, int** m) {
        // Standard typemaps
        if (SWIG_IsOK(SWIG_AsVal_int(p, m ? *m : 0))) return true;
        // No match
        return false;
      }
    };
  } // namespace casadi

  int to_int(GUESTOBJECT *p, void *mv, int offs) {
    int *m = static_cast<int*>(mv);
    if (m) m += offs;

    // Refactored code
    if (casadi::casptr(p, m ? &m : 0)) return true;

    // Scalar IMatrix
    if (is_a(p, type_IMatrix())) {
      casadi::Matrix<int> *temp;
      SWIG_ConvertPtr(p, (void **) &temp, type_IMatrix(), 0 );
      if (temp->isScalar(true)) {
        if (m) *m = temp->getIntValue();
        return true;
      }
      return false;
    }
#ifdef SWIGPYTHON
    // Object has __int__ attribute
    if (PyObject_HasAttrString(p,"dtype")) {
      PyObject *r = PyObject_GetAttrString(p,"dtype");
      if (!PyObject_HasAttrString(r,"kind")) {
        Py_DECREF(r);
        return false;
      }
      PyObject *k = PyObject_GetAttrString(r,"kind");
      if (!PyObject_HasAttrString(p,"__int__")) {
        Py_DECREF(k);
        Py_DECREF(r);
        return false;
      }
      char name[] = "__int__";
      PyObject *mm = PyObject_CallMethod(p, name,0);
      if (!mm) {
        PyErr_Clear();
        Py_DECREF(k);
        Py_DECREF(r);
        return false;
      }
      char *kk = PyString_AsString(k);
      bool result =  kk[0]=='i';
      Py_DECREF(k);
      Py_DECREF(r);
      if (result) {
        if (m) *m = PyLong_AsLong(mm);
      }
      Py_DECREF(mm);
      return result;
    }
#endif //SWIGPYTHON
    // Failure if reached this point
    return false;
  }
 }

%casadi_typemaps(int, SWIG_TYPECHECK_INTEGER, int)

%fragment("to"{double}, "header", fragment="fwd") {
  // Traits specialization for double
  namespace casadi {
    template <> struct traits_casptr<double> {
      static bool casptr(GUESTOBJECT *p, double** m) {
        // Standard typemaps
        if (SWIG_IsOK(SWIG_AsVal_double(p, m ? *m : 0))) return true;
        // No match
        return false;
      }
    };
  } // namespace casadi

  int to_double(GUESTOBJECT *p, void *mv, int offs) {
    double *m = static_cast<double*>(mv);
    if (m) m += offs;

    // Refactored code
    if (casadi::casptr(p, m ? &m : 0)) return true;

    // Scalar DMatrix
    if (is_a(p, type_DMatrix())) {
      casadi::Matrix<double> *ptr;
      SWIG_ConvertPtr(p, (void **) &ptr, type_DMatrix(), 0 );
      if (ptr->isScalar(true)) {
        if (m) *m = ptr->getValue();
        return true;
      }
      return false;
    }

    // Scalar IMatrix
    if (is_a(p, type_IMatrix())) {
      casadi::Matrix<int> *ptr;
      SWIG_ConvertPtr(p, (void **) &ptr, type_IMatrix(), 0 );
      if (ptr->isScalar(true)) {
        if (m) *m = ptr->getIntValue();
        return true;
      }
      return false;
    }
#ifdef SWIGPYTHON
    // Has dtype attribute
    if (PyObject_HasAttrString(p,"dtype")) {
      PyObject *r = PyObject_GetAttrString(p,"dtype");
      if (!PyObject_HasAttrString(r,"kind")) {
        Py_DECREF(r);
        return false;
      }
      PyObject *k = PyObject_GetAttrString(r,"kind");
      if (!PyObject_HasAttrString(p,"__float__")) {
        Py_DECREF(k);
        Py_DECREF(r);
        return false;
      }
      char name[] = "__float__";
      PyObject *mm = PyObject_CallMethod(p, name,NULL);
      if (!mm) {
        PyErr_Clear();
        Py_DECREF(k);
        Py_DECREF(r);
        return false;
      }
      char *kk = PyString_AsString(k);
      bool result = kk[0]=='f' || kk[0]=='i';
      Py_DECREF(k);
      Py_DECREF(r); 
      if (result && m) *m = PyFloat_AsDouble(mm);
      Py_DECREF(mm);
      return result;
    }
#endif //SWIGPYTHON
    // Failure if reached this point
    return false;
  }
 }

%casadi_typemaps(double, SWIG_TYPECHECK_DOUBLE, double)
%casadi_typemaps_constref(double, SWIG_TYPECHECK_DOUBLE, double)

#ifdef SWIGPYTHON
%fragment("from"{GenericType}, "header", fragment="fwd") {
  GUESTOBJECT * from_GenericType(const casadi::GenericType &a) {
    GUESTOBJECT *p = 0;
    if (a.isBool()) {
      p=PyBool_FromLong(a.toBool());
    } else if (a.isInt()) {
      p=PyInt_FromLong(a.toInt());
    } else if (a.isDouble()) {
      p=PyFloat_FromDouble(a.toDouble());
    } else if (a.isString()) {
      p=PyString_FromString(a.toString().c_str());
    } else if (a.isIntVector()) {
      p = swig::from(a.toIntVector());
    } else if (a.isIntVectorVector()) {
      p = swig::from(a.toIntVectorVector());
    } else if (a.isDoubleVector()) {
      p = swig::from( a.toDoubleVector());
    }  else if (a.isStringVector()) {
      p = swig::from(a.toStringVector());
    } else if (a.isDict()) {
      p = from_Dict(a.toDict());
    } else if (a.isFunction()) {
      p = swig::from( a.toFunction());
    } else if (a.isNull()) {
      p = Py_None;
    }
    return p;
  }
}

%typemap(out, noblock=1, fragment="from"{GenericType}) casadi::GenericType {
  if(!($result = from_GenericType($1))) SWIG_exception_fail(SWIG_TypeError,"GenericType not yet implemented");
}

%typemap(out, fragment="from"{GenericType}) std::vector< casadi::GenericType > {
  PyObject* ret = PyList_New(0);
  std::vector< casadi::GenericType > & in = $1;
  for (int k=0 ; k < in.size(); ++k) {
    PyObject* rete;
    if (!(rete = from_GenericType(in[k]))) SWIG_exception_fail(SWIG_TypeError,"GenericType not yet implemented");
    PyList_Append(ret, rete);
  }
  $result = ret;
}

%fragment("from"{Dict}, "header", fragment="fwd") {
  GUESTOBJECT * from_Dict(const casadi::GenericType::Dict &a) {
    PyObject *p = PyDict_New();
    casadi::GenericType::Dict::const_iterator end = a.end();
    for (casadi::GenericType::Dict::const_iterator it = a.begin(); it != end; ++it) {
      PyObject * e = from_GenericType(it->second);
      if (!e) {
        Py_DECREF(p);
        return 0;
      }
      PyDict_SetItemString(p,(it->first).c_str(),e);
      Py_DECREF(e);
    }
    return p;
  }
}
%typemap(out, noblock=1, fragment="from"{Dict}) const casadi::GenericType::Dict&  {
  if(!($result = from_Dict(*$1))) SWIG_exception_fail(SWIG_TypeError,"GenericType not yet implemented");
}
#endif // SWIGPYTHON

%fragment("to"{string}, "header", fragment="fwd") {
  // Traits specialization for string
  namespace casadi {
    template <> struct traits_casptr<std::string> {
      static bool casptr(GUESTOBJECT *p, std::string** m) {
        // Treat Null
        if (is_null(p)) return false;

        // String already?
        if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(m),
                                      $descriptor(std::string*), 0))) {
          return true;
        }

#ifdef SWIGPYTHON
        if (PyString_Check(p)) {
          if (m) (*m)->clear();
          if (m) (*m)->append(PyString_AsString(p));
          return true;
        }
#endif // SWIGPYTHON

        // No match
        return false;
      }
    };
  } // namespace casadi

  int to_string(GUESTOBJECT *p, void *mv, int offs) {
    std::string *m = static_cast<std::string*>(mv);
    if (m) m += offs;

    // Call refactored version
    std::string *m_orig = m;
    if (casptr(p, m ? &m : 0)) {
      if (m!=m_orig) *m_orig=*m;
      return true;
    }

    return false;
  }
}

%fragment("to"{Function}, "header", fragment="fwd") {
  int to_Function(GUESTOBJECT *p, void *mv, int offs) {
    casadi::Function *m = static_cast<casadi::Function*>(mv);
    if (m) m += offs;
    casadi::Function *t = 0;
    int res = swig::asptr(p, &t);
    if(SWIG_CheckState(res) && t) {
      if(m) *m=*t;
      if (SWIG_IsNewObj(res)) delete t;
      return true;
    } else {
      return false;
    }
  }
}

%fragment("to"{Map}, "header", fragment="fwd") {
  template<typename M> int to_Map(GUESTOBJECT *p, std::map<std::string, M> *m) {
#ifdef SWIGPYTHON
    if (PyDict_Check(p)) {
      PyObject *key, *value;
      Py_ssize_t pos = 0;
      while (PyDict_Next(p, &pos, &key, &value)) {
        if (!PyString_Check(key)) return false;
        if (m) {
          M& v=(*m)[std::string(PyString_AsString(key))];
          if (!to_Casadi(value, &v)) return false;
        } else {
          if (!to_Casadi(value, static_cast<M*>(0))) return false;
        }
      }
      return true;
    }
#endif // SWIGPYTHON
    return false;
  }
 }

%fragment("to"{PairMap}, "header", fragment="fwd") {
  template<typename M> int to_PairMap(GUESTOBJECT *p, std::pair<std::map<std::string, M>, std::vector<std::string> > *m) {
#ifdef SWIGPYTHON
    if (PyTuple_Check(p) && PyTuple_Size(p)==2) {
      // Convert map
      GUESTOBJECT *p_first = PyTuple_GetItem(p,0);
      std::map<std::string, M> *m_first=0;
      if (m) m_first=&m->first;
      if (!to_Map(p_first, m_first)) return false;
      // Convert vector of strings
      GUESTOBJECT *p_second = PyTuple_GetItem(p,1);
      if (!PyList_Check(p_second)) return false;
      Py_ssize_t n = PyList_Size(p_second);
      if (m) m->second.resize(n);
      for (Py_ssize_t i=0; i!=n; ++i) {
        GUESTOBJECT *s_i = PyList_GetItem(p_second, i);
        if (!PyString_Check(s_i)) return false;
        if (m) m->second.at(i) = std::string(PyString_AsString(s_i));
      }
      // Success
      return true;
    }
#endif // SWIGPYTHON
    return false;
  }
 }


%fragment("to"{IVector}, "header", fragment="fwd,make_vector", fragment="to"{int}) {
  int to_IVector(GUESTOBJECT *p, void *mv, int offs) {
    std::vector<int> *m = static_cast<std::vector<int>*>(mv);
    if (m) m += offs;
#ifdef SWIGPYTHON
    std::vector< int > *mp = 0;
    if (SWIG_ConvertPtr(p, (void **) &mp, $descriptor(std::vector<int> *), 0) != -1) {
      if (m) *m=*mp;
      return true;
    } else if (is_array(p)) {
      if (!(array_numdims(p)==1 && array_type(p)!=NPY_OBJECT)) {
        //SWIG_Error_return(SWIG_TypeError, "std::vector<int>: array must be 1D and of a numeric type");
        return false;
      }
      int size = array_size(p,0);
      if (!array_is_native(p)) {
        //SWIG_Error_return(SWIG_TypeError, "std::vector<int>: array byte order should be native.");
        return false;
      }
      
      // Make sure we have a contigous array with int datatype
      int array_is_new_object;
      PyArrayObject* array = obj_to_array_contiguous_allow_conversion(p,NPY_INT,&array_is_new_object);
      if (!array) { // Trying LONG
        array = obj_to_array_contiguous_allow_conversion(p,NPY_LONG,&array_is_new_object);
        if (!array) { 
          //PyErr_Print() ; SWIG_Error_return(SWIG_TypeError, "std::vector<int>: no luck converting numpy array to int. Better don't use unsigned datatypes.");
          return false;
        }
        long* temp=(long*) array_data(array);
        if (m) {
          m->resize(size);
          for (int k=0; k<size; k++) (*m)[k]=temp[k];
        }
        return true;
      }
      int *d=(int*) array_data(array);

      if (m) m->assign( d, d+size );

                  
      // Free memory
      if (array_is_new_object)
        Py_DECREF(array); 
      return true;
    }
    return make_vector(p, m, to_int);
#endif // SWIGPYTHON
    return false;
  }
}

// TODO: Remove completely?
#ifdef SWIGPYTHON
%casadi_typemaps_constref(IVector, PRECEDENCE_IVector, std::vector<int>)
#endif

%fragment("to"{DVector}, "header", fragment="fwd,make_vector", fragment="to"{double}) {
  int to_DVector(GUESTOBJECT *p, void *mv, int offs) {
    std::vector<double> *m = static_cast<std::vector<double>*>(mv);
    if (m) m += offs;
#ifdef SWIGPYTHON
    std::vector< double > *mp = 0;
    if (SWIG_ConvertPtr(p, (void **) &mp, $descriptor(std::vector<double> *), 0) != -1) {
      if (m) *m=*mp;
      return true;
    } else if (is_array(p)) {
      if (!(array_numdims(p)==1 && array_type(p)!=NPY_OBJECT)) {
        return false;
        //SWIG_Error_return(SWIG_TypeError, "std::vector<int>: array must be 1D and of a numeric type");
      }
      int size = array_size(p,0);
      if (!array_is_native(p)) {
        return false;
        //SWIG_Error_return(SWIG_TypeError, "std::vector<double>: array byte order should be native.");
      }
      // Make sure we have a contigous array with double datatype
      int array_is_new_object;
      PyArrayObject* array = obj_to_array_contiguous_allow_conversion(p,NPY_DOUBLE,&array_is_new_object);
      if (!array) { 
        //PyErr_Print() ; SWIG_Error_return(SWIG_TypeError, "asMatrixDouble: no luck converting numpy array to double");
        return false;
      }
      double* d=(double*) array_data(array);
    
      if (m) m->assign( d, d+size );
    
                  
      // Free memory
      if (array_is_new_object)
        Py_DECREF(array); 
      return true;
    }
    return make_vector(p, m, to_double);
#endif // SWIGPYTHON
    return false;
  }
 }
// TODO: Remove completely?
#ifdef SWIGPYTHON
%casadi_typemaps_constref(DVector, PRECEDENCE_DVector, std::vector<double>)
#endif

%my_creator_typemap(PRECEDENCE_CREATOR, casadi::implicitFunctionCreator);
%my_creator_typemap(PRECEDENCE_CREATOR, casadi::linearSolverCreator);

#ifdef SWIGPYTHON
%fragment("to"{DerivativeGenerator}, "header", fragment="fwd", fragment="to"{Function}) {
  namespace casadi {
    Function DerivativeGeneratorPythonInternal::call(Function& fcn, int ndir, void* user_data) {
      casadi_assert(p_!=0);
      PyObject * ndir_py = PyInt_FromLong(ndir);
      PyObject * fcn_py = SWIG_NewPointerObj((new Function(static_cast< const Function& >(fcn))),
                                             $descriptor(casadi::Function *), SWIG_POINTER_OWN |  0 );
      if(!fcn_py) {
        Py_DECREF(ndir_py);
        throw CasadiException("DerivativeGeneratorPythonInternal: failed to convert Function to python");
      }
      PyObject *r = PyObject_CallFunctionObjArgs(p_, fcn_py, ndir_py, NULL);
      Py_DECREF(ndir_py);
      Py_DECREF(fcn_py);
      if (r) {
        Function ret;  
        if(!to_Function(r, &ret)) {
          Py_DECREF(r);
          throw CasadiException("DerivativeGeneratorPythonInternal: return type was not Function.");
        }
        Py_DECREF(r);
        return ret;
      } else {
        PyErr_Print();
        throw CasadiException("DerivativeGeneratorPythonInternal: python method execution raised an Error.");
      }
    }
  } // namespace casadi

  int to_DerivativeGenerator(GUESTOBJECT *p, void *mv, int offs) {
    casadi::DerivativeGenerator *m = static_cast<casadi::DerivativeGenerator*>(mv);
    if (m) m += offs;
    casadi::DerivativeGenerator *mp = 0;
    if (SWIG_ConvertPtr(p, (void **) &mp, $descriptor(casadi::DerivativeGenerator *), 0) != -1) {
      if (m) *m=*mp;
      return true;
    }
    PyObject* return_type = getReturnType(p);
    if (!return_type) return false;
    PyObject* function = getCasadiObject("Function");
    if (!function) { Py_DECREF(return_type); return false; }
    bool res = PyClass_IsSubclass(return_type,function);
    Py_DECREF(return_type);
    Py_DECREF(function);
    if (res) {
      if (m) *m = casadi::DerivativeGeneratorPython(p);
      return true;
    }
    return false;
  }
 }
%casadi_typemaps_constref(DerivativeGenerator, PRECEDENCE_DERIVATIVEGENERATOR, casadi::DerivativeGenerator)

%fragment("to"{CustomEvaluate}, "header", fragment="fwd") {
  namespace casadi {
    void CustomEvaluatePythonInternal::call(CustomFunction& fcn, void* user_data) {
      casadi_assert(p_!=0);
      PyObject * fcn_py = SWIG_NewPointerObj((new CustomFunction(static_cast< const CustomFunction& >(fcn))),
                                             $descriptor(casadi::CustomFunction *), SWIG_POINTER_OWN |  0 );
      if(!fcn_py) throw CasadiException("CustomEvaluatePythonInternal: failed to convert CustomFunction to python");
      PyObject *r = PyObject_CallFunctionObjArgs(p_, fcn_py, NULL);
      Py_DECREF(fcn_py);
      if (!r) {
        PyErr_Print();
        throw CasadiException("CustomEvaluatePythonInternal: Python method execution raised an Error.");
      }
      Py_DECREF(r);
    }
  } // namespace casadi

  int to_CustomEvaluate(GUESTOBJECT *p, void *mv, int offs) {
    casadi::CustomEvaluate *m = static_cast<casadi::CustomEvaluate*>(mv);
    if (m) m += offs;
    casadi::CustomEvaluate *mp = 0;
    if (SWIG_ConvertPtr(p, (void **) &mp, $descriptor(casadi::CustomEvaluate *), 0) != -1) {
      if (m) *m=*mp;
      return true;
    }
    PyObject* return_type = getReturnType(p);
    bool res = (return_type==Py_None) || !return_type;
    if (return_type) Py_DECREF(return_type);
    if (res) {
      if (m) *m = casadi::CustomEvaluatePython(p);
    }
    return res;
  }
 }
%casadi_typemaps_constref(CustomEvaluate, PRECEDENCE_CUSTOMEVALUATE, casadi::CustomEvaluate)

%fragment("to"{Callback}, "header", fragment="fwd") {
  namespace casadi {
    int CallbackPythonInternal::call(Function& fcn, void* user_data) {
      casadi_assert(p_!=0);
      PyObject * fcn_py = SWIG_NewPointerObj((new Function(static_cast< const Function& >(fcn))),
                                             $descriptor(casadi::CustomFunction *), SWIG_POINTER_OWN |  0 );
      if(!fcn_py) throw CasadiException("CallbackPythonInternal: failed to convert CustomFunction to python");
      PyObject *r = PyObject_CallFunctionObjArgs(p_, fcn_py, NULL);
      Py_DECREF(fcn_py);
      if (!r) {
        PyErr_Print();
        throw CasadiException("CallbackPythonInternal: python method execution raised an Error.");
      }
      int ret = 0;
      if ( to_int(r, 0)) to_int(r, &ret);
      Py_DECREF(r);
      return ret;
    }  
  } // namespace casadi

  int to_Callback(GUESTOBJECT *p, void *mv, int offs) {
    casadi::Callback *m = static_cast<casadi::Callback*>(mv);
    if (m) m += offs;
    casadi::Callback *mp = 0;
    if (p != Py_None && SWIG_ConvertPtr(p, (void **) &mp, $descriptor(casadi::Callback *), 0) != -1) {
      if (m) *m=*mp;
      return true;
    }
    PyObject* return_type = getReturnType(p);
    if (!return_type) return false;
    bool res = PyType_IsSubtype((PyTypeObject *)return_type,&PyInt_Type) || PyType_IsSubtype((PyTypeObject *)return_type,&PyLong_Type);
    Py_DECREF(return_type);
    if (res) {
      if (m) *m = casadi::CallbackPython(p);
      return true;
    }
    return false;
  }
 }
%casadi_typemaps_constref(Callback, PRECEDENCE_CALLBACK, casadi::Callback)
#endif

%fragment("to"{GenericType}, "header", fragment="fwd", fragment="to"{string},
          fragment="to"{IVector}, fragment="to"{DVector}, fragment="to"{Function}) {
  int to_GenericType(GUESTOBJECT *p, void *mv, int offs) {
    casadi::GenericType *m = static_cast<casadi::GenericType*>(mv);
    if (m) m += offs;
#ifndef SWIGPYTHON
    return false;
#else // SWIGPYTHON
    casadi::GenericType *mp = 0;
    if (p != Py_None && SWIG_ConvertPtr(p, (void **) &mp, $descriptor(casadi::GenericType *), 0) != -1) {
      if (m) *m=*mp;
      return true;
    }
    if (p==Py_None) {
      if (m) *m=casadi::GenericType();
    } else if (PyBool_Check(p)) {
      if (m) *m=casadi::GenericType((bool) PyInt_AsLong(p));
    } else if (PyInt_Check(p)) {
      if (m) *m=casadi::GenericType((int) PyInt_AsLong(p));
    } else if (PyFloat_Check(p)) {
      if (m) *m=casadi::GenericType(PyFloat_AsDouble(p));
    } else if (to_string(p, 0)) {
      std::string temp;
      if (!to_string(p, &temp)) return false;
      if (m) *m = casadi::GenericType(temp);
    } else if (to_IVector(p, 0)) {
      std::vector<int> temp;
      if (!to_IVector(p, &temp)) return false;
      if (m) *m = casadi::GenericType(temp);
    } else if (to_DVector(p, 0)) {
      std::vector<double> temp;
      if (!to_DVector(p, &temp)) return false;
      if (m) *m = casadi::GenericType(temp);
    } else if (make_vector(p, static_cast<std::vector<std::string>*>(0), to_string)) {
      std::vector<std::string> temp;
      if (!make_vector(p, &temp, to_string)) return false;
      if (m) *m = casadi::GenericType(temp);
    } else if (PyType_Check(p) && PyObject_HasAttrString(p,"creator")) {
      PyObject *c = PyObject_GetAttrString(p,"creator");
      if (!c) return false;
      PyObject* gt = getCasadiObject("GenericType");
      if (!gt) return false;

      PyObject* args = PyTuple_New(1);
      PyTuple_SetItem(args,0,c);
    
      PyObject* g = PyObject_CallObject(gt,args);
    
      Py_DECREF(args);
      Py_DECREF(gt);
    
      if (g) {
        int result = to_GenericType(g, m);
        Py_DECREF(g);
        return result;
      }
      if (!g) { PyErr_Clear();  return false;}
    
    } else if (to_Function(p, 0)) {
      casadi::Function temp;
      if (!to_Function(p, &temp)) return false;
      if (m) *m = casadi::GenericType(temp);
    } else if (to_Dict(p, 0) || to_DerivativeGenerator(p, 0) || to_Callback(p, 0)) {
      PyObject* gt = getCasadiObject("GenericType");
      if (!gt) return false;

      PyObject* args = PyTuple_New(1);
      Py_INCREF(p); // Needed because PyTuple_SetItem steals the reference
      PyTuple_SetItem(args,0,p);
    
      PyObject* g = PyObject_CallObject(gt,args);
    
      Py_DECREF(args);
      Py_DECREF(gt);
    
      if (g) {
        int result = to_GenericType(g, m);
        Py_DECREF(g);
        return result;
      }
      if (!g) {
        PyErr_Clear();
        return false;
      }
    } else {
      return false;
    }
    return true;
#endif // SWIGPYTHON
  }
 }

#ifdef SWIGPYTHON
%fragment("to"{Dict}, "header", fragment="fwd") {
  int to_Dict(GUESTOBJECT *p, void *mv, int offs) {
    casadi::GenericType::Dict *m = static_cast<casadi::GenericType::Dict*>(mv);
    if (m) m += offs;
#ifndef SWIGPYTHON
    return false;
#else // SWIGPYTHON
    casadi::GenericType::Dict *mp = 0;
    if (p != Py_None && SWIG_ConvertPtr(p, (void **) &mp, $descriptor(casadi::GenericType::Dict *), 0) != -1) {
      if (m) *m=*mp;
      return true;
    }
    if (!PyDict_Check(p)) return false;
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    casadi::GenericType gt;
    while (PyDict_Next(p, &pos, &key, &value)) {
      if (!PyString_Check(key)) return false;
      if (!to_GenericType(value, &gt)) return false;
      if (m) (*m)[std::string(PyString_AsString(key))] = gt;
    }
    return true;
#endif // SWIGPYTHON
  }
 }
%casadi_typemaps_constref(Dict, PRECEDENCE_DICT, casadi::GenericType::Dict)
#endif

%fragment("to"{SX}, "header", fragment="fwd") {
  // Traits specialization for SX
  namespace casadi {
    template <> struct traits_casptr<SX> {
      static bool casptr(GUESTOBJECT *p, SX** m) {
        // Treat Null
        if (is_null(p)) return false;

        // SX already?
        if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(m),
                                      $descriptor(casadi::Matrix<casadi::SXElement>*), 0))) {
          return true;
        }

        // Object is an DMatrix
        {
          // Pointer to object
          DMatrix *m2;
          if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(&m2),
                                        $descriptor(casadi::Matrix<double>*), 0))) {
            if (m) **m=*m2;
            return true;
          }
        }

        // Object is an IMatrix
        {
          // Pointer to object
          IMatrix *m2;
          if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(&m2),
                                        $descriptor(casadi::Matrix<int>*), 0))) {
            if (m) **m=*m2;
            return true;
          }
        }

        // First convert to double
        if (false) {
          double tmp;
          if (SWIG_IsOK(SWIG_AsVal_double(p, m ? &tmp : 0))) {
            if (m) **m=tmp;
            return true;
          }
        }

        // First convert to integer
        if (false) {
          int tmp;
          if (SWIG_IsOK(SWIG_AsVal_int(p, m ? &tmp : 0))) {
            if (m) **m=tmp;
            return true;
          }
        }

        // Try first converting to a temporary DMatrix
        {
          DMatrix tmp, *mt=&tmp;
          if(casadi::casptr(p, m ? &mt : 0)) {
            if (m) **m = *mt;
            return true;
          }
        }

#ifdef SWIGPYTHON
        // Numpy arrays will be cast to dense SX
        if (is_array(p)) {
          if (array_type(p) != NPY_OBJECT) return false;
          if (array_numdims(p)>2 || array_numdims(p)<1) return false;
          int nrows = array_size(p,0); // 1D array is cast into column vector
          int ncols  = array_numdims(p)==2 ? array_size(p,1) : 1;
          PyArrayIterObject* it = (PyArrayIterObject*)PyArray_IterNew(p);
          casadi::SX mT;
          if (m) mT = casadi::SX::zeros(ncols, nrows);
          int k=0;
          casadi::SX tmp;
          PyObject *pe;
          while (it->index < it->size) { 
            pe = *((PyObject**) PyArray_ITER_DATA(it));
            if (!to_SX(pe, &tmp) || !tmp.isScalar()) {
              Py_DECREF(it);
              return false;
            }
            if (m) mT(k++) = tmp;
            PyArray_ITER_NEXT(it);
          }
          Py_DECREF(it);
          if (m) **m = mT.T();
          return true;
        }
        // Object has __SX__ method
        if (PyObject_HasAttrString(p,"__SX__")) {
          char cmd[] = "__SX__";
          PyObject *cr = PyObject_CallMethod(p, cmd, 0);
          if (!cr) return false;
          int flag = casptr(cr, m);
          Py_DECREF(cr);
          return flag;
        }
#endif // SWIGPYTHON

        // No match
        return false;
      }
    };
  } // namespace casadi

  int to_SX(GUESTOBJECT *p, void *mv, int offs) {
    casadi::SX *m = static_cast<casadi::SX*>(mv);
    if (m) m += offs;

    // Call refactored version
    casadi::SX *m_orig = m;
    if (casptr(p, m ? &m : 0)) {
      if (m!=m_orig) *m_orig=*m;
      return true;
    }

    return false;
  }
 }
%casadi_typemaps_constref(SX, PRECEDENCE_SX, casadi::Matrix<casadi::SXElement>)
%casadi_typemaps_genericmatrix(SX, PRECEDENCE_SX, casadi::Matrix<casadi::SXElement>)
%casadi_typemaps_vector(SX, PRECEDENCE_SXVector, casadi::Matrix<casadi::SXElement>)
%casadi_typemaps_vector2(SX, PRECEDENCE_SXVectorVector, casadi::Matrix<casadi::SXElement>)

%fragment("to"{MX}, "header", fragment="fwd") {
  // Traits specialization for MX
  namespace casadi {
    template <> struct traits_casptr<MX> {
      static bool casptr(GUESTOBJECT *p, MX** m) {
        // Treat Null
        if (is_null(p)) return false;

        // MX already?
        if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(m),
                                      $descriptor(casadi::MX*), 0))) {
          return true;
        }

        // Object is an DMatrix
        {
          // Pointer to object
          DMatrix *m2;
          if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(&m2),
                                        $descriptor(casadi::Matrix<double>*), 0))) {
            if (m) **m=*m2;
            return true;
          }
        }

        // Try first converting to a temporary DMatrix
        {
          DMatrix tmp, *mt=&tmp;
          if(casadi::casptr(p, m ? &mt : 0)) {
            if (m) **m = *mt;
            return true;
          }
        }

#ifdef SWIGPYTHON
        if (PyObject_HasAttrString(p,"__MX__")) {
          char cmd[] = "__MX__";
          PyObject *cr = PyObject_CallMethod(p, cmd, 0);
          if (!cr) return false;
          int flag = casptr(cr, m);
          Py_DECREF(cr);
          return flag;
        }
#endif // SWIGPYTHON

        // No match
        return false;
      }
    };
  } // namespace casadi

  int to_MX(GUESTOBJECT *p, void *mv, int offs) {
    MX *m = static_cast<MX*>(mv);
    if (m) m += offs;

    // Call refactored version
    casadi::MX *m_orig = m;
    if (casptr(p, m ? &m : 0)) {
      if (m!=m_orig) *m_orig=*m;
      return true;
    }

    // Failure if reached this point
    return false;
  }
 }
%casadi_typemaps_constref(MX, PRECEDENCE_MX, casadi::MX)
%casadi_typemaps_genericmatrix(MX, PRECEDENCE_MX, casadi::MX)

 /* Maps are treated as dictionaries in the target language
    First instantiate the templates (no proxy classes needed)
  */
%template() std::map<std::string, casadi::MX >;
%template() std::map<std::string, casadi::Matrix<casadi::SXElement> >;
%template() std::map<std::string, casadi::Matrix<double> >;
%template() std::map<std::string, casadi::Sparsity >;
%template() std::map<std::string, std::vector<casadi::Sparsity > >;

/* Map conversion from dictionaries (default routines handle opposite direction) */
%define %typemaps_map(xName, xPrec, xType...)
%casadi_typecheck_typemap_constref(Map, xPrec, std::map<std::string, xType >)
%casadi_in_typemap_constref2(Map, std::map<std::string, xType >)
%casadi_freearg_typemap(const std::map<std::string, xType >&)
%enddef
%typemaps_map(Map, PRECEDENCE_MX, casadi::MX)
%typemaps_map(Map, PRECEDENCE_DMatrix, casadi::Matrix<double>)
%typemaps_map(Map, PRECEDENCE_SX, casadi::Matrix<casadi::SXElement>)

 /* Pair of above maps and vector of strings
  */
%template() std::pair<std::map<std::string, casadi::MX >, std::vector<std::string> >;
%template() std::pair<std::map<std::string, casadi::Matrix<casadi::SXElement> >, std::vector<std::string> >;
%template() std::pair<std::map<std::string, casadi::Matrix<double> >, std::vector<std::string> >;
%template() std::pair<std::map<std::string, casadi::Sparsity >, std::vector<std::string> >;

/* Pair of Map conversion from dictionaries (default routines handle opposite direction) */
%define %typemaps_pairmap(xPrec, xType...)
%casadi_typecheck_typemap_constref(PairMap, xPrec, std::pair<std::map<std::string, xType >, std::vector<std::string> >)
%casadi_in_typemap_constref2(PairMap, std::pair<std::map<std::string, xType >, std::vector<std::string> >)
%casadi_freearg_typemap(const std::pair<std::map<std::string, xType >, std::vector<std::string> >&)
%enddef
%typemaps_pairmap(PRECEDENCE_MX, casadi::MX)
%typemaps_pairmap(PRECEDENCE_DMatrix, casadi::Matrix<double>)
%typemaps_pairmap(PRECEDENCE_SX, casadi::Matrix<casadi::SXElement>)

%fragment("to"{DMatrix}, "header", fragment="fwd,make_vector", fragment="to"{double}) {
  // Traits specialization for DMatrix
  namespace casadi {
    template <> struct traits_casptr<DMatrix> {
      static bool casptr(GUESTOBJECT *p, DMatrix** m) {
        // Treat Null
        if (is_null(p)) return false;

        // DMatrix already?
        if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(m),
                                      $descriptor(casadi::Matrix<double>*), 0))) {
          return true;
        }

        // Object is an IMatrix
        {
          // Pointer to object
          IMatrix *m2;
          if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(&m2),
                                        $descriptor(casadi::Matrix<int>*), 0))) {
            if (m) **m=*m2;
            return true;
          }
        }

        // First convert to double
        {
          double tmp;
          if (SWIG_IsOK(SWIG_AsVal_double(p, m ? &tmp : 0))) {
            if (m) **m=tmp;
            return true;
          }
        }

        // First convert to integer
        {
          int tmp;
          if (SWIG_IsOK(SWIG_AsVal_int(p, m ? &tmp : 0))) {
            if (m) **m=tmp;
            return true;
          }
        }

        // First convert to a vector<double>
        if (false) {
          std::vector<double> tmp, *tmp_ptr=m ? &tmp : 0;
          if (SWIG_IsOK(swig::asval(p, tmp_ptr))) {
            if (m) **m = tmp.empty() ? DMatrix() : DMatrix(tmp);
            return true;
          }
        }

#ifdef SWIGPYTHON
    // Object has __DMatrix__ method
    if (PyObject_HasAttrString(p,"__DMatrix__")) {
      char name[] = "__DMatrix__";
      PyObject *cr = PyObject_CallMethod(p, name, 0);
      if (!cr) return false;
      int result = to_DMatrix(cr, m ? *m : 0);
      Py_DECREF(cr);
      return result;
    }
    // Numpy arrays will be cast to dense Matrix<double>
    if (is_array(p)) { 
      if (array_numdims(p)==0) {
        double d;
        int result = to_double(p, &d);
        if (!result) return result;
        if (m) **m = casadi::Matrix<double>(d);
        return result;
      }
      if (array_numdims(p)>2 || array_numdims(p)<1) {
        return false;
      }
      int nrows = array_size(p,0); // 1D array is cast into column vector
      int ncols  = 1;
      if (array_numdims(p)==2) ncols=array_size(p,1); 
      int size=nrows*ncols; // number of elements in the dense matrix
      if (!array_is_native(p)) return false;
      // Make sure we have a contigous array with double datatype
      int array_is_new_object;
      PyArrayObject* array = obj_to_array_contiguous_allow_conversion(p,NPY_DOUBLE,&array_is_new_object);
      if (!array) return false;
    
      double* d=(double*) array_data(array);
      std::vector<double> v(d,d+size);
    
      if (m) {
        **m = casadi::Matrix<double>::zeros(nrows,ncols);
        (*m)->set(d, true);
      }
           
      // Free memory
      if (array_is_new_object) Py_DECREF(array); 
      return true;
    }
    // scipy's csc_matrix will be cast to sparse DMatrix
    if(PyObjectHasClassName(p, "csc_matrix")) {
    
      // Get the dimensions of the csc_matrix
      PyObject * shape = PyObject_GetAttrString( p, "shape"); // need's to be decref'ed
      if (!shape) return false;
      if(!PyTuple_Check(shape) || PyTuple_Size(shape)!=2) {
        Py_DECREF(shape);
        return false;
      }
      int nrows=PyInt_AsLong(PyTuple_GetItem(shape,0));
      int ncols=PyInt_AsLong(PyTuple_GetItem(shape,1));
      Py_DECREF(shape);
  
      bool ret= false;
    
      PyObject * narray=0;
      PyObject * row=0;
      PyObject * colind=0;
      PyArrayObject* array=0;
      PyArrayObject* array_row=0;
      PyArrayObject* array_colind=0;
    
      int array_is_new_object=0;
      int row_is_new_object=0;
      int colind_is_new_object=0;
    
      // Fetch data
      narray=PyObject_GetAttrString( p, "data"); // need's to be decref'ed
      if (!narray || !is_array(narray) || array_numdims(narray)!=1) goto cleanup;
      array = obj_to_array_contiguous_allow_conversion(narray,NPY_DOUBLE,&array_is_new_object);
      if (!array) goto cleanup;

      // Construct the 'row' vector needed for initialising the correct sparsity
      row = PyObject_GetAttrString(p,"indices"); // need's to be decref'ed
      if (!row || !is_array(row) || array_numdims(row)!=1) goto cleanup;
      array_row = obj_to_array_contiguous_allow_conversion(row,NPY_INT,&row_is_new_object);
      if (!array_row) goto cleanup;

      // Construct the 'colind' vector needed for initialising the correct sparsity
      colind = PyObject_GetAttrString(p,"indptr"); // need's to be decref'ed
      if (!colind || !is_array(colind) || array_numdims(colind)!=1) goto cleanup;
      array_colind = obj_to_array_contiguous_allow_conversion(colind,NPY_INT,&colind_is_new_object);
      if (!array_colind) goto cleanup;
      {
        int size=array_size(array,0); // number on non-zeros
        double* d=(double*) array_data(array);
        std::vector<double> v(d,d+size);
      
        int* rowd=(int*) array_data(array_row);
        std::vector<int> rowv(rowd,rowd+size);
      
        int* colindd=(int*) array_data(array_colind);
        std::vector<int> colindv(colindd,colindd+(ncols+1));
      
        if (m) **m = casadi::Matrix<double>(casadi::Sparsity(nrows,ncols,colindv,rowv), v, false);
      
        ret = true;
      }
    
    cleanup: // yes that's right; goto.
      // Rather that than a pyramid of conditional memory-deallocation
      // TODO(jaeandersson): Create a helper struct and put the below in the destructor
      if (array_is_new_object && array) Py_DECREF(array);
      if (narray) Py_DECREF(narray);
      if (row_is_new_object && array_row) Py_DECREF(array_row);
      if (row) Py_DECREF(row);
      if (colind_is_new_object && array_colind) Py_DECREF(array_colind);
      if (colind) Py_DECREF(colind);
      return ret;
    }
    if(PyObject_HasAttrString(p,"tocsc")) {
      char name[] = "tocsc";
      PyObject *cr = PyObject_CallMethod(p, name,0);
      if (!cr) return false;
      int result = to_DMatrix(cr, m ? *m : 0);
      Py_DECREF(cr);
      return result;
    }


    /* { */
    /*   std::vector<double> tmp, *tmp_ptr=m ? &tmp : 0; */
    /*   if (SWIG_IsOK(swig::asval(p, tmp_ptr))) { */
    /*     if (m) *m = tmp.empty() ? DMatrix() : DMatrix(tmp); */
    /*     return true; */
    /*   } */
    /* } */


    {
      std::vector <double> t;
      int res = make_vector(p, &t, to_double);
      if (t.size()>0) {
        if (m) **m = casadi::Matrix<double>(t);
      } else {
        if (m) **m = casadi::Matrix<double>(0,0);
      }
      return res;
    }
#endif // SWIGPYTHON
#ifdef SWIGMATLAB
    // MATLAB double matrix (sparse or dense)
    if (mxIsDouble(p) && mxGetNumberOfDimensions(p)==2) {
      if (m) {
        **m = casadi::DMatrix(getSparsity(p));
        double* data = static_cast<double*>(mxGetData(p));
        (*m)->setNZ(data);
      }
      return true;
    }
#endif // SWIGMATLAB

        // No match
        return false;
      }
    };
  } // namespace casadi

  int to_DMatrix(GUESTOBJECT *p, void *mv, int offs) {
    casadi::DMatrix *m = static_cast<casadi::DMatrix*>(mv);
    if (m) m += offs;

    // Call refactored version
    casadi::DMatrix *m_orig = m;
    if (casptr(p, m ? &m : 0)) {
      if (m!=m_orig) *m_orig=*m;
      return true;
    }

    return false;
  }
}
%casadi_typemaps_constref(DMatrix, PRECEDENCE_DMatrix, casadi::Matrix<double>)
%casadi_typemaps_genericmatrix(DMatrix, PRECEDENCE_DMatrix, casadi::Matrix<double>)
%casadi_typemaps_vector(MX, PRECEDENCE_MXVector, casadi::MX)

%fragment("to"{IMatrix}, "header", fragment="fwd,make_vector") {
  // Traits specialization for IMatrix
  namespace casadi {
    template <> struct traits_casptr<IMatrix> {
      static bool casptr(GUESTOBJECT *p, IMatrix** m) {
        // Treat Null
        if (is_null(p)) return false;

        // IMatrix already?
        if (SWIG_IsOK(SWIG_ConvertPtr(p, reinterpret_cast<void**>(m),
                                      $descriptor(casadi::Matrix<int>*), 0))) {
          return true;
        }

        // First convert to integer
        {
          int tmp;
          if (SWIG_IsOK(SWIG_AsVal_int(p, m ? &tmp : 0))) {
            if (m) **m=tmp;
            return true;
          }
        }

#ifdef SWIGPYTHON
        // User-supplied conversion
        if (false && PyObject_HasAttrString(p,"__IMatrix__")) {
          char cmd[] = "__IMatrix__";
          PyObject *cr = PyObject_CallMethod(p, cmd, 0);
          if (!cr) return false;
          // Call recursively
          int flag = casptr(cr, m);
          Py_DECREF(cr);
          return flag;
        }
#endif // SWIGPYTHON

        // First convert to DMatrix and check if all entries integers
        if (false) {
          DMatrix tmp;
          if (to_DMatrix(p, &tmp) && tmp.isInteger()) {
            if (m) **m=IMatrix(tmp.sparsity(), tmp.nonzeros_int());
            return true;
          }
        }

        // No match
        return false;
      }
    };
  } // namespace casadi


  int to_IMatrix(GUESTOBJECT *p, void *mv, int offs) {
    casadi::IMatrix *m = static_cast<casadi::IMatrix*>(mv);
    // Use operator=
    if (m) m += offs;

    // Call refactored version
    casadi::IMatrix *m_orig = m;
    if (casptr(p, m ? &m : 0)) {
      if (m!=m_orig) *m_orig=*m;
      return true;
    }

#ifdef SWIGPYTHON
    // Numpy arrays will be cast to dense Matrix<int>
    if (is_array(p)) {
      if (array_numdims(p)==0) {
        int d;
        int result = to_int(p, &d);
        if (!result) return result;
        if (m) *m = casadi::Matrix<int>(d);
        return result;
      }
      if (array_numdims(p)>2 || array_numdims(p)<1) return false;
      int nrows = array_size(p,0); // 1D array is cast into column vector
      int ncols  = 1;
      if (array_numdims(p)==2) ncols=array_size(p,1); 
      if (!array_is_native(p)) return false;
      // Make sure we have a contigous array with int datatype
      int array_is_new_object=0;
      PyArrayObject* array = obj_to_array_contiguous_allow_conversion(p,NPY_INT,&array_is_new_object);
      if (!array) { // Trying LONG
        PyErr_Clear();
        array = obj_to_array_contiguous_allow_conversion(p,NPY_LONG,&array_is_new_object);
        if (!array) {
          PyErr_Clear();
          return false;
        }
        long* temp=(long*) array_data(array);
        if (m) {
          *m = casadi::Matrix<int>::zeros(ncols,nrows);
          for (int k=0;k<nrows*ncols;k++) m->data()[k]=temp[k];
          *m = m->T();                  
        }
        // Free memory
        if (array_is_new_object) Py_DECREF(array);
        return true;
      }
    
      if (m) {
        int* d= static_cast<int*>(array_data(array));
        *m = casadi::Matrix<int>::zeros(nrows, ncols);
        for (int rr=0; rr<nrows; ++rr) {
          for (int cc=0; cc<ncols; ++cc) {
            (*m)(rr, cc) = *d++;
          }
        }
      }
           
      // Free memory
      if (array_is_new_object) Py_DECREF(array);
      return true;
    }
    if (PyObject_HasAttrString(p,"__IMatrix__")) {
      char name[] = "__IMatrix__";
      PyObject *cr = PyObject_CallMethod(p, name,0);
      if (!cr) { return false; }
      int result = to_IMatrix(cr, m);
      Py_DECREF(cr);
      return result;
    }
    {
      std::vector <int> t;
      int res = make_vector(p, &t, to_int);
      if (m) *m = casadi::Matrix<int>(t);
      return res;
    }
    return true;
#endif // SWIGPYTHON
#ifdef SWIGMATLAB
    // In MATLAB, it is common to use floating point values to represent integers
    if (mxIsDouble(p) && mxGetNumberOfDimensions(p)==2) {
      double* data = static_cast<double*>(mxGetData(p));

      // Check if all integers
      bool all_integers=true;
      size_t sz = getNNZ(p);
      for (size_t i=0; i<sz; ++i) {
        if (data[i] != int(data[i])) {
          all_integers = false;
          break;
        }
      }

      // If successful
      if (all_integers) {
        if (m) {
          *m = casadi::IMatrix(getSparsity(p));
          for (size_t i=0; i<sz; ++i) {
            m->at(i) = int(data[i]);
          }
        }
        return true;
      }
    }
#endif // SWIGMATLAB
    return false;
  }
 }
%casadi_typemaps_constref(IMatrix, PRECEDENCE_IMatrix, casadi::Matrix<int>)
%casadi_typemaps_genericmatrix(IMatrix, PRECEDENCE_IMatrix, casadi::Matrix<int>)
%casadi_typemaps_vector2(MX, PRECEDENCE_MXVectorVector, casadi::MX)
%casadi_typemaps_vector(DMatrix, PRECEDENCE_DMatrixVector, casadi::Matrix<double>)
%casadi_typemaps_vector(IMatrix, PRECEDENCE_IMatrixVector, casadi::Matrix<int>)
%casadi_typemaps_vector2(DMatrix, PRECEDENCE_DMatrixVectorVector, casadi::Matrix<double>)
%casadi_typemaps_vector2(IMatrix, PRECEDENCE_IMatrixVectorVector, casadi::Matrix<int>)
%casadi_typemaps_vector(IVector, PRECEDENCE_IVectorVector, std::vector<int>)

%define %my_value_output_typemaps(Type,...)
%value_output_typemap(%arg(swig::from), %arg(SWIG_Traits_frag(Type)), %arg(Type));
%enddef

// These make OUTPUT behave like expected for non std container types
%my_value_output_typemaps(casadi::Matrix< casadi::SXElement >);
%my_value_output_typemaps(casadi::Matrix< double >);
%my_value_output_typemaps(casadi::Matrix< int >);
%my_value_output_typemaps(casadi::MX);
%my_value_output_typemaps(casadi::Sparsity);

#ifdef SWIGPYTHON
%outputRefOwn(casadi::Sparsity)

%outputRefNew(std::vector< int >)
%outputRefNew(std::vector< double >)

%outputRefOwn(casadi::Matrix< double >)
%outputRefOwn(casadi::Matrix< casadi::SXElement >)
#endif // SWIGPYTHON

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

%template() std::map<std::string,casadi::GenericType>;

%traits_swigtype(casadi::Function);
%fragment(SWIG_Traits_frag(casadi::Function));

#endif


// These dummy things would go away when we properly use fragments
// %traits_swigtype

%{
namespace std {
void dummy(std::vector< std::vector<double> > foo1,
           std::vector<double> &foo2,
           std::vector<casadi::MX> &foo3,
           casadi::MX foo4,
           casadi::Matrix<double> foo5,
           casadi::Sparsity foo6,
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
           casadi::Dict foo18,
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
  void dummy(std::vector< std::vector<double> > foo1,
             std::vector<double> &foo2,
             std::vector<casadi::MX> &foo3,
             casadi::MX foo4,
             casadi::Matrix<double> foo5,
             casadi::Sparsity foo6,
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
             casadi::Dict foo18,
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
%rename(row) getRow;
%rename(colind) getColind;
%rename(sparsity) getSparsity;

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
%rename(__float__) getValue;
%rename(__int__) getIntValue;

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
            return self.get(False, s[0], s[1])
          return self.get(False, s)

    def __setitem__(self,s,val):
        with internalAPI():
          if isinstance(s,tuple) and len(s)==2:
            return self.set(val, False, s[0], s[1])  
          return self.set(val, False, s)
        
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
    const Type getitem(const Slice& rr) const { Type m; $self->get(m, true, rr); return m;}
    const Type getitem(const Matrix<int>& rr) const { Type m; $self->get(m, true, rr); return m;}
    const Type getitem(const Sparsity& sp) const { Type m; $self->get(m, true, sp); return m;}
    const Type getitem(const Slice& rr, const Slice& cc) const { Type m; $self->get(m, true, rr, cc); return m;}
    const Type getitem(const Slice& rr, const Matrix<int>& cc) const { Type m; $self->get(m, true, rr, cc); return m;}
    const Type getitem(const Matrix<int>& rr, const Slice& cc) const { Type m; $self->get(m, true, rr, cc); return m;}
    const Type getitem(const Matrix<int>& rr, const Matrix<int>& cc) const { Type m; $self->get(m, true, rr, cc); return m;}

    // Set a submatrix (index-1)
    void setitem(const Type& m, const Slice& rr) { $self->set(m, true, rr);}
    void setitem(const Type& m, const Matrix<int>& rr) { $self->set(m, true, rr);}
    void setitem(const Type& m, const Sparsity& sp) { $self->set(m, true, sp);}
    void setitem(const Type& m, const Slice& rr, const Slice& cc) { $self->set(m, true, rr, cc);}
    void setitem(const Type& m, const Slice& rr, const Matrix<int>& cc) { $self->set(m, true, rr, cc);}
    void setitem(const Type& m, const Matrix<int>& rr, const Slice& cc) { $self->set(m, true, rr, cc);}
    void setitem(const Type& m, const Matrix<int>& rr, const Matrix<int>& cc) { $self->set(m, true, rr, cc);}

    // Get nonzeros (index-1)
    const Type getitemcurl(const Slice& rr) const { Type m; $self->getNZ(m, true, rr); return m;}
    const Type getitemcurl(const Matrix<int>& rr) const { Type m; $self->getNZ(m, true, rr); return m;}

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

%include <casadi/core/printable_object.hpp>

#ifdef SWIGPYTHON
%rename(SWIG_STR) getDescription;
#endif // SWIGPYTHON

%template(PrintSharedObject) casadi::PrintableObject<casadi::SharedObject>;
%template(PrintSlice)        casadi::PrintableObject<casadi::Slice>;
%template(PrintIMatrix)      casadi::PrintableObject<casadi::Matrix<int> >;
%template(PrintDMatrix)      casadi::PrintableObject<casadi::Matrix<double> >;
//%template(PrintSX)           casadi::PrintableObject<casadi::Matrix<casadi::SXElement> >;
%template(PrintNlpBuilder)     casadi::PrintableObject<casadi::NlpBuilder>;
%template(PrintVariable)        casadi::PrintableObject<casadi::Variable>;
%template(PrintDaeBuilder)     casadi::PrintableObject<casadi::DaeBuilder>;

%include <casadi/core/shared_object.hpp>
%include <casadi/core/std_vector_tools.hpp>

VECTOR_TOOLS_TEMPLATES(int)
VECTOR_TOOLS_TEMPLATES(double)
%define VECTOR_REPR(type)
%extend std::vector< type >{
  std::string SWIG_REPR(){ return casadi::getRepresentation(*$self); }
  std::string SWIG_STR(){ return casadi::getDescription(*$self); }
};
%enddef
%include <casadi/core/weak_ref.hpp>
%include <casadi/core/casadi_types.hpp>
%include <casadi/core/generic_type.hpp>

%casadi_typemaps_constref(GenericType, PRECEDENCE_GENERICTYPE, casadi::GenericType)
%casadi_typemaps_vector(GenericType, PRECEDENCE_GENERICTYPE, casadi::GenericType)

%include <casadi/core/options_functionality.hpp>

namespace casadi {
  %extend OptionsFunctionality {
    void setOption(const std::string &name, const std::string& val){$self->setOption(name,val);} 
    void setOption(const std::string &name, const std::vector<int>& val){$self->setOption(name,val);} 
    void setOption(const std::string &name, const std::vector<double>& val){$self->setOption(name,val);} 
    void setOption(const std::string &name, double val){$self->setOption(name,val);}
    void setOption(const std::string &name, int val){$self->setOption(name,val);} 
    void setOption(const std::string &name, bool val){$self->setOption(name,val);}  
    void setOption(const std::string &name, const std::vector< std::vector<int> >& val){$self->setOption(name,val);}  
  }
} // namespace casadi

%include <casadi/core/casadi_calculus.hpp>

%include <casadi/core/matrix/sparsity_interface.hpp>

%template(SpSparsity) casadi::SparsityInterface<casadi::Sparsity>;
%include <casadi/core/matrix/sparsity.hpp>

// Logic for pickling
#ifdef SWIGPYTHON
namespace casadi{
%extend Sparsity {
  %pythoncode %{
    def __setstate__(self, state):
        if state:
          self.__init__(state["nrow"],state["ncol"],state["colind"],state["row"])
        else:
          self.__init__()

    def __getstate__(self):
        if self.isNull(): return {}
        return {"nrow": self.size1(), "ncol": self.size2(), "colind": numpy.array(self.colind(),dtype=int), "row": numpy.array(self.row(),dtype=int)}
  %}  
}

} // namespace casadi
#endif // SWIGPYTHON

%include <casadi/core/matrix/slice.hpp>

namespace casadi{

#ifndef SWIGXML

%fragment("to"{Slice}, "header", fragment="fwd") {
   int to_Slice(GUESTOBJECT *p, void *mv, int offs) {
    casadi::Slice *m = static_cast<casadi::Slice*>(mv);
    if (m) m += offs;
    // CasADi-Slice already
    if (is_a(p, $descriptor(casadi::Slice *))) {
      casadi::Slice *mp;
      if (SWIG_ConvertPtr(p, (void **) &mp, $descriptor(casadi::Slice *), 0) == -1)
        return false;
      if(m) *m = *mp;
      return true;
    }
#ifdef SWIGPYTHON
    // Python int
    if (PyInt_Check(p)) {
      if (m) {
        m->start_ = PyInt_AsLong(p);
        m->stop_ = m->start_+1;
        if (m->stop_==0) m->stop_ = std::numeric_limits<int>::max();
      }
      return true;
    }
    // Python slice
    if (PySlice_Check(p)) {
      PySliceObject *r = (PySliceObject*)(p);
      if (m) {
        m->start_ = (r->start == Py_None || PyInt_AsLong(r->start) < std::numeric_limits<int>::min()) 
          ? std::numeric_limits<int>::min() : PyInt_AsLong(r->start);
        m->stop_  = (r->stop ==Py_None || PyInt_AsLong(r->stop)> std::numeric_limits<int>::max())
          ? std::numeric_limits<int>::max() : PyInt_AsLong(r->stop);
        if(r->step !=Py_None) m->step_  = PyInt_AsLong(r->step);
      }
      return true;
    }
#endif // SWIGPYTHON
#ifdef SWIGMATLAB
    if (mxIsChar(p) && mxGetM(p)==1 && mxGetN(p)==1) {
      char ch;
      if(mxGetString(p, &ch,(mwSize)sizeof(&ch))) return SWIG_TypeError;
      if (ch==':') {
        if (m) *m = casadi::Slice();
        return true;
      }
    }
#endif // SWIGMATLAB
    // Failure if reached this point
    return false;
  }
}
%casadi_typemaps_constref(Slice, PRECEDENCE_SLICE, casadi::Slice)

#endif // SWIGXML

} // namespace casadi
%include <casadi/core/matrix/generic_expression_tools.hpp>

// map the template name to the instantiated name
%define GET_INST(DataType, function_name)
%template(function_name) casadi::function_name< DataType >;
%enddef

// Define template instantiations
%define GENERIC_EXPRESSION_TOOLS_TEMPLATES(DataType)
GET_INST(DataType, logic_and)
GET_INST(DataType, logic_or)
GET_INST(DataType, logic_not)
%enddef

#ifndef SWIGMATLAB
GENERIC_EXPRESSION_TOOLS_TEMPLATES(int)
GENERIC_EXPRESSION_TOOLS_TEMPLATES(double)
GENERIC_EXPRESSION_TOOLS_TEMPLATES(casadi::Matrix<int>)
GENERIC_EXPRESSION_TOOLS_TEMPLATES(casadi::Matrix<double>)
GENERIC_EXPRESSION_TOOLS_TEMPLATES(casadi::SX)
GENERIC_EXPRESSION_TOOLS_TEMPLATES(casadi::MX)
#endif // SWIGMATLAB

%template(SpIMatrix)        casadi::SparsityInterface<casadi::Matrix<int> >;
%template(SpDMatrix)        casadi::SparsityInterface<casadi::Matrix<double> >;
%template(SpSX)             casadi::SparsityInterface<casadi::Matrix<casadi::SXElement> >;
%template(SpMX)             casadi::SparsityInterface<casadi::MX>;


%include <casadi/core/matrix/generic_matrix.hpp>

%template(GenIMatrix)        casadi::GenericMatrix<casadi::Matrix<int> >;
%template(GenDMatrix)        casadi::GenericMatrix<casadi::Matrix<double> >;
%template(GenSX)             casadi::GenericMatrix<casadi::Matrix<casadi::SXElement> >;
%template(GenMX)             casadi::GenericMatrix<casadi::MX>;

namespace casadi{

#ifdef SWIGPYTHON

/**

Accepts: 2D numpy.ndarray, numpy.matrix (contiguous, native byte order, datatype double)   - DENSE
         1D numpy.ndarray, numpy.matrix (contiguous, native byte order, datatype double)   - SPARSE
         2D scipy.csc_matrix
*/

  %typemap(in,numinputs=1) (double * val,int len,int stride1, int stride2,SparsityType sp)  {
    PyObject* p = $input;
    $3 = 0;
    $4 = 0;
    if (is_array(p)) {
      if (!(array_is_native(p) && array_type(p)==NPY_DOUBLE))
        SWIG_exception_fail(SWIG_TypeError, "Array should be native & of datatype double");
			  
      if (!(array_is_contiguous(p))) {
        if (PyArray_CHKFLAGS((PyArrayObject *) p,NPY_ALIGNED)) {
          $3 = PyArray_STRIDE((PyArrayObject *) p,0)/sizeof(double);
          $4 = PyArray_STRIDE((PyArrayObject *) p,1)/sizeof(double);
        } else {
          SWIG_exception_fail(SWIG_TypeError, "Array should be contiguous or aligned");
        }
      }
	    
      if (!(array_size(p,0)==arg1->size1() && array_size(p,1)==arg1->size2()) ) {
        std::stringstream s;
        s << "SWIG::typemap(in) (double *val,int len,SparsityType sp) " << std::endl;
        s << "Array is not of correct shape.";
        s << "Expecting shape (" << arg1->size1() << "," << arg1->size2() << ")" << ", but got shape ("
          << array_size(p,0) << "," << array_size(p,1) <<") instead.";
        const std::string tmp(s.str());
        const char* cstr = tmp.c_str();
        SWIG_exception_fail(SWIG_TypeError,  cstr);
      }
      $5 = casadi::SP_DENSETRANS;
      $2 = array_size(p,0)*array_size(p,1);
      $1 = (double*) array_data(p);
    } else if (PyObjectHasClassName(p,"csc_matrix")) {
      $5 = casadi::SP_SPARSE;
      PyObject * narray=PyObject_GetAttrString( p, "data"); // narray needs to be decref'ed
      if (!(array_is_contiguous(narray) && array_is_native(narray) && array_type(narray)==NPY_DOUBLE))
        SWIG_exception_fail(SWIG_TypeError, "csc_matrix should be contiguous, native & of datatype double");
      $2 = array_size(narray,0);
      if (!(array_size(narray,0)==arg1->nnz() ) ) {
        std::stringstream s;
        s << "SWIG::typemap(in) (double *val,int len,SparsityType sp) " << std::endl;
        s << "csc_matrix does not have correct number of non-zero elements.";
        s << "Expecting " << arg1->nnz() << " non-zeros, but got " << array_size(narray,0) << " instead.";
        const std::string tmp(s.str());
        const char* cstr = tmp.c_str();
        Py_DECREF(narray);
        SWIG_exception_fail(SWIG_TypeError,  cstr);
      }
      $1 = (double*) array_data(narray);
      Py_DECREF(narray);
    } else {
      SWIG_exception_fail(SWIG_TypeError, "Unrecognised object");
    }
  }

%typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) (double * val,int len,int stride1, int stride2,SparsityType sp) {
  PyObject* p = $input;
  if (((is_array(p) && array_numdims(p) == 2)  && array_type(p)!=NPY_OBJECT) || PyObjectHasClassName(p,"csc_matrix")) {
    $1=1;
  } else {
    $1=0;
  }
}

%typemap(in,numinputs=1) (double * val,int len,int stride1, int stride2)  {
    PyObject* p = $input;
    $3 = 0;
    $4 = 0;
    if (!(array_is_native(p) && array_type(p)==NPY_DOUBLE))
      SWIG_exception_fail(SWIG_TypeError, "Array should be native & of datatype double");
			  
    if (!(array_is_contiguous(p))) {
      if (PyArray_CHKFLAGS((PyArrayObject *) p,NPY_ALIGNED)) {
        $3 = PyArray_STRIDE((PyArrayObject *) p,0)/sizeof(double);
        $4 = PyArray_STRIDE((PyArrayObject *) p,1)/sizeof(double);
      } else {
        SWIG_exception_fail(SWIG_TypeError, "Array should be contiguous or aligned");
      }
    }
	    
    if (!(array_size(p,0)==arg1->nnz()) ) {
      std::stringstream s;
      s << "SWIG::typemap(in) (double *val,int len,SparsityType sp) " << std::endl;
      s << "Array is not of correct size. Should match number of non-zero elements.";
      s << "Expecting " << array_size(p,0) << " non-zeros, but got " << arg1->nnz() <<" instead.";
      const std::string tmp(s.str());
      const char* cstr = tmp.c_str();
      SWIG_exception_fail(SWIG_TypeError,  cstr);
    }
    $2 = array_size(p,0);
    $1 = (double*) array_data(p);
  }

%typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) (double * val,int len,int stride1, int stride2) {
  PyObject* p = $input;
  if ((is_array(p) && array_numdims(p) == 1) && array_type(p)!=NPY_OBJECT) {
    $1=1;
  } else {
    $1=0;
  }
}

#endif // SWIGPYTHON

} // namespace casadi

%include <casadi/core/matrix/generic_expression.hpp>

%template(ExpIMatrix)        casadi::GenericExpression<casadi::Matrix<int> >;
%template(ExpDMatrix)        casadi::GenericExpression<casadi::Matrix<double> >;
%template(ExpSX)             casadi::GenericExpression<casadi::Matrix<casadi::SXElement> >;
%template(ExpMX)             casadi::GenericExpression<casadi::MX>;

// FIXME: Placing in printable_object.i does not work
%template(PrintSX)           casadi::PrintableObject<casadi::Matrix<casadi::SXElement> >;

%include <casadi/core/matrix/matrix.hpp>

%template(IMatrix)           casadi::Matrix<int>;
%template(DMatrix)           casadi::Matrix<double>;

%extend casadi::Matrix<double> {
   %template(DMatrix) Matrix<int>;
};

namespace casadi{
  %extend Matrix<double> {

    void assign(const casadi::Matrix<double>&rhs) { (*$self)=rhs; }
    %matrix_convertors
    %matrix_helpers(casadi::Matrix<double>)

  }
  %extend Matrix<int> {

    void assign(const casadi::Matrix<int>&rhs) { (*$self)=rhs; }
    %matrix_convertors
    %matrix_helpers(casadi::Matrix<int>)

  }
}

// Extend DMatrix with SWIG unique features
namespace casadi{
  %extend Matrix<double> {
    // Convert to a dense matrix
    GUESTOBJECT* full() const {
#ifdef SWIGPYTHON
      npy_intp dims[2] = {$self->size1(), $self->size2()};
      PyObject* ret = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
      double* d = static_cast<double*>(array_data(ret));
      $self->get(d, true); // Row-major
      return ret;
#elif defined(SWIGMATLAB)
      mxArray *p  = mxCreateDoubleMatrix($self->size1(), $self->size2(), mxREAL);
      double* d = static_cast<double*>(mxGetData(p));
      $self->get(d); // Column-major
      return p;
#else
      return 0;
#endif
    }

#ifdef SWIGMATLAB
    // Convert to a sparse matrix
    GUESTOBJECT* zz_sparse() const {
      mxArray *p  = mxCreateSparse($self->size1(), $self->size2(), $self->nnz(), mxREAL);
      $self->getNZ(static_cast<double*>(mxGetData(p)));
      std::copy($self->colind(), $self->colind()+$self->size2()+1, mxGetJc(p));
      std::copy($self->row(), $self->row()+$self->size2()+1, mxGetIr(p));
      return p;
    }
#endif
  }
} // namespace casadi


#ifdef SWIGPYTHON
namespace casadi{
%extend Matrix<double> {
/// Create a 2D contiguous NP_DOUBLE numpy.ndarray

PyObject* arrayView() {
  if ($self->nnz()!=$self->numel()) 
    throw  casadi::CasadiException("Matrix<double>::arrayview() can only construct arrayviews for dense DMatrices.");
  npy_intp dims[2];
  dims[0] = $self->size2();
  dims[1] = $self->size1();
  std::vector<double> &v = $self->data();
  PyArrayObject* temp = (PyArrayObject*) PyArray_New(&PyArray_Type, 2, dims, NPY_DOUBLE, NULL, &v[0], 0, NPY_ARRAY_CARRAY, NULL);
  PyObject* ret = PyArray_Transpose(temp,NULL);
  Py_DECREF(temp); 
  return ret;
}
    
%pythoncode %{
  def toArray(self,shared=False):
    import numpy as n
    if shared:
      if not self.isDense():
        raise Expection("toArray(shared=True) only possible for dense arrays.")
      return self.arrayView()
    else:
      if isinstance(self,IMatrix):
        return n.array(self.get(True),n.int).reshape(self.shape)
      else:    
        return n.array(self.get(True)).reshape(self.shape)
%}

%python_array_wrappers(999.0)

// The following code has some trickery to fool numpy ufunc.
// Normally, because of the presence of __array__, an ufunctor like nump.sqrt
// will unleash its activity on the output of __array__
// However, we wish DMatrix to remain a DMatrix
// So when we receive a call from a functor, we return a dummy empty array
// and return the real result during the postprocessing (__array_wrap__) of the functor.
%pythoncode %{
  def __array_custom__(self,*args,**kwargs):
    if "dtype" in kwargs and not(isinstance(kwargs["dtype"],n.double)):
      return n.array(self.toArray(),dtype=kwargs["dtype"])
    else:
      return self.toArray()
%}

%pythoncode %{
  def toCsc_matrix(self):
    import numpy as n
    import warnings
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      from scipy.sparse import csc_matrix
    return csc_matrix( (self.nonzeros(),self.row(),self.colind()), shape = self.shape, dtype=n.double )

  def tocsc(self):
    return self.toCsc_matrix()

%}

%pythoncode %{
  def __nonzero__(self):
    if self.numel()!=1:
      raise Exception("Only a scalar can be cast to a float")
    if self.nnz()==0:
      return 0
    return float(self)!=0
%}

%pythoncode %{
  def __abs__(self):
    return abs(float(self))
%}

#ifdef SWIGPYTHON
binopsrFull(casadi::Matrix<double>)
binopsFull(const casadi::SX & b,,casadi::SX,casadi::SX)
binopsFull(const casadi::MX & b,,casadi::MX,casadi::MX)
#endif // SWIGPYTHON

}; // extend Matrix<double>

%extend Matrix<int> {

  %python_array_wrappers(998.0)

  %pythoncode %{
    def toArray(self):
      import numpy as n
      r = n.zeros((self.size1(),self.size2()))
      d = self.nonzeros_int()
      for j in range(self.size2()):
        for k in range(self.colind(j),self.colind(j+1)):
          i = self.row(k)
          r[i,j] = d[k]
      return r
  %}
  
#ifdef SWIGPYTHON
  binopsrFull(casadi::Matrix<int>)
  binopsFull(const casadi::SX & b,,casadi::SX,casadi::SX)
  binopsFull(const casadi::Matrix<double> & b,,casadi::Matrix<double>,casadi::Matrix<double>)
  binopsFull(const casadi::MX & b,,casadi::MX,casadi::MX)
#endif // SWIGPYTHON
  %pythoncode %{
    def __abs__(self):
      return abs(int(self))
  %}
} // extend Matrix<int>


// Logic for pickling

%extend Matrix<int> {

  %pythoncode %{
    def __setstate__(self, state):
        sp = Sparsity.__new__(Sparsity)
        sp.__setstate__(state["sparsity"])
        self.__init__(sp,state["data"])

    def __getstate__(self):
        return {"sparsity" : self.sparsity().__getstate__(), "data": numpy.array(self.nonzeros_int(),dtype=int)}
  %} 
}

%extend Matrix<double> {

  %pythoncode %{
    def __setstate__(self, state):
        sp = Sparsity.__new__(Sparsity)
        sp.__setstate__(state["sparsity"])
        self.__init__(sp,state["data"])

    def __getstate__(self):
        return {"sparsity" : self.sparsity().__getstate__(), "data": numpy.array(self.nonzeros(),dtype=float)}
  %}
  
}


} // namespace casadi
#endif // SWIGPYTHON

%include <casadi/core/sx/sx_element.hpp>

#ifdef SWIGPYTHON
%extend casadi::Sparsity{
    %pythoncode %{
        @property
        def shape(self):
            return (self.size1(),self.size2())
            
        @property
        def T(self):
            return _casadi.transpose(self)
            
        def __array__(self,*args,**kwargs):
            return DMatrix.ones(self).toArray()
    %}
};

#endif // SWIGPYTHON

VECTOR_REPR(casadi::Matrix<casadi::SXElement>)

#ifdef SWIGPYTHON
%pythoncode %{

try:
  import numpy

  def constpow(x,y):
    pass

  constpow=numpy.frompyfunc(constpow,2,1)

  fmin_backup = fmin
  fmax_backup = fmax

  def fmin(x,y):
    pass
    
  def fmax(x,y):
    pass

  _min_ufunc = numpy.frompyfunc(fmin,2,1)
  _max_ufunc = numpy.frompyfunc(fmax,2,1)

  fmin = fmin_backup
  fmax = fmax_backup

  _defaultmin = min
  def min(*args,**kwargs):
    if len(args)==2 and len(kwargs)==0 and (hasattr(args[0],'fmin') or hasattr(args[1],'fmin')):
      return _min_ufunc(*args)
    else:
      return _defaultmin(*args,**kwargs)
      
  _defaultmax = max
  def max(*args,**kwargs):
    if len(args)==2 and len(kwargs)==0 and (hasattr(args[0],'fmax') or hasattr(args[1],'fmax')):
      return _max_ufunc(*args)
    else:
      return _defaultmax(*args,**kwargs)
except:
  pass
%}
#endif // SWIGPYTHON

namespace casadi {
%extend Matrix<SXElement>{
    
    %matrix_convertors
    %matrix_helpers(casadi::Matrix<casadi::SXElement>)
       
    #ifdef SWIGPYTHON
    %pythoncode %{
    def toArray(self):
      import numpy as n
      r = n.array((),dtype=object)
      r.resize(self.size1(),self.size2())
      for j in range(self.size2()):
        for el in range(self.colind(j),self.colind(j+1)):
          i=self.row(el)
          r[i,j] = self.nz[el]
      return r
    %}
    
  %python_array_wrappers(1001.0)
  #endif // SWIGPYTHON 
  
#ifdef SWIGPYTHON
  binopsrFull(casadi::Matrix<casadi::SXElement>)  
#endif // SWIGPYTHON

};

} // namespace casadi

#ifdef SWIGPYTHON
#include <arrayobject.h>

// Template instantiations
%template()    std::vector<PyObject*>;


#endif // SWIGPYTHON


%template(SX)             casadi::Matrix<casadi::SXElement>;

%extend casadi::Matrix<casadi::SXElement> {
   %template(SX) Matrix<int>;
   %template(SX) Matrix<double>;
};


%include <casadi/core/mx/mx.hpp>

// Template instantiations
%template(Pair_MX_MXVector) std::pair<casadi::MX, std::vector<casadi::MX> >;

%extend casadi::MX{
  
  %matrix_helpers(casadi::MX)
  
  #ifdef SWIGPYTHON
  %python_array_wrappers(1002.0)
  
  %pythoncode %{
  def __array_custom__(self,*args,**kwargs):
    import numpy as np
    if np.__version__=="1.8.1": #1083
      return np.array(np.nan)
    raise Exception("MX cannot be converted to an array. MX.__array__ purely exists to allow ufunc/numpy goodies")
    
  def __iter__(self):
    return self.nz.__iter__()
    
  %}
  #endif //SWIGPYTHON
  
#ifdef SWIGPYTHON
  binopsrFull(casadi::MX)
#endif // SWIGPYTHON
};

VECTOR_REPR(casadi::MX)
VECTOR_REPR(std::vector<casadi::MX>)

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
      res = f([f2.getInput(i) for i in range(f2.nIn())])
      if not isinstance(res,list):
        res = [res]
      for i in range(f2.nOut()):
        f2.setOutput(res[i],i)
    Fun = CustomFunction(fcustom,inputs,outputs)
    Fun.setOption("name","CustomFunction")
    return Fun
  return wrap
  
def PyFunction(obj,inputs,outputs):
    @pyevaluate
    def fcustom(f):
      obj.evaluate([f.input(i) for i in range(f.nIn())],[f.output(i) for i in range(f.nOut())])
      
    Fun = CustomFunction(fcustom,inputs,outputs)
    Fun.setOption("name","CustomFunction")
    if hasattr(obj,'getDerForward'):
      @pyderivativegenerator
      def derivativewrap(f,nfwd):
        return obj.getDerForward(f,nfwd)
      Fun.setOption("custom_forward",derivativewrap)
      
    if hasattr(obj,'getDerReverse'):
      @pyderivativegenerator
      def derivativewrap(f,adj):
        return obj.getDerReverse(f,adj)
      Fun.setOption("custom_reverse",derivativewrap)
      
    if hasattr(obj,'fwd'):
      @pyderivativegenerator
      def derivativewrapFwd(f,nfwd):
        num_in = f.nIn()
        num_out = f.nOut()
        
        @pyevaluate
        def der(f2):
          all_inputs = [f2.input(i) for i in range(f2.nIn())]
          all_outputs = [f2.output(i) for i in range(f2.nOut())]
          inputs=all_inputs[:num_in]
          outputs=all_inputs[num_in:num_in+num_out]
          fwd_seeds=zip(*[iter(all_inputs[num_in+num_out:])]*num_in)
          fwd_sens=zip(*[iter(all_outputs)]*num_out)
          obj.fwd(inputs,outputs,fwd_seeds,fwd_sens)
          
        DerFun = CustomFunction(der,inputs+outputs+nfwd*inputs,nfwd*outputs)
        DerFun.setOption("name","CustomFunction_derivative")
        DerFun.init()
        return DerFun
 
      Fun.setOption("custom_forward",derivativewrapFwd)
    
    if hasattr(obj,'adj'):
      @pyderivativegenerator
      def derivativewrapAdj(f,nadj):
        num_in = f.nIn()
        num_out = f.nOut()
        
        @pyevaluate
        def der(f2):
          all_inputs = [f2.input(i) for i in range(f2.nIn())]
          all_outputs = [f2.output(i) for i in range(f2.nOut())]
          inputs=all_inputs[:num_in]
          outputs=all_inputs[num_in:num_in+num_out]
          adj_seeds=zip(*[iter(all_inputs[num_in+num_out:])]*num_out)
          adj_sens=zip(*[iter(all_outputs)]*num_in)
          obj.adj(inputs,outputs,adj_seeds,adj_sens)
          
        DerFun = CustomFunction(der,inputs+outputs+nadj*outputs,nadj*inputs)
        DerFun.setOption("name","CustomFunction_derivative")
        DerFun.init()
        return DerFun
 
      Fun.setOption("custom_reverse",derivativewrapAdj)
    return Fun
  
%}
#endif

%include <casadi/core/function/io_interface.hpp>

%template(IOInterfaceFunction) casadi::IOInterface<casadi::Function>;

%extend casadi::IOInterface<casadi::Function> {
  casadi::Matrix<double> getInput(int iind=0) const             { static_cast<const casadi::Function*>($self)->assertInit(); return $self->input(iind);}
  casadi::Matrix<double> getInput(const std::string &iname) const             { return $self->input($self->inputIndex(iname)); }
  casadi::Matrix<double> getOutput(int oind=0) const            { static_cast<const casadi::Function*>($self)->assertInit(); return $self->output(oind);}
}

%include <casadi/core/function/io_scheme.hpp>

%include <casadi/core/function/function.hpp>
%feature("copyctor", "0") casadi::CodeGenerator;
%include <casadi/core/function/code_generator.hpp>
%template() std::pair<casadi::Function,casadi::Function>;

%include <casadi/core/matrix/matrix_tools.hpp>

// map the template name to the instantiated name
%define MTT_INST(DataType, function_name)
%template(function_name) casadi::function_name <DataType >;
%enddef

// Define template instantiations
%define MATRIX_TOOLS_TEMPLATES(DataType)
  MTT_INST(DataType, solve)
  MTT_INST(DataType, pinv)
%enddef

#ifndef SWIGMATLAB
MATRIX_TOOLS_TEMPLATES(int)
MATRIX_TOOLS_TEMPLATES(double)
MATRIX_TOOLS_TEMPLATES(casadi::SXElement)
#endif // SWIGMATLAB

%define SPARSITY_INTERFACE_DECL(MatType...)
MatType horzcat(const std::vector< MatType > &v);
std::vector<MatType > horzsplit(const MatType &v, const std::vector<int>& offset);
std::vector<MatType > horzsplit(const MatType &v, int incr=1);
MatType vertcat(const std::vector< MatType > &v);
std::vector<MatType > vertsplit(const MatType &v, const std::vector<int>& offset);
std::vector<MatType > vertsplit(const MatType &v, int incr=1);
MatType blockcat(const std::vector< std::vector<MatType > > &v);
std::vector< std::vector< MatType > >
blocksplit(const MatType& x, const std::vector<int>& vert_offset, const std::vector<int>& horz_offset);
std::vector< std::vector< MatType > > blocksplit(const MatType& x, int vert_incr=1, int horz_incr=1);
MatType diagcat(const std::vector< MatType > &v);
std::vector< MatType > diagsplit(const MatType& x, const std::vector<int>& output_offset1,
                                 const std::vector<int>& output_offset2);
std::vector< MatType > diagsplit(const MatType& x, const std::vector<int>& output_offset);
std::vector< MatType > diagsplit(const MatType& x, int incr=1);
std::vector< MatType > diagsplit(const MatType& x, int incr1, int incr2);
MatType veccat(const std::vector< MatType >& x);
MatType vecNZcat(const std::vector< MatType >& x);
MatType mul(const MatType &x, const MatType &y);
MatType mul(const MatType &x, const MatType &y, const MatType &z);
MatType mul(const std::vector< MatType > &args);
MatType transpose(const MatType &X);
MatType vec(const MatType& a);
MatType vecNZ(const MatType& a);
MatType reshape(const MatType& a, const Sparsity& sp);
MatType reshape(const MatType& a, int nrow, int ncol);
MatType reshape(const MatType& a, std::pair<int, int> rc);
int sprank(const MatType& A);
MatType triu(const MatType& a, bool includeDiagonal=true);
MatType tril(const MatType& a, bool includeDiagonal=true);
int norm_0_mul(const MatType &x, const MatType &y);
MatType kron(const MatType& a, const MatType& b);
MatType repmat(const MatType &A, int n, int m=1);
MatType repmat(const MatType &A, const std::pair<int, int>& rc);
%enddef

%define GENERIC_MATRIX_DECL(MatType...)
MatType quad_form(const MatType &X, const MatType &A);
MatType quad_form(const MatType &X);
MatType sum_square(const MatType &X);
MatType linspace(const MatType &a, const MatType &b, int nsteps);
MatType cross(const MatType &a, const MatType &b, int dim = -1);
MatType det(const MatType& A);
MatType inv(const MatType& A);
MatType trace(const MatType& a);
bool isEqual(const MatType& x, const MatType& y, int depth=0);
MatType tril2symm(const MatType &a);
MatType triu2symm(const MatType &a);
MatType norm_F(const MatType &x);
MatType norm_2(const MatType &x);
MatType norm_1(const MatType &x);
MatType norm_inf(const MatType &x);
MatType sumAll(const MatType &x);
MatType sumCols(const MatType &x);
MatType sumRows(const MatType &x);
MatType inner_prod(const MatType &x, const MatType &y);
MatType outer_prod(const MatType &x, const MatType &y);
MatType nullspace(const MatType& A);
MatType polyval(const MatType& p, const MatType& x);
MatType diag(const MatType &A);
MatType unite(const MatType& A, const MatType& B);
MatType densify(const MatType& x);
MatType simplify(const MatType &x);
MatType if_else(const MatType &cond, const MatType &if_true, const MatType &if_false,
                bool short_circuit=true);
MatType conditional(const MatType& ind, const std::vector< MatType > &x,
                    const MatType &x_default, bool short_circuit=true);
bool dependsOn(const MatType& f, const MatType &arg);
%enddef

%define MATRIX_DECL(MatType...)
MatType adj(const MatType& A);
MatType getMinor(const MatType &x, int i, int j);
MatType cofactor(const MatType &x, int i, int j);
void qr(const MatType& A, MatType& OUTPUT, MatType& OUTPUT);
//MatType all(const MatType &x);
//MatType any(const MatType &x);
MatType project(const MatType& A, const Sparsity& sp);
MatType sparsify(const MatType& A, double tol=0);
MatType norm_inf_mul(const MatType &x, const MatType &y);
%enddef

%define GENERIC_MATRIX_TOOLS_TEMPLATES(MatType...)
SPARSITY_INTERFACE_DECL(MatType)
GENERIC_MATRIX_DECL(MatType)
%enddef

%define GENERIC_MATRIX_TOOLS_TEMPLATES_MATRIX(DataType...)
GENERIC_MATRIX_TOOLS_TEMPLATES(casadi::Matrix<DataType>)
MATRIX_DECL(casadi::Matrix<DataType>)
%enddef

#ifndef SWIGMATLAB
SPARSITY_INTERFACE_DECL(casadi::Sparsity)
GENERIC_MATRIX_TOOLS_TEMPLATES_MATRIX(int)
GENERIC_MATRIX_TOOLS_TEMPLATES_MATRIX(double)
GENERIC_MATRIX_TOOLS_TEMPLATES_MATRIX(casadi::SXElement)
GENERIC_MATRIX_TOOLS_TEMPLATES(casadi::MX)
#endif // SWIGMATLAB

%include <casadi/core/sx/sx_tools.hpp>
%include <casadi/core/mx/mx_tools.hpp>

%include <casadi/core/function/sx_function.hpp>
%include <casadi/core/function/mx_function.hpp>
%include <casadi/core/function/linear_solver.hpp>
%include <casadi/core/function/implicit_function.hpp>
%include <casadi/core/function/integrator.hpp>
%template() std::vector<casadi::Integrator>;
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
%include <casadi/core/function/switch.hpp>
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
%include <casadi/core/misc/nlp_builder.hpp>
%include <casadi/core/misc/variable.hpp>
%include <casadi/core/misc/dae_builder.hpp>
%include <casadi/core/misc/xml_file.hpp>
