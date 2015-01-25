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


#ifndef SWIGXML
%include "typemaps.i"
#endif

/// Generic typemap structure

  /// Data structure in the target language holding data
#ifdef SWIGPYTHON
#define GUESTOBJECT PyObject
#elif defined(SWIGMATLAB)
#define GUESTOBJECT mxArray
#else
#define GUESTOBJECT void
#endif

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

%fragment("conv_constref", "header") {
  template<typename T>
  bool conv_constref(GUESTOBJECT *p, T* &ptr, T &m, swig_type_info *type, int (*f)(GUESTOBJECT *p, void *mv, int offs)) {
    if (SWIG_ConvertPtr(p, (void **)&ptr, type, 0) == -1) {
      if (!f(p, &m, 0)) return false;
      ptr = &m;
    }
    return true;
  }
 }

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
  $1 = to_##xName($input, 0);
 }
%enddef

%define %casadi_typecheck_typemap_constref(xName, xPrec, xType...)
%typemap(typecheck, noblock=1, fragment="to"{xName}, precedence=xPrec) const xType& {
  $1 = to_##xName($input, 0);
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

%fragment("try_copy", "header") {
  template<typename FromType>
  struct CopyTraits {
    template<typename ToType>
    static bool try_copy(GUESTOBJECT *p, swig_type_info *type, ToType* m) {
      FromType *mp = 0;
      if (SWIG_ConvertPtr(p, (void **) &mp, type, 0) != -1) {
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
  int to_int(GUESTOBJECT *p, void *mv, int offs=0);
  int to_double(GUESTOBJECT *p, void *mv, int offs=0);
  int to_Dictionary(GUESTOBJECT *p, void *mv, int offs=0);
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

  GUESTOBJECT * from_GenericType(const casadi::GenericType &a);
  GUESTOBJECT * from_Dictionary(const casadi::GenericType::Dictionary &a);

  swig_type_info * type_SX() { return $descriptor(casadi::Matrix<casadi::SXElement> *); }
  swig_type_info * type_DMatrix() { return $descriptor(casadi::Matrix<double> *); }
  swig_type_info * type_IMatrix() { return $descriptor(casadi::Matrix<int> *); }
  swig_type_info * type_MX() { return $descriptor(casadi::MX *); }
}
