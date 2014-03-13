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

%{
#include "symbolic/matrix/sparsity.hpp"
#include "symbolic/matrix/matrix.hpp"
#include <sstream>
#include "symbolic/casadi_exception.hpp"

// to allow for typechecking
#include "symbolic/sx/sx_element.hpp"

// to typecheck for MX
#include "symbolic/mx/mx.hpp"
%}

#ifndef SWIGXML
%include "typemaps.i"
#endif

/// Generic typemap structure
%inline %{

#ifdef SWIGPYTHON
#define GUESTOBJECT PyObject * p
#endif // SWIGPYTHON

#ifdef  SWIGOCTAVE
#define GUESTOBJECT const octave_value& p
#endif // SWIGOCTAVE

#ifdef  SWIGOCTAVE
/** Check if Guest object is of a particular SWIG type */
bool istype(GUESTOBJECT, swig_type_info *type) {
  if (p.is_cell() && p.rows() == 1 && p.columns() == 1)
    return istype(p.cell_value()(0),type);
  if (!p.is_defined())
    return false;
  if (p.type_id() != octave_swig_ref::static_type_id())
    return false;
  octave_swig_ref *osr = static_cast < octave_swig_ref *>(p.internal_rep());
  octave_swig_type *ost = osr->get_ptr();
  void *vptr = ost->cast(type,0, 0);
  if (!vptr)
    return false;
  return true;
}
#endif // SWIGOCTAVE

#ifdef  SWIGPYTHON
/** Check if Guest object is of a particular SWIG type */
bool istype(GUESTOBJECT, swig_type_info *type) {
	void *dummy = 0 ; 
  return SWIG_IsOK(SWIG_ConvertPtr(p, &dummy, type, 0)); 
}
#endif // SWIGPYTHON

template<typename DataType>
class meta {
  public:
    /// Check if Python object is of type DataType
    static bool isa(GUESTOBJECT) {
      #ifdef SWIGPYTHON
      if (p == Py_None) return false;
      #endif // SWIGPYTHON
      #ifdef SWIGOCTAVE
      if (p.is_null_value()) return false;
      #endif // SWIGOCTAVE
      return istype(p,*meta<DataType>::name);
    };
    /// Convert Python object to pointer of type DataType
    static bool get_ptr(GUESTOBJECT,DataType*& m) {
      void *pd = 0 ;
      int res = SWIG_ConvertPtr(p, &pd,*meta<DataType>::name, 0 );
      if (!SWIG_IsOK(res)) {
        return false;
      }
      m = reinterpret_cast< DataType*  >(pd);
      return true;
    };
    /// Convert Guest object to type DataType
    /// This function must work when isa(GUESTOBJECT) too
    static int as(GUESTOBJECT,DataType& m) {
        DataType *t = (DataType *)(0);
        int res = swig::asptr(p, &t);
        bool succes = SWIG_CheckState(res) && t;
        if (succes) m=*t;
        if (succes && SWIG_IsNewObj(res)) delete t;
        return succes;
    }
    /// Check if Guest object could ultimately be converted to type DataType
    /// may return true when isa(GUESTOBJECT), but this is not required.
    static bool couldbe(GUESTOBJECT) { 
        int res = swig::asptr(p, (DataType**)(0));
        return SWIG_CheckState(res);
    }
    static swig_type_info** name;
    static char expected_message[];
    
    // Vector specific stuff
    
    #ifdef SWIGPYTHON
    static bool couldbe_sequence(PyObject * p) {
      if(PySequence_Check(p) && !PyString_Check(p) && !meta< CasADi::Matrix<CasADi::SXElement> >::isa(p) && !meta< CasADi::MX >::isa(p) && !meta< CasADi::Matrix<int> >::isa(p) && !meta< CasADi::Matrix<double> >::isa(p) &&!PyObject_HasAttrString(p,"__DMatrix__") && !PyObject_HasAttrString(p,"__SX__") && !PyObject_HasAttrString(p,"__MX__")) {
        PyObject *it = PyObject_GetIter(p);
        if (!it) return false;
        PyObject *pe;
        while ((pe = PyIter_Next(it))) {                                // Iterate over the sequence inside the sequence
          if (!meta<DataType>::couldbe(pe)) {
            Py_DECREF(pe);Py_DECREF(it);return false;
          }
          Py_DECREF(pe);
        }
        Py_DECREF(it);
        if (PyErr_Occurred()) { PyErr_Clear(); return false; }
        return true;
      } else {
        return false;
      }
    }
    #endif // SWIGPYTHON
    
   #ifdef SWIGOCTAVE
    static bool couldbe_sequence(const octave_value& p) {
      if (!p.is_cell()) return false;
      int nrow = p.rows();
      int ncol = p.columns();
      if (nrow!=1 && ncol!=1) return false;
      for(int i=0; i<p.length(); ++i){
        if (!meta<DataType>::couldbe(p.cell_value()(i))) return false;
      }
      return true;
    }
    #endif // SWIGOCTAVE
    
    // Assumes that p is a PYTHON sequence
    #ifdef SWIGPYTHON
    static int as_vector(PyObject * p, std::vector<DataType> &m) {
      PyObject *it = PyObject_GetIter(p);
      PyObject *pe;
      m.resize(PySequence_Size(p));
      int i=0;
      while ((pe = PyIter_Next(it))) {                                // Iterate over the sequence inside the sequence
        bool result=meta<DataType>::as(pe,m[i++]);
        if (!result) {
          Py_DECREF(pe);Py_DECREF(it);
          return false;
        }
        Py_DECREF(pe);
      }
      Py_DECREF(it);
      return true;
    }
    #endif // SWIGPYTHON
    
    // Assumes that p is an octave cell
    #ifdef SWIGOCTAVE
    static int as_vector(const octave_value& p , std::vector<DataType> &m) {
      if (!p.is_cell()) return false;
      int nrow = p.rows();
      int ncol = p.columns();
      if (nrow!=1 && ncol!=1) return false;
      m.resize(p.length());
      
      Cell c = p.cell_value();
      for(int i=0; i<p.length(); ++i){
        // Get the octave object
        const octave_value& obj_i = c(i);
        
        if (!(obj_i.is_real_matrix() && obj_i.is_empty())) {
          bool ret = meta<DataType>::as(obj_i,m[i]);
          if(!ret) return false;
        }
      }
      return true;
    }
    #endif // SWIGOCTAVE
    
    
    #ifdef SWIGPYTHON
    // Would love to make this const T&, but looks like not allowed
    static bool toPython(const DataType &, PyObject *&p);
    #endif //SWIGPYTHON
};

%}

%inline %{
#define NATIVERETURN(Type, m) if (meta<Type>::isa(p)) { Type *mp; int result = meta<Type>::get_ptr(p,mp); if (!result) return false; m=*mp; return true;}
%}


%define %my_generic_const_typemap(Precedence,Type...) 
%typemap(in) const Type & (Type m) {
  if (meta< Type >::isa($input)) { // Type object get passed on as-is, and fast.
    int result = meta< Type >::get_ptr($input,$1);
    if (!result)
      SWIG_exception_fail(SWIG_TypeError,"Type cast failed");
  } else {
    bool result=meta< Type >::as($input,m);
    if (!result)
      SWIG_exception_fail(SWIG_TypeError,meta< Type >::expected_message);
    $1 = &m;
  }
}

%typemap(typecheck,precedence=Precedence) const Type & { $1 = meta< Type >::isa($input) || meta< Type >::couldbe($input); }
%typemap(freearg) const Type  & {}

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
%typemap(out) const Type & {
   $result = SWIG_NewPointerObj((new Type(*$1)), *meta< Type >::name, SWIG_POINTER_OWN |  0 );
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
%typemap(out) Type & {
   $result = SWIG_NewPointerObj($1, *meta< Type >::name, 0 |  0 );
   PySetParent($result, obj0);
}
%typemap(out) const Type & {
   $result = swig::from(static_cast< Type * >($1));
   PySetParent($result, obj0);
}
%extend Type {
%pythoncode%{
    def __del__(self):
      if not(_casadi_main_module is None):
         _casadi_main_module.PyDECREFParent(self)

%}

}
%enddef

%inline %{
/// std::vector< Type >
#define meta_vector(Type) \
template<> char meta< std::vector< Type > >::expected_message[] = "Expecting sequence(Type)"; \
 \
template <> \
int meta< std::vector< Type > >::as(GUESTOBJECT,std::vector< Type > &m) { \
  NATIVERETURN(std::vector< Type >,m) \
  return meta< Type >::as_vector(p,m); \
} \
 \
template <> \
bool meta< std::vector< Type > >::couldbe(GUESTOBJECT) { \
  return meta< std::vector< Type > >::isa(p) ||  meta< Type >::couldbe_sequence(p); \
}
%}


/// std::pair< TypeA, TypeB >
#ifdef SWIGPYTHON
%define %meta_pair(TypeA,TypeB) 
%inline %{
template <>
int meta< std::pair< TypeA, TypeB > >::as(PyObject * p,std::pair< TypeA, TypeB > &m) {
  if(!PySequence_Check(p)) return false;
  if(PySequence_Size(p)!=2) return false;
  PyObject * first =  PySequence_GetItem(p,0);
  PyObject * second = PySequence_GetItem(p,1);
  bool result = meta< TypeA  >::as(p,m.first) && meta< TypeB  >::as(p,m.second);
  
  Py_DECREF(first);Py_DECREF(second);
  return result;   
}

template <>
bool meta< std::pair< TypeA, TypeB > >::couldbe(PyObject * p) {
  if (meta< std::pair< TypeA, TypeB > >::isa(p)) return true; 
  if(!PySequence_Check(p)) return false;
  if(PySequence_Size(p)!=2) return false;
  PyObject * first =  PySequence_GetItem(p,0);
  PyObject * second = PySequence_GetItem(p,1);
  
  bool success = (meta< TypeA >::isa(first) || meta< TypeA >::couldbe(first)) &&
                 (meta< TypeB >::isa(second) || meta< TypeB >::couldbe(second));
  Py_DECREF(first);Py_DECREF(second);              
  return success;
}

template <>
bool meta < std::pair< TypeA, TypeB > >::toPython(const std::pair< TypeA, TypeB > &m, PyObject *&p) {
  p = PyTuple_New(2);
  PyObject *first = 0;
  first  = SWIG_NewPointerObj((new TypeA(static_cast< const TypeA& >(m.first ))), *meta< TypeA >::name , SWIG_POINTER_OWN |  0 );
  PyObject *second = 0;
  second = SWIG_NewPointerObj((new TypeB(static_cast< const TypeB& >(m.second))), *meta< TypeB >::name , SWIG_POINTER_OWN |  0 );
  
  if (first==0 || second==0) return false;
  
  PyTuple_SetItem(p, 0, first);
  PyTuple_SetItem(p, 1, second);
  
  return true;
}
%}
%enddef
#endif //SWIGPYTHON

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

bool PyIsSequence(PyObject* p) {
  return PySequence_Check(p) && !meta< CasADi::Matrix<CasADi::SXElement> >::isa(p) && !meta< CasADi::MX >::isa(p);
}

%}
#endif // SWIGPYTHON

%{
#define SWIG_Error_return(code, msg)  { std::cerr << "Error occured in CasADi SWIG interface code:" << std::endl << "  "<< msg << std::endl;SWIG_Error(code, msg); return 0; }
%}
