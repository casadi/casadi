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


%{


// Returns a new reference
/* PyObject * getReturnType(PyObject* p) { */
/*   if (!p) return 0; */
/*   if (!PyCallable_Check(p)) return 0; */
/*   int ret = PyObject_HasAttrString( p, "func_annotations"); */
/*   if (!ret) return 0; */
/*   PyObject * func_annotations = PyObject_GetAttrString(p,"func_annotations"); */
/*   if (!PyDict_Check(func_annotations)) { */
/*     Py_DECREF(func_annotations); */
/*     return 0; */
/*   } */
/*   PyObject * return_type = PyDict_GetItemString(func_annotations, "return"); // Borrowed */
/*   Py_INCREF(return_type); // Make a new reference */
/*   Py_DECREF(func_annotations); */
/*   if (return_type==0) { */
/*     return 0; */
/*   } */
/*   if (!PyType_Check(return_type) && return_type!=Py_None) { */
/*     Py_DECREF(return_type); */
/*     return 0; */
/*   } */
/*   return return_type; */
/* } */

// Returns a new reference
/* PyObject * getCasadiObject(const std::string &s) { */
/*   PyObject* pPyObjectModuleName = PyString_FromString("casadi"); */
/*   if (!pPyObjectModuleName) { PyErr_Clear(); return 0; } */
/*   PyObject* pObjectModule = PyImport_Import(pPyObjectModuleName); */
/*   Py_DECREF(pPyObjectModuleName); */
/*   if (!pObjectModule) { PyErr_Clear(); return 0; } */
/*   PyObject* pObjectDict = PyModule_GetDict(pObjectModule); // Borrowed */
/*   Py_DECREF(pObjectModule); */
/*   if (!pObjectDict) { PyErr_Clear(); return 0; } */
/*   PyObject* ret = PyDict_GetItemString(pObjectDict,  s.c_str()); // Borrowed */
/*   if (!ret) { PyErr_Clear(); return 0; } */
/*   Py_INCREF(ret); // New reference */
/*   return ret; */
/* } */

#include "casadi/core/functor_internal.hpp"
#include "casadi/core/function/custom_function.hpp"

namespace casadi {
  //using namespace casadi;

  /* class FunctorPythonInternal { */
  /*   public: */
  /*     FunctorPythonInternal(PyObject *p) : p_(p) { Py_INCREF(p_); } */
  /*     ~FunctorPythonInternal() { Py_DECREF(p_); } */
  /*   protected: */
  /*     PyObject *p_; */
  /* }; */
  
  /* class DerivativeGeneratorPythonInternal : public DerivativeGeneratorInternal, FunctorPythonInternal { */
  /*   friend class DerivativeGeneratorPython; */
    
  /* DerivativeGeneratorPythonInternal(PyObject *p) : FunctorPythonInternal(p) {} */
  /*   virtual Function call(Function& fcn, int nfwd, int nadj, void* user_data); */
  /*   virtual DerivativeGeneratorPythonInternal* clone() const { return new DerivativeGeneratorPythonInternal(p_); } */
  /* }; */
   
  /* class CustomEvaluatePythonInternal : public CustomEvaluateInternal, FunctorPythonInternal { */
  /*   friend class CustomEvaluatePython; */
    
  /*   CustomEvaluatePythonInternal(PyObject *p) : FunctorPythonInternal(p) {} */
  /*   virtual void call(CustomFunction& fcn, void* user_data); */
  /*   virtual CustomEvaluatePythonInternal* clone() const { return new CustomEvaluatePythonInternal(p_); } */
  /* }; */
  
  /* class DerivativeGeneratorPython : public DerivativeGenerator { */
  /* public: */
  /*   DerivativeGeneratorPython(PyObject *p) { assignNode(new DerivativeGeneratorPythonInternal(p)); } */
  /* }; */
  
  /* class CallbackPythonInternal : public CallbackInternal, FunctorPythonInternal { */
  /*   friend class CallbackPython; */
    
  /*   CallbackPythonInternal(PyObject *p) : FunctorPythonInternal(p) {} */
  /*   virtual int call(Function& fcn, void* user_data); */
  /*   virtual CallbackPythonInternal* clone() const { return new CallbackPythonInternal(p_); } */
  /* }; */
  
  /* class CustomEvaluatePython : public CustomEvaluate { */
  /*   public: */
  /*     CustomEvaluatePython(PyObject *p) { assignNode(new CustomEvaluatePythonInternal(p)); } */
  /* }; */

  /* class CallbackPython : public Callback { */
  /*   public: */
  /*     CallbackPython(PyObject *p) { assignNode(new CallbackPythonInternal(p)); } */
  /* }; */
  
}
%}

%wrapper %{  

namespace casadi {

  /* Function DerivativeGeneratorPythonInternal::call(Function& fcn, int nfwd, int nadj, void* user_data) { */
  /*   casadi_assert(p_!=0); */
  /*   PyObject * nfwd_py = PyInt_FromLong(nfwd); */
  /*   PyObject * nadj_py = PyInt_FromLong(nadj); */
  /*   PyObject * fcn_py = SWIG_NewPointerObj((new Function(static_cast< const Function& >(fcn))), SWIGTYPE_p_casadi__Function, SWIG_POINTER_OWN |  0 ); */
  /*   if(!fcn_py) { */
  /*     Py_DECREF(nfwd_py); */
  /*     Py_DECREF(nadj_py); */
  /*     throw CasadiException("DerivativeGeneratorPythonInternal: failed to convert Function to python"); */
  /*   } */
    
  /*   PyObject *r = PyObject_CallFunctionObjArgs(p_, fcn_py, nfwd_py, nadj_py, NULL); */
  /*   Py_DECREF(nfwd_py); */
  /*   Py_DECREF(nadj_py); */
  /*   Py_DECREF(fcn_py); */
    
  /*   if (r) { */
  /*     Function ret;   */
  /*     int result = meta< Function >::as(r,ret); */
  /*     if(!result) { Py_DECREF(r); throw CasadiException("DerivativeGeneratorPythonInternal: return type was not Function."); } */
  /*     Py_DECREF(r); */
  /*     return ret; */
  /*   } else { */
  /*     PyErr_Print(); */
  /*     throw CasadiException("DerivativeGeneratorPythonInternal: python method execution raised an Error."); */
  /*   } */
  /* } */

  /* void CustomEvaluatePythonInternal::call(CustomFunction& fcn, void* user_data) { */
  /*   casadi_assert(p_!=0); */
  /*   PyObject * fcn_py = SWIG_NewPointerObj((new CustomFunction(static_cast< const CustomFunction& >(fcn))), SWIGTYPE_p_casadi__CustomFunction, SWIG_POINTER_OWN |  0 ); */
  /*   if(!fcn_py) { */
  /*     throw CasadiException("CustomEvaluatePythonInternal: failed to convert CustomFunction to python"); */
  /*   } */
    
  /*   PyObject *r = PyObject_CallFunctionObjArgs(p_, fcn_py, NULL); */
  /*   Py_DECREF(fcn_py); */
  /*   if (!r) { */
  /*    PyErr_Print(); */
  /*    throw CasadiException("CustomEvaluatePythonInternal: python method execution raised an Error."); */
  /*   } */
    
  /*   Py_DECREF(r); */

  /* } */
  
  /* int CallbackPythonInternal::call(Function& fcn, void* user_data) { */
  /*   casadi_assert(p_!=0); */

  /*   PyObject * fcn_py = SWIG_NewPointerObj((new Function(static_cast< const Function& >(fcn))), SWIGTYPE_p_casadi__CustomFunction, SWIG_POINTER_OWN |  0 ); */
  /*   if(!fcn_py) { */
  /*     throw CasadiException("CallbackPythonInternal: failed to convert CustomFunction to python"); */
  /*   } */
    
  /*   PyObject *r = PyObject_CallFunctionObjArgs(p_, fcn_py, NULL); */

  /*   Py_DECREF(fcn_py); */
  /*   if (!r) { */
  /*    PyErr_Print(); */
  /*    throw CasadiException("CallbackPythonInternal: python method execution raised an Error."); */
  /*   } */
  /*   int ret = 0; */
  /*   if ( meta< int >::couldbe(r)) { */
  /*     /\*int res = *\/ meta< int >::as(r,ret);     */
  /*   } */
    
  /*   Py_DECREF(r); */
  /*   return ret; */

  /* } */
  
}

%}

%inline %{
/// int
template<> char meta< int >::expected_message[] = "Expecting integer";

template <>
int meta< int >::as(GUESTOBJECT * p, int &m) {
  NATIVERETURN(int,m)
  if (meta< casadi::Matrix<int> >::isa(p)) {
    casadi::Matrix<int> *temp = meta< casadi::Matrix<int> >::get_ptr(p);
    if (temp->numel()==1 && temp->size()==1) {
      m = temp->data()[0];
      return true;
    }
    return false;
  } else {
    return false;
  }
}

/// std::vector<double>
template <> char meta< std::vector< double > >::expected_message[] = "Expecting sequence(double)";
/* template <> */
/* int meta< std::vector< double > >::as(GUESTOBJECT * p,std::vector<double > &m) { */
/*   NATIVERETURN(std::vector< double >,m) */
/*   return meta< double >::as_vector(p,m); */
/* } */

/// std::vector<int>
template <> char meta< std::vector< int > >::expected_message[] = "Expecting sequence(integer) or 1D numpy.array of ints"; 
/* template <> */
/* int meta< std::vector< int > >::as(GUESTOBJECT * p,std::vector< int > &m) { */
/*   NATIVERETURN(std::vector< int >,m) */
/*   return meta< int >::as_vector(p,m); */
/* } */

/// double
template<> char meta< double >::expected_message[] = "Expecting double";

template <>
int meta< double >::as(GUESTOBJECT * p, double &m) {
  NATIVERETURN(double,m)
  if (meta< casadi::Matrix<double> >::isa(p)) {
    casadi::Matrix<double> *temp = meta< casadi::Matrix<double> >::get_ptr(p);
    if (temp->numel()==1 && temp->size()==1) {
      m = temp->data()[0];
      return true;
    }
    return false;
  } else if (meta< casadi::Matrix<int> >::isa(p)) {
    casadi::Matrix<int> *temp = meta< casadi::Matrix<int> >::get_ptr(p);
    if (temp->numel()==1 && temp->size()==1) {
      m = temp->data()[0];
      return true;
    }
    return false;
  } else {
    return false;
  }
}

/// std::string
template<> char meta< std::string >::expected_message[] = "Expecting string";

/* template <> */
/* int meta< std::string >::as(PyObject * p, std::string &m) { */
/*   NATIVERETURN(std::string,m) */
/*   if (!PyString_Check(p)) return false; */
/*   m.clear(); */
/*   m.append(PyString_AsString(p)); */
/*   return true; */
/* } */

/* meta_vector(std::string); */

// Forward declarations
/* template<> int meta< casadi::GenericType::Dictionary >::as(PyObject * p,casadi::GenericType::Dictionary &s); */
/* template<> bool meta< casadi::GenericType::Dictionary >::toPython(const casadi::GenericType::Dictionary &a, PyObject *&p); */

/// casadi::DerivativeGenerator
 template<> char meta< casadi::DerivativeGenerator >::expected_message[] = "Expecting sparsity generator";

 /* template <> */
 /*   int meta< casadi::DerivativeGenerator >::as(PyObject * p, casadi::DerivativeGenerator &s) { */
 /*   NATIVERETURN(casadi::DerivativeGenerator, s) */
 /*   PyObject* return_type = getReturnType(p); */
 /*   if (!return_type) return false; */
 /*   PyObject* function = getCasadiObject("Function"); */
 /*   if (!function) { Py_DECREF(return_type); return false; } */
 /*   bool res = PyClass_IsSubclass(return_type,function); */
 /*   Py_DECREF(return_type);Py_DECREF(function); */
 /*   if (res) { */
 /*     s = casadi::DerivativeGeneratorPython(p); */
 /*     return true; */
 /*   } */
 /*   return false; */
 /* } */

/// casadi::CustomEvaluate
template<> char meta< casadi::CustomEvaluate >::expected_message[] = "Expecting CustomFunction wrapper generator";

/* template <> */
/* int meta< casadi::CustomEvaluate >::as(PyObject * p, casadi::CustomEvaluate &s) { */
/*   NATIVERETURN(casadi::CustomEvaluate, s) */
/*   PyObject* return_type = getReturnType(p); */
/*   bool res = (return_type==Py_None) || !return_type; */
/*   if (return_type) Py_DECREF(return_type); */
/*   if (res) { */
/*     s = casadi::CustomEvaluatePython(p); */
/*   } */
/*   return res; */
/* } */

/// casadi::Callback
template<> char meta< casadi::Callback >::expected_message[] = "Expecting Callback";

/* template <> */
/* int meta< casadi::Callback >::as(PyObject * p, casadi::Callback &s) { */
/*   NATIVERETURN(casadi::Callback, s) */
/*   PyObject* return_type = getReturnType(p); */
/*   if (!return_type) return false; */
/*   bool res = PyType_IsSubtype((PyTypeObject *)return_type,&PyInt_Type) || PyType_IsSubtype((PyTypeObject *)return_type,&PyLong_Type); */
/*   Py_DECREF(return_type); */
/*   if (res) { */
/*     s = casadi::CallbackPython(p); */
/*     return true; */
/*   } */
/*   return false; */
/* } */

/// casadi::GenericType
template<> char meta< casadi::GenericType >::expected_message[] = "Expecting any type (None might be an exception)";

/* template <> */
/* int meta< casadi::GenericType >::as(PyObject * p,casadi::GenericType &s) { */
/*   NATIVERETURN(casadi::GenericType, s) */
/*   if (p==Py_None) { */
/*     s=casadi::GenericType(); */
/*   } else if (PyBool_Check(p)) { */
/*     s=casadi::GenericType((bool) PyInt_AsLong(p)); */
/*   } else if (PyInt_Check(p)) { */
/*     s=casadi::GenericType((int) PyInt_AsLong(p)); */
/*   } else if (PyFloat_Check(p)) { */
/*     s=casadi::GenericType(PyFloat_AsDouble(p)); */
/*   } else if (meta< std::string >::couldbe(p)) { */
/*     std::string temp; */
/*     int ret = meta< std::string >::as(p,temp);  */
/*     if (!ret) return false; */
/*     s = casadi::GenericType(temp); */
/*   } else if (meta< std::vector<int> >::couldbe(p)) { */
/*     std::vector<int> temp; */
/*     int ret = meta< std::vector<int> >::as(p,temp);  */
/*     if (!ret) return false; */
/*     s = casadi::GenericType(temp); */
/*   } else if (meta< std::vector<double> >::couldbe(p)) { */
/*     std::vector<double> temp; */
/*     int ret = meta< std::vector<double> >::as(p,temp);  */
/*     if (!ret) return false; */
/*     s = casadi::GenericType(temp); */
/*   } else if (meta< std::vector<std::string> >::couldbe(p)) { */
/*     std::vector<std::string> temp; */
/*     int ret = meta< std::vector<std::string> >::as(p,temp);  */
/*     if (!ret) return false; */
/*     s = casadi::GenericType(temp); */
/*   } else if (PyType_Check(p) && PyObject_HasAttrString(p,"creator")) { */
/*     PyObject *c = PyObject_GetAttrString(p,"creator"); */
/*     if (!c) return false; */
/*     PyObject* gt = getCasadiObject("GenericType"); */
/*     if (!gt) return false; */

/*     PyObject* args = PyTuple_New(1); */
/*     PyTuple_SetItem(args,0,c); */
    
/*     PyObject* g = PyObject_CallObject(gt,args); */
    
/*     Py_DECREF(args); */
/*     Py_DECREF(gt); */
    
/*     if (g) { */
/*       int result = meta< casadi::GenericType >::as(g,s); */
/*       Py_DECREF(g); */
/*       return result; */
/*     } */
/*     if (!g) { PyErr_Clear();  return false;} */
    
/*   } else if (meta< casadi::Function >::couldbe(p)) { */
/*     casadi::Function temp; */
/*     int ret = meta< casadi::Function >::as(p,temp);  */
/*     if (!ret) return false; */
/*     s = casadi::GenericType(temp); */
/*   } else if (meta< casadi::GenericType::Dictionary >::couldbe(p) || meta< casadi::DerivativeGenerator >::couldbe(p) || meta< casadi::Callback >::couldbe(p)) { */
/*     PyObject* gt = getCasadiObject("GenericType"); */
/*     if (!gt) return false; */

/*     PyObject* args = PyTuple_New(1); */
/*     Py_INCREF(p); // Needed because PyTuple_SetItem steals the reference */
/*     PyTuple_SetItem(args,0,p); */
    
/*     PyObject* g = PyObject_CallObject(gt,args); */
    
/*     Py_DECREF(args); */
/*     Py_DECREF(gt); */
    
/*     if (g) { */
/*       int result = meta< casadi::GenericType >::as(g,s); */
/*       Py_DECREF(g); */
/*       return result; */
/*     } */
/*     if (!g) { PyErr_Clear();  return false;} */

    
/*   } else { */
/*     return false; */
/*   } */
/*   return true; */
/* } */

/* template <> */
/* bool meta< casadi::GenericType >::toPython(const casadi::GenericType &a, PyObject *&p) { */
/*   if (a.isBool()) { */
/*     p=PyBool_FromLong(a.toBool()); */
/*   } else if (a.isInt()) { */
/*     p=PyInt_FromLong(a.toInt()); */
/*   } else if (a.isDouble()) { */
/*     p=PyFloat_FromDouble(a.toDouble()); */
/*   } else if (a.isString()) { */
/*     p=PyString_FromString(a.toString().c_str()); */
/*   } else if (a.isIntVector()) { */
/*     p = swig::from(a.toIntVector()); */
/*   } else if (a.isDoubleVector()) { */
/*     p = swig::from( a.toDoubleVector()); */
/*   }  else if (a.isStringVector()) { */
/*     p = swig::from(a.toStringVector()); */
/*   } else if (a.isDictionary()) { */
/*     meta< casadi::GenericType::Dictionary >::toPython(a.toDictionary(),p); */
/*   } else if (a.isNull()) { */
/*     p = Py_None; */
/*   } else { */
/*     return false; */
/*   } */
/*   return true; */
/* } */

/// casadi::GenericType::Dictionary
template<> char meta< casadi::GenericType::Dictionary >::expected_message[] = "Expecting dictionary of GenericTypes";

/* template <> */
/* int meta< casadi::GenericType::Dictionary >::as(PyObject * p,casadi::GenericType::Dictionary &s) { */
/*   NATIVERETURN(casadi::GenericType::Dictionary, s) */
/*   if (!PyDict_Check(p)) */
/*     return false; */
/*   PyObject *key, *value; */
/*   Py_ssize_t pos = 0; */
/*   casadi::GenericType gt; */
/*   while (PyDict_Next(p, &pos, &key, &value)) { */
/*     if (!PyString_Check(key)) */
/*       return false; */
/*     bool ret=meta< casadi::GenericType >::as(value,gt); */
/*     if (!ret) */
/*       return false; */
/*     s[std::string(PyString_AsString(key))] = gt; */
/*   } */

/*   return true; */
/* } */

/* template <> */
/* bool meta< casadi::GenericType::Dictionary >::toPython(const casadi::GenericType::Dictionary &a, PyObject *&p) { */
/*   p = PyDict_New(); */
/*   //casadi::Dictionary::const_iterator end = a.end();  */
/*   casadi::GenericType::Dictionary::const_iterator end = a.end(); */
/*   for (casadi::GenericType::Dictionary::const_iterator it = a.begin(); it != end; ++it) */
/*   { */
/*     PyObject * e; */
/*     bool ret=meta< casadi::GenericType >::toPython(it->second,e); */
/*     if (!ret) { */
/*       Py_DECREF(p); */
/*       return false; */
/*     } */
/*     PyDict_SetItemString(p,(it->first).c_str(),e); */
/*     Py_DECREF(e); */
/*   } */
/*   return true; */
/* } */


// Explicit intialization of these two member functions, so we can use them in meta< casadi::SXElement >
template<> int meta< casadi::SX >::as(GUESTOBJECT *p, casadi::SX &);
//template<> bool meta< casadi::SX >::couldbe(GUESTOBJECT *p);

/// casadi::SX
template<> char meta< casadi::SXElement >::expected_message[] = "Expecting SXElement or number";

template <>
int meta< casadi::SXElement >::as(GUESTOBJECT *p, casadi::SXElement &s) {
  NATIVERETURN(casadi::SXElement, s)
  if (meta< double >::couldbe(p)) {
    double res;
    int result = meta< double >::as(p,res);
    if (!result)
      return false;
    s=casadi::SXElement(res);
  } else {
    return false;
  }
  return true;
}

/// casadi::Matrix<int>
template<> char meta< casadi::Matrix<int> >::expected_message[] = "Expecting IMatrix";

template <>
int meta< casadi::Matrix<int> >::as(GUESTOBJECT * p,casadi::Matrix<int> &m) {
  NATIVERETURN(casadi::Matrix<int>,m)
  if (meta< int >::couldbe(p)) {
    int t;
    int res = meta< int >::as(p,t);
    m = t;
    return res;
  } else {
    //SWIG_Error(SWIG_TypeError, "asDMatrix: unrecognised type. Should have been caught by typemap(typecheck)");
    return false;
  }
  return true;
}

/* meta_vector(casadi::Matrix<int>) */
/* meta_vector(std::vector< casadi::Matrix<int> >) */

/// casadi::Matrix<double>
template<> char meta< casadi::Matrix<double> >::expected_message[] = "Expecting DMatrix";

template <>
int meta< casadi::Matrix<double> >::as(GUESTOBJECT * p,casadi::Matrix<double> &m) {
  NATIVERETURN(casadi::Matrix<double>,m)
  NATIVERETURN(casadi::Matrix<int>,m)
  if (meta< double >::couldbe(p)) {
    double t;
    int res = meta< double >::as(p,t);
    m = t;
    return res;
  } else {
    //SWIG_Error(SWIG_TypeError, "asDMatrix: unrecognised type. Should have been caught by typemap(typecheck)");
    return false;
  }
  return true;
}

/* meta_vector(casadi::Matrix<double>) */
/* meta_vector(std::vector< casadi::Matrix<double> >) */

/// casadi::SX
template<> char meta< casadi::SX >::expected_message[] = "Expecting one of: SX, SXElement, number";

template <>
int meta< casadi::SX >::as(GUESTOBJECT * p,casadi::SX &m) {
  NATIVERETURN(casadi::SX, m)
  NATIVERETURN(casadi::SXElement, m)
  casadi::DMatrix mt;
  if(meta< casadi::Matrix<double> >::as(p,mt)) {
    m = casadi::SX(mt);
  } else {
    //SWIG_Error(SWIG_TypeError, "asSX: unrecognised type. Should have been caught by typemap(typecheck)");
    return false;
  }
  return true;
 }

meta_vector(std::vector<casadi::SXElement>);
meta_vector(casadi::SXElement);
meta_vector(casadi::Matrix< casadi::SXElement >);
meta_vector(std::vector< casadi::Matrix< casadi::SXElement > >);

/// casadi::MX
template<> char meta< casadi::MX >::expected_message[] = "Expecting (MX, numberarray)";

template <>
int meta< casadi::MX >::as(GUESTOBJECT * p,casadi::MX &m) {
  NATIVERETURN(casadi::MX,m)
  casadi::DMatrix mt;
  if(meta< casadi::Matrix<double> >::as(p,mt)) {
    m = casadi::MX(mt);
    return true;
  }
  return false;
}

meta_vector(casadi::MX);
meta_vector(std::vector< casadi::MX >);
%}


