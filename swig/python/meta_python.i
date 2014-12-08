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
  //using namespace casadi;

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
    virtual Function call(Function& fcn, int nfwd, int nadj, void* user_data);
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

%wrapper %{  

namespace casadi {

  Function DerivativeGeneratorPythonInternal::call(Function& fcn, int nfwd, int nadj, void* user_data) {
    casadi_assert(p_!=0);
    PyObject * nfwd_py = PyInt_FromLong(nfwd);
    PyObject * nadj_py = PyInt_FromLong(nadj);
    PyObject * fcn_py = SWIG_NewPointerObj((new Function(static_cast< const Function& >(fcn))), SWIGTYPE_p_casadi__Function, SWIG_POINTER_OWN |  0 );
    if(!fcn_py) {
      Py_DECREF(nfwd_py);
      Py_DECREF(nadj_py);
      throw CasadiException("DerivativeGeneratorPythonInternal: failed to convert Function to python");
    }
    
    PyObject *r = PyObject_CallFunctionObjArgs(p_, fcn_py, nfwd_py, nadj_py, NULL);
    Py_DECREF(nfwd_py);
    Py_DECREF(nadj_py);
    Py_DECREF(fcn_py);
    
    if (r) {
      Function ret;  
      if(!meta< Function >::toCpp(r, &ret, *meta< Function >::name)) {
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

  void CustomEvaluatePythonInternal::call(CustomFunction& fcn, void* user_data) {
    casadi_assert(p_!=0);
    PyObject * fcn_py = SWIG_NewPointerObj((new CustomFunction(static_cast< const CustomFunction& >(fcn))), SWIGTYPE_p_casadi__CustomFunction, SWIG_POINTER_OWN |  0 );
    if(!fcn_py) {
      throw CasadiException("CustomEvaluatePythonInternal: failed to convert CustomFunction to python");
    }
    
    PyObject *r = PyObject_CallFunctionObjArgs(p_, fcn_py, NULL);
    Py_DECREF(fcn_py);
    if (!r) {
     PyErr_Print();
     throw CasadiException("CustomEvaluatePythonInternal: python method execution raised an Error.");
    }
    
    Py_DECREF(r);

  }
  
  int CallbackPythonInternal::call(Function& fcn, void* user_data) {
    casadi_assert(p_!=0);

    PyObject * fcn_py = SWIG_NewPointerObj((new Function(static_cast< const Function& >(fcn))), SWIGTYPE_p_casadi__CustomFunction, SWIG_POINTER_OWN |  0 );
    if(!fcn_py) {
      throw CasadiException("CallbackPythonInternal: failed to convert CustomFunction to python");
    }
    
    PyObject *r = PyObject_CallFunctionObjArgs(p_, fcn_py, NULL);

    Py_DECREF(fcn_py);
    if (!r) {
     PyErr_Print();
     throw CasadiException("CallbackPythonInternal: python method execution raised an Error.");
    }
    int ret = 0;
    if ( meta< int >::toCpp(r, 0, *meta< int >::name)) {
      /*int res = */ meta< int >::toCpp(r, &ret, *meta< int >::name);    
    }
    
    Py_DECREF(r);
    return ret;

  }
  
}

%}

%inline %{
/// int
template <>
  int meta< int >::toCpp(PyObject * p, int *m, swig_type_info *type) {
  int *mp = 0;
  if (SWIG_ConvertPtr(p, (void **) &mp, type, 0) != -1) {
    if (m) *m=*mp;
    return true;
  } else if (PyInt_Check(p) || PyLong_Check(p) || PyBool_Check(p)) {
    PyObject *r = PyNumber_Long(p);
    if (!r) return false;
    if (m) *m = PyLong_AsLong(r);
    Py_DECREF(r);
    return true;
  } else if (PyObject_HasAttrString(p,"dtype")) {
    PyObject *r = PyObject_GetAttrString(p,"dtype");
    if (!PyObject_HasAttrString(r,"kind")) { Py_DECREF(r); return false;}
    PyObject *k = PyObject_GetAttrString(r,"kind");
    if (!PyObject_HasAttrString(p,"__int__")) { Py_DECREF(k);Py_DECREF(r); return false;}
    char name[] = "__int__";
    PyObject *mm = PyObject_CallMethod(p, name,0);
    if (!mm) { PyErr_Clear(); Py_DECREF(k); Py_DECREF(r); return false; }
    char *kk = PyString_AsString(k);
    bool result =  kk[0]=='i';
    Py_DECREF(k);
    Py_DECREF(r);
    if (result) {
      if (m) *m = PyLong_AsLong(mm);
    }
    Py_DECREF(mm);
    return result;
  } else if (is_a(p, *meta< casadi::Matrix<int> >::name)) {
    casadi::Matrix<int> *temp;
    SWIG_ConvertPtr(p, (void **) &temp, *meta< casadi::Matrix<int> >::name, 0 );
    if (temp->numel()==1 && temp->size()==1) {
      if (m) *m = temp->data()[0];
      return true;
    }
    return false;
  } else {
    return false;
  }
}

// Forward declarations
GUESTOBJECT * fromCpp(const casadi::GenericType::Dictionary &a);

/// casadi::DerivativeGenerator
template <>
int meta< casadi::DerivativeGenerator >::toCpp(PyObject * p, casadi::DerivativeGenerator *m, swig_type_info *type) {
  casadi::DerivativeGenerator *mp = 0;
  if (SWIG_ConvertPtr(p, (void **) &mp, type, 0) != -1) {
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

GUESTOBJECT * fromCpp(const casadi::GenericType &a) {
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
  } else if (a.isDoubleVector()) {
    p = swig::from( a.toDoubleVector());
  }  else if (a.isStringVector()) {
    p = swig::from(a.toStringVector());
  } else if (a.isDictionary()) {
    p = fromCpp(a.toDictionary());
  } else if (a.isNull()) {
    p = Py_None;
  }
  return p;
}

PyObject * fromCpp(const casadi::GenericType::Dictionary &a) {
  PyObject *p = PyDict_New();
  casadi::GenericType::Dictionary::const_iterator end = a.end();
  for (casadi::GenericType::Dictionary::const_iterator it = a.begin(); it != end; ++it) {
    PyObject * e = fromCpp(it->second);
    if (!e) {
      Py_DECREF(p);
      return 0;
    }
    PyDict_SetItemString(p,(it->first).c_str(),e);
    Py_DECREF(e);
  }
  return p;
}

%}
