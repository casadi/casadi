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
