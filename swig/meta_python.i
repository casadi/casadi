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

%pythoncode %{
def attach_return_type(f,t):
  if not(hasattr(f,'func_annotations')):
    f.func_annotations = {}
  if not(isinstance(getattr(f,'func_annotations'),dict)):
    raise Exception("Cannot annotate this python Method to be a sparsitygenerator. Method has func_annotations attribute with unknown type.")
  f.func_annotations["return"] = t
  return f

def pyderivativegenerator(f):
  return attach_return_type(f,FX)

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


%{


// Returns a new reference
PyObject * getReturnType(PyObject* p) {
  if (!p) return 0;
  if (!PyCallable_Check(p)) return 0;
  int ret = PyObject_HasAttrString( p, "func_annotations");
  if (!ret) return 0;
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
  if (!pPyObjectModuleName) { PyErr_Clear(); return 0; }
  PyObject* pObjectModule = PyImport_Import(pPyObjectModuleName);
  Py_DECREF(pPyObjectModuleName);
  if (!pObjectModule) { PyErr_Clear(); return 0; }
  PyObject* pObjectDict = PyModule_GetDict(pObjectModule); // Borrowed
  Py_DECREF(pObjectModule);
  if (!pObjectDict) { PyErr_Clear(); return 0; }
  PyObject* ret = PyDict_GetItemString(pObjectDict,  s.c_str()); // Borrowed
  if (!ret) { PyErr_Clear(); return 0; }
  Py_INCREF(ret); // New reference
  return ret;
}

#include "symbolic/functor_internal.hpp"
#include "symbolic/fx/custom_function.hpp"

namespace CasADi {
  //using namespace CasADi;

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
    virtual FX call(FX& fcn, int nfwd, int nadj, void* user_data);
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
    virtual int call(FX& fcn, void* user_data);
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

namespace CasADi {

  FX DerivativeGeneratorPythonInternal::call(FX& fcn, int nfwd, int nadj, void* user_data) {
    casadi_assert(p_!=0);
    PyObject * nfwd_py = PyInt_FromLong(nfwd);
    PyObject * nadj_py = PyInt_FromLong(nadj);
    PyObject * fcn_py = SWIG_NewPointerObj((new FX(static_cast< const FX& >(fcn))), SWIGTYPE_p_CasADi__FX, SWIG_POINTER_OWN |  0 );
    if(!fcn_py) {
      Py_DECREF(nfwd_py);
      Py_DECREF(nadj_py);
      throw CasadiException("DerivativeGeneratorPythonInternal: failed to convert FX to python");
    }
    
    PyObject *r = PyObject_CallFunctionObjArgs(p_, fcn_py, nfwd_py, nadj_py, NULL);
    Py_DECREF(nfwd_py);
    Py_DECREF(nadj_py);
    Py_DECREF(fcn_py);
    
    if (r) {
      FX ret;  
      int result = meta< FX >::as(r,ret);
      if(!result) { Py_DECREF(r); throw CasadiException("DerivativeGeneratorPythonInternal: return type was not FX."); }
      Py_DECREF(r);
      return ret;
    } else {
      PyErr_Print();
      throw CasadiException("DerivativeGeneratorPythonInternal: python method execution raised an Error.");
    }
  }

  void CustomEvaluatePythonInternal::call(CustomFunction& fcn, void* user_data) {
    casadi_assert(p_!=0);
    PyObject * fcn_py = SWIG_NewPointerObj((new CustomFunction(static_cast< const CustomFunction& >(fcn))), SWIGTYPE_p_CasADi__CustomFunction, SWIG_POINTER_OWN |  0 );
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
  
  int CallbackPythonInternal::call(FX& fcn, void* user_data) {
    casadi_assert(p_!=0);

    PyObject * fcn_py = SWIG_NewPointerObj((new FX(static_cast< const FX& >(fcn))), SWIGTYPE_p_CasADi__CustomFunction, SWIG_POINTER_OWN |  0 );
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
    if ( meta< int >::couldbe(r)) {
      int res = meta< int >::as(r,ret);    
    }
    
    Py_DECREF(r);
    return ret;

  }
  
}

%}

%inline %{
/// int
template<> char meta< int >::expected_message[] = "Expecting integer";

template <>
int meta< int >::as(PyObject * p, int &m) {
  NATIVERETURN(int,m)
  if (PyInt_Check(p) || PyLong_Check(p) || PyBool_Check(p)) {
    PyObject *r = PyNumber_Long(p);
    if (!r) return false;
    m = PyLong_AsLong(r);
    Py_DECREF(r);
    return true;
  } else if (PyObject_HasAttrString(p,"__int__")) {
    char name[] = "__int__";
    PyObject *r = PyObject_CallMethod(p, name,0);
    if (!r) { Py_DECREF(r); return false;}
    m = PyLong_AsLong(r);
    Py_DECREF(r);
    return true;
  } else {
    return false;
  }
}


template <> bool meta< int >::couldbe(PyObject * p) {
 if (PyObject_HasAttrString(p,"dtype")) {
   PyObject *r = PyObject_GetAttrString(p,"dtype");
   if (!PyObject_HasAttrString(r,"kind")) { Py_DECREF(r); return false;}
   PyObject *k = PyObject_GetAttrString(r,"kind");
   if (!PyObject_HasAttrString(p,"__int__")) { Py_DECREF(k);Py_DECREF(r); return false;}
   char name[] = "__int__";
   PyObject *m = PyObject_CallMethod(p, name,0);
   if (!m) { PyErr_Clear(); Py_DECREF(k); Py_DECREF(r); return false; }
   char *kk = PyString_AsString(k);
   bool result =  kk[0]=='i';
   Py_DECREF(k); Py_DECREF(r); Py_DECREF(m);
   return result;
   
 }
 if (meta< CasADi::Matrix<int> >::isa(p)) {
   CasADi::Matrix<int> *temp;
   meta< CasADi::Matrix<int> >::get_ptr(p,temp);
   if (temp->numel()==1 && temp->size()==1) return true;
 }
 return PyInt_Check(p) || PyLong_Check(p) || PyBool_Check(p);
}


/// std::vector<double>
template<> char meta< std::vector< double > >::expected_message[] = "Expecting sequence(double)"; 
template <>
int meta< std::vector< double > >::as(PyObject * p,std::vector<double > &m) {
  NATIVERETURN(std::vector< double >,m)
  if (is_array(p)) {
    if (!(array_numdims(p)==1 && array_type(p)!=NPY_OBJECT))
      SWIG_Error_return(SWIG_TypeError, "std::vector<int>: array must be 1D and of a numeric type");
    int size = array_size(p,0);
    if (!array_is_native(p))
      SWIG_Error_return(SWIG_TypeError, "std::vector<double>: array byte order should be native.");
    // Make sure we have a contigous array with double datatype
    int array_is_new_object;
    PyArrayObject* array = obj_to_array_contiguous_allow_conversion(p,NPY_DOUBLE,&array_is_new_object);
    if (!array) { PyErr_Print() ; SWIG_Error_return(SWIG_TypeError, "asMatrixDouble: no luck converting numpy array to double"); }
    double* d=(double*) array_data(array);
    
    m.assign( d, d+size );
    
                  
    // Free memory
    if (array_is_new_object)
      Py_DECREF(array); 
    return true;
  }
  return meta< double >::as_vector(p,m);
}

template <> 
bool meta< std::vector< double > >::couldbe(PyObject * p) {
  return meta< std::vector< double > >::isa(p) || (is_array(p) && array_numdims(p)==1 && array_type(p)!=NPY_OBJECT) || (meta< double >::couldbe_sequence(p) && !is_array(p) && !(meta< CasADi::Matrix<CasADi::SX> >::isa(p) || meta< CasADi::Matrix<double> >::isa(p)) );
}

/// std::vector<int>
template<> char meta< std::vector< int > >::expected_message[] = "Expecting sequence(integer) or 1D numpy.array of ints"; 
template <>
int meta< std::vector< int > >::as(PyObject * p,std::vector< int > &m) {
  NATIVERETURN(std::vector< int >,m)
  if (is_array(p)) {
    if (!(array_numdims(p)==1 && array_type(p)!=NPY_OBJECT))
      SWIG_Error_return(SWIG_TypeError, "std::vector<int>: array must be 1D and of a numeric type");
    int size = array_size(p,0);
    if (!array_is_native(p))
      SWIG_Error_return(SWIG_TypeError, "std::vector<int>: array byte order should be native.");
      
    // Make sure we have a contigous array with int datatype
    int array_is_new_object;
    PyArrayObject* array = obj_to_array_contiguous_allow_conversion(p,NPY_INT,&array_is_new_object);
    if (!array) { // Trying LONG
      array = obj_to_array_contiguous_allow_conversion(p,NPY_LONG,&array_is_new_object);
      if (!array) { PyErr_Print() ; SWIG_Error_return(SWIG_TypeError, "std::vector<int>: no luck converting numpy array to int. Better don't use unsigned datatypes."); }
      long* temp=(long*) array_data(array);
      m.resize(size);
      for (int k=0;k<size;k++) m[k]=temp[k];
      return true;
    }
    int *d=(int*) array_data(array);

    m.assign( d, d+size );

                  
    // Free memory
    if (array_is_new_object)
      Py_DECREF(array); 
    return true;
  }
  return meta< int >::as_vector(p,m);
}

template <> 
bool meta< std::vector< int > >::couldbe(PyObject * p) {
  return meta< std::vector< int > >::isa(p) || (is_array(p) && array_numdims(p)==1) || (meta< int >::couldbe_sequence(p) && !is_array(p) && !(meta< CasADi::Matrix<CasADi::SX> >::isa(p) || meta< CasADi::Matrix<double> >::isa(p)));
}


/// double
template<> char meta< double >::expected_message[] = "Expecting double";

template <>
int meta< double >::as(PyObject * p, double &m) {
  NATIVERETURN(double,m)
  if (PyInt_Check(p) || PyLong_Check(p) || PyBool_Check(p) || PyFloat_Check(p)) {
    PyObject *r = PyNumber_Float(p);
    if (!r) return false;
    m = PyFloat_AsDouble(r);
    Py_DECREF(r);
    return true;
  } else if (PyObject_HasAttrString(p,"__float__")) {
    char name[] = "__float__";
    PyObject *r = PyObject_CallMethod(p, name,0);
    if (!r) { PyErr_Clear();return false;}
    m = PyFloat_AsDouble(r);
    Py_DECREF(r);
    return true;
  } else {
    return false;
  }
}
   
template <> bool meta< double >::couldbe(PyObject * p) {
 if (PyInt_Check(p) || PyLong_Check(p) || PyBool_Check(p) || PyFloat_Check(p)) return true;
 
 if (PyObject_HasAttrString(p,"dtype")) {
   PyObject *r = PyObject_GetAttrString(p,"dtype");
   if (!PyObject_HasAttrString(r,"kind")) { Py_DECREF(r); return false;}
   PyObject *k = PyObject_GetAttrString(r,"kind");
   if (!PyObject_HasAttrString(p,"__float__")) { Py_DECREF(k);Py_DECREF(r); return false;}
   char name[] = "__float__";
   PyObject *m = PyObject_CallMethod(p, name,NULL);
   if (!m) {   PyErr_Clear(); Py_DECREF(k); Py_DECREF(r); return false; }
   char *kk = PyString_AsString(k);
   bool result = kk[0]=='f' || kk[0]=='i';
   Py_DECREF(k); Py_DECREF(r); Py_DECREF(m);
   return result;
 } 
 if (meta< CasADi::Matrix<double> >::isa(p)) {
   CasADi::Matrix<double> *temp;
   meta< CasADi::Matrix<double> >::get_ptr(p,temp);
   if (temp->numel()==1 && temp->size()==1) return true;
 }

  return false;
}

/// std::string
template<> char meta< std::string >::expected_message[] = "Expecting string";

template <>
int meta< std::string >::as(PyObject * p, std::string &m) {
  m.clear();
  if (!PyString_Check(p)) return false;
  m.append(PyString_AsString(p));
  return true;
}
   
template <> bool meta< std::string >::couldbe(PyObject * p) {
    return PyString_Check(p);
}

meta_vector(std::string);

// Forward declarations
template<> int meta< CasADi::GenericType::Dictionary >::as(PyObject * p,CasADi::GenericType::Dictionary &s);
template<> bool meta< CasADi::GenericType::Dictionary >::couldbe(PyObject * p);
template<> bool meta< CasADi::GenericType::Dictionary >::toPython(const CasADi::GenericType::Dictionary &a, PyObject *&p);

/// CasADi::DerivativeGenerator
 template<> char meta< CasADi::DerivativeGenerator >::expected_message[] = "Expecting sparsity generator";

 template <>
   int meta< CasADi::DerivativeGenerator >::as(PyObject * p, CasADi::DerivativeGenerator &s) {
   NATIVERETURN(CasADi::DerivativeGenerator, s)
     s = CasADi::DerivativeGeneratorPython(p);
   return true;
 }

 template<> bool meta< CasADi::DerivativeGenerator >::couldbe(PyObject * p) {
   PyObject* return_type = getReturnType(p);
   if (!return_type) return false;
   PyObject* fx = getCasadiObject("FX");
   if (!fx) { Py_DECREF(return_type); return false; }
   bool res = PyClass_IsSubclass(return_type,fx);
   Py_DECREF(return_type);Py_DECREF(fx);
   return res;
 }


/// CasADi::CustomEvaluate
template<> char meta< CasADi::CustomEvaluate >::expected_message[] = "Expecting CustomFunction wrapper generator";

template <>
int meta< CasADi::CustomEvaluate >::as(PyObject * p, CasADi::CustomEvaluate &s) {
  NATIVERETURN(CasADi::CustomEvaluate, s)
  s = CasADi::CustomEvaluatePython(p);
  return true;
}

template<> bool meta< CasADi::CustomEvaluate >::couldbe(PyObject * p) {
  PyObject* return_type = getReturnType(p);
  if (!return_type) return true;
  bool res = return_type==Py_None;
  Py_DECREF(return_type);
  return res;
}

/// CasADi::Callback
template<> char meta< CasADi::Callback >::expected_message[] = "Expecting Callback";

template <>
int meta< CasADi::Callback >::as(PyObject * p, CasADi::Callback &s) {
  NATIVERETURN(CasADi::Callback, s)
  s = CasADi::CallbackPython(p);
  return true;
}

template<> bool meta< CasADi::Callback >::couldbe(PyObject * p) {
  PyObject* return_type = getReturnType(p);
  if (!return_type) return false;
  bool res = PyType_IsSubtype((PyTypeObject *)return_type,&PyInt_Type) || PyType_IsSubtype((PyTypeObject *)return_type,&PyLong_Type);
  Py_DECREF(return_type);
  return res;
}


/// CasADi::GenericType
template<> char meta< CasADi::GenericType >::expected_message[] = "Expecting any type (None might be an exception)";

template <>
int meta< CasADi::GenericType >::as(PyObject * p,CasADi::GenericType &s) {
  NATIVERETURN(CasADi::GenericType, s)
  if (p==Py_None) {
    s=CasADi::GenericType();
  } else if (PyBool_Check(p)) {
    s=CasADi::GenericType((bool) PyInt_AsLong(p));
  } else if (PyInt_Check(p)) {
    s=CasADi::GenericType((int) PyInt_AsLong(p));
  } else if (PyFloat_Check(p)) {
    s=CasADi::GenericType(PyFloat_AsDouble(p));
  } else if (meta< std::string >::couldbe(p)) {
    std::string temp;
    int ret = meta< std::string >::as(p,temp); 
    if (!ret) return false;
    s = CasADi::GenericType(temp);
  } else if (meta< std::vector<int> >::couldbe(p)) {
    std::vector<int> temp;
    int ret = meta< std::vector<int> >::as(p,temp); 
    if (!ret) return false;
    s = CasADi::GenericType(temp);
  } else if (meta< std::vector<double> >::couldbe(p)) {
    std::vector<double> temp;
    int ret = meta< std::vector<double> >::as(p,temp); 
    if (!ret) return false;
    s = CasADi::GenericType(temp);
  } else if (meta< std::vector<std::string> >::couldbe(p)) {
    std::vector<std::string> temp;
    int ret = meta< std::vector<std::string> >::as(p,temp); 
    if (!ret) return false;
    s = CasADi::GenericType(temp);
  } else if (PyType_Check(p) && PyObject_HasAttrString(p,"creator")) {
    PyObject *c = PyObject_GetAttrString(p,"creator");
    int result = meta< CasADi::GenericType >::as(c,s);
    Py_DECREF(c);
    return result;
  } else if (meta< CasADi::FX >::couldbe(p)) {
    CasADi::FX temp;
    int ret = meta< CasADi::FX >::as(p,temp); 
    if (!ret) return false;
    s = CasADi::GenericType(temp);
  } else {
    PyObject* gt = getCasadiObject("GenericType");
    if (!gt) return false;

    PyObject* args = PyTuple_New(1);
    Py_INCREF(p); // Needed because PyTuple_SetItem steals the reference
    PyTuple_SetItem(args,0,p);
    
    PyObject* g = PyObject_CallObject(gt,args);
    
    Py_DECREF(args);
    Py_DECREF(gt);
    
    if (g) {
      int result = meta< CasADi::GenericType >::as(g,s);
      Py_DECREF(g);
      return result;
    }
    if (!g) { PyErr_Clear();  return false;}

    
  }
  return true;
}

template <>
bool meta< CasADi::GenericType >::couldbe(PyObject * p) {
  return meta< CasADi::GenericType >::isa(p) || PyBool_Check(p) ||  PyInt_Check(p) || PyFloat_Check(p) || PyString_Check(p) || meta< std::vector<double> >::couldbe(p) || meta< CasADi::FX >::couldbe(p) || meta< std::vector<std::string> >::couldbe(p) || PyType_Check(p) || meta< CasADi::GenericType::Dictionary >::couldbe(p) || meta< CasADi::DerivativeGenerator >::couldbe(p) || meta< CasADi::Callback >::couldbe(p);

  
  }

template <>
bool meta< CasADi::GenericType >::toPython(const CasADi::GenericType &a, PyObject *&p) {
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
    meta< CasADi::GenericType::Dictionary >::toPython(a.toDictionary(),p);
  } else if (a.isNull()) {
    p = Py_None;
  } else {
    return false;
  }
  return true;
}

/// CasADi::GenericType::Dictionary
template<> char meta< CasADi::GenericType::Dictionary >::expected_message[] = "Expecting dictionary of GenericTypes";

template <>
int meta< CasADi::GenericType::Dictionary >::as(PyObject * p,CasADi::GenericType::Dictionary &s) {
  NATIVERETURN(CasADi::GenericType::Dictionary, s)
  if (!PyDict_Check(p))
    return false;
  PyObject *key, *value;
  Py_ssize_t pos = 0;
  CasADi::GenericType gt;
  while (PyDict_Next(p, &pos, &key, &value)) {
    if (!PyString_Check(key))
      return false;
    bool ret=meta< CasADi::GenericType >::as(value,gt);
    if (!ret)
      return false;
    s[std::string(PyString_AsString(key))] = gt;
  }

  return true;
}

template <>
bool meta< CasADi::GenericType::Dictionary >::couldbe(PyObject * p) {
  return PyDict_Check(p);
}

template <>
bool meta< CasADi::GenericType::Dictionary >::toPython(const CasADi::GenericType::Dictionary &a, PyObject *&p) {
  p = PyDict_New();
  //CasADi::Dictionary::const_iterator end = a.end(); 
  CasADi::GenericType::Dictionary::const_iterator end = a.end();
  for (CasADi::GenericType::Dictionary::const_iterator it = a.begin(); it != end; ++it)
  {
    PyObject * e;
    bool ret=meta< CasADi::GenericType >::toPython(it->second,e);
    if (!ret) {
      Py_DECREF(p);
      return false;
    }
    PyDict_SetItemString(p,(it->first).c_str(),e);
    Py_DECREF(e);
  }
  return true;
}


// Explicit intialization of these two member functions, so we can use them in meta< CasADi::SX >
template<> int meta< CasADi::Matrix<CasADi::SX> >::as(GUESTOBJECT,CasADi::Matrix<CasADi::SX> &);
template<> bool meta< CasADi::Matrix<CasADi::SX> >::couldbe(GUESTOBJECT);

/// CasADi::SX
template<> char meta< CasADi::SX >::expected_message[] = "Expecting SX or number";

template <>
int meta< CasADi::SX >::as(PyObject * p,CasADi::SX &s) {
  NATIVERETURN(CasADi::SX, s)
  if (meta< double >::couldbe(p)) {
    double res;
    int result = meta< double >::as(p,res);
    if (!result)
      return false;
    s=CasADi::SX(res);
  } else {
    return false;
  }
  return true;
}

template <>
bool meta< CasADi::SX >::couldbe(PyObject * p) {
  return (meta< CasADi::SX >::isa(p) || meta< double >::couldbe(p));
}


/// CasADi::Matrix<int>
template<> char meta< CasADi::Matrix<int> >::expected_message[] = "Expecting numpy.array2D, numpy.matrix, csr_matrix, IMatrix";

template <>
int meta< CasADi::Matrix<int> >::as(PyObject * p,CasADi::Matrix<int> &m) {
  NATIVERETURN(CasADi::Matrix<int>,m)
  if (meta< int >::couldbe(p)) {
    int t;
    int res = meta< int >::as(p,t);
    m = t;
    return res;
  } else if ( meta< int >::couldbe_sequence(p)) {
    std::vector <int> t;
    int res = meta< int >::as_vector(p,t);
    m = CasADi::Matrix<int>(t,t.size(),1);
    return res;
  } else if (PyObject_HasAttrString(p,"__IMatrix__")) {
    char name[] = "__IMatrix__";
    PyObject *cr = PyObject_CallMethod(p, name,0);
    if (!cr) { return false; }
    int result = meta< CasADi::Matrix<int> >::as(cr,m);
    Py_DECREF(cr);
    return result;
  } else {
    SWIG_Error(SWIG_TypeError, "asDMatrix: unrecognised type. Should have been caught by typemap(typecheck)");
    return false;
  }
  return true;
}

// Disallow 1D numpy arrays. Allowing them may introduce conflicts with other typemaps or overloaded methods
template <>
bool meta< CasADi::Matrix<int> >::couldbe(PyObject * p) {
  return meta< int >::couldbe(p) || meta< int >::couldbe_sequence(p) || meta< CasADi::Matrix<int> >::isa(p) || PyObject_HasAttrString(p,"__IMatrix__");
}

meta_vector(CasADi::Matrix<int>)
meta_vector(std::vector< CasADi::Matrix<int> >)

/// CasADi::Matrix<double>
template<> char meta< CasADi::Matrix<double> >::expected_message[] = "Expecting numpy.array2D, numpy.matrix, csr_matrix, DMatrix";

template <>
int meta< CasADi::Matrix<double> >::as(PyObject * p,CasADi::Matrix<double> &m) {
  NATIVERETURN(CasADi::Matrix<double>,m)
  NATIVERETURN(CasADi::Matrix<int>,m)
  if (is_array(p)) { // Numpy arrays will be cast to dense Matrix<double>
    if (array_numdims(p)==0) {
      double d;
      int result = meta< double >::as(p,d);
      if (!result) return result;
      m = CasADi::Matrix<double>(d);
      return result;
    }
    if (array_numdims(p)>2 || array_numdims(p)<1) {
      std::stringstream s;
      s << "SWIG::typemapDMatrixHelper:";
      s << "Number of dimensions must be 1 or 2.";
      s << "Got " << array_numdims(p) << " instead.";
      const std::string tmp(s.str());
      const char* cstr = tmp.c_str();
      SWIG_Error_return(SWIG_TypeError,  cstr);
    }
    int nrows = array_size(p,0); // 1D array is cast into column vector
    int ncols  = 1;
    if (array_numdims(p)==2)
      ncols=array_size(p,1); 
    int size=nrows*ncols; // number of elements in the dense matrix
    if (!array_is_native(p))
      SWIG_Error_return(SWIG_TypeError, "asMatrixDouble: array byte order should be native.");
    // Make sure we have a contigous array with double datatype
    int array_is_new_object;
    PyArrayObject* array = obj_to_array_contiguous_allow_conversion(p,NPY_DOUBLE,&array_is_new_object);
    if (!array) { PyErr_Print() ; SWIG_Error_return(SWIG_TypeError, "asMatrixDouble: no luck converting numpy array to double"); }
    
    double* d=(double*) array_data(array);
    std::vector<double> v(d,d+size);
    
    m = CasADi::Matrix<double>(v, nrows, ncols);
                  
    // Free memory
    if (array_is_new_object)
      Py_DECREF(array); 
  } else if(PyObjectHasClassName(p,"csr_matrix")) { // scipy's csr_matrix will be cast to sparse Matrix<double>
    PyObject * narray=PyObject_GetAttrString( p, "data"); // need's to be decref'ed
    if (!(is_array(narray) && array_numdims(narray)==1))
      SWIG_Error_return(SWIG_TypeError, "asMatrixDouble: data should be numpy array");
    int array_is_new_object;
    PyArrayObject* array = obj_to_array_contiguous_allow_conversion(narray,NPY_DOUBLE,&array_is_new_object);
    if (!array) { PyErr_Print() ; SWIG_Error_return(SWIG_TypeError, "asMatrixDouble: no luck converting numpy array to double"); }
    int size=array_size(array,0); // number on non-zeros
    double* d=(double*) array_data(array);
    std::vector<double> v(d,d+size);

    // Get the dimensions of the csr_matrix
    PyObject * shape = PyObject_GetAttrString( p, "shape"); // need's to be decref'ed
    int nrows=PyInt_AsLong(PyTuple_GetItem(shape,0));
    int ncols=PyInt_AsLong(PyTuple_GetItem(shape,1));
		
    // Construct the 'col' vector needed for initialising the correct sparsity
    PyObject * col = PyObject_GetAttrString(p,"indices"); // need's to be decref'ed
    if (!(is_array(col) && array_numdims(col)==1 && array_type(col)==NPY_INT)) { PyErr_Print(); SWIG_Error_return(SWIG_TypeError, "asMatrixDouble: data.indices should be numpy array");}
    
    int* cold=(int*) array_data(col);
    std::vector<int> colv(cold,cold+size);
    
    // Construct the 'rowind' vector needed for initialising the correct sparsity
    PyObject * rowind = PyObject_GetAttrString(p,"indptr"); // need's to be decref'ed
    if (!(is_array(rowind) && array_numdims(rowind)==1 && array_type(rowind)==NPY_INT)) { PyErr_Print();   SWIG_Error_return(SWIG_TypeError, "asMatrixDouble: data.indptr should be numpy array");}
    int* rowindd=(int*) array_data(rowind);
    std::vector<int> rowindv(rowindd,rowindd+(nrows+1));
    
    m = CasADi::Matrix<double>(nrows,ncols,colv,rowindv, v);
    
    Py_DECREF(narray);Py_DECREF(shape);Py_DECREF(col);Py_DECREF(rowind);
    
    if (array_is_new_object)
      Py_DECREF(array);
  } else if (meta< double >::couldbe(p)) {
    double t;
    int res = meta< double >::as(p,t);
    m = t;
    return res;
  } else if ( meta< double >::couldbe_sequence(p)) {
    std::vector <double> t;
    int res = meta< double >::as_vector(p,t);
    if (t.size()>0) {
      m = CasADi::Matrix<double>(t,t.size(),1);
    } else {
      m = CasADi::Matrix<double>(t,t.size(),0);
    }
    return res;
  } else if (PyObject_HasAttrString(p,"__DMatrix__")) {
    char name[] = "__DMatrix__";
    PyObject *cr = PyObject_CallMethod(p, name,0);
    if (!cr) { return false; }
    int result = meta< CasADi::Matrix<double> >::as(cr,m);
    Py_DECREF(cr);
    return result;
  } else {
    SWIG_Error(SWIG_TypeError, "asDMatrix: unrecognised type. Should have been caught by typemap(typecheck)");
    return false;
  }
  return true;
}

// Disallow 1D numpy arrays. Allowing them may introduce conflicts with other typemaps or overloaded methods
template <>
bool meta< CasADi::Matrix<double> >::couldbe(PyObject * p) {
  return meta< double >::couldbe(p) || (((is_array(p) && array_numdims(p)==2) && array_type(p)!=NPY_OBJECT) || PyObjectHasClassName(p,"csr_matrix") || PyObjectHasClassName(p,"DMatrix")) || meta< double >::couldbe_sequence(p) || meta< CasADi::Matrix<int> >::couldbe(p) || PyObject_HasAttrString(p,"__DMatrix__");
}

meta_vector(CasADi::Matrix<double>)
meta_vector(std::vector< CasADi::Matrix<double> >)

/// CasADi::Matrix<CasADi::SX>
template<> char meta< CasADi::Matrix<CasADi::SX> >::expected_message[] = "Expecting one of: numpy.ndarray(SX/number) , SXMatrix, SX, number, sequence(SX/number)";

template <>
int meta< CasADi::Matrix<CasADi::SX> >::as(PyObject * p,CasADi::Matrix<CasADi::SX> &m) {
  NATIVERETURN(CasADi::Matrix<CasADi::SX>, m)
  NATIVERETURN(CasADi::SX, m)
  if(meta< CasADi::Matrix<double> >::couldbe(p)) {
    CasADi::DMatrix mt;
    bool result=meta< CasADi::Matrix<double> >::as(p,mt);
    if (!result)
      return false;
    m = CasADi::SXMatrix(mt);
  } else if (is_array(p)) { // Numpy arrays will be cast to dense Matrix<SX>
		if (array_type(p)!=NPY_OBJECT) {
			SWIG_Error(SWIG_TypeError, "asSXMatrix: numpy.ndarray must be of dtype object");
			return false;
		}
		if (array_numdims(p)>2 || array_numdims(p)<1) {
			SWIG_Error(SWIG_TypeError, "asSXMatrix: Number of dimensions must be 1 or 2.");
			return false;
		}
		int nrows = array_size(p,0); // 1D array is cast into column vector
		int ncols  = 1;
		if (array_numdims(p)==2)
			ncols=array_size(p,1); 
		int size=nrows*ncols; // number of elements in the dense matrix
		std::vector<CasADi::SX> v(size);
    PyArrayIterObject* it = (PyArrayIterObject*)PyArray_IterNew(p);
    PyObject *pe;
    int i=0;
		while (it->index < it->size) { 
		  pe = *((PyObject**) PyArray_ITER_DATA(it));
      bool result=meta< CasADi::SX >::as(pe,v[i++]);
      if (!result)
        return false;
		  PyArray_ITER_NEXT(it);
		}
    Py_DECREF(it);
		m = CasADi::Matrix< CasADi::SX >(v, nrows, ncols);
  } else if (PyObject_HasAttrString(p,"__SXMatrix__")) {
    char name[] = "__SXMatrix__";
    PyObject *cr = PyObject_CallMethod(p, name,0);
    if (!cr) { return false; }
    int result = meta< CasADi::Matrix<CasADi::SX> >::as(cr,m);
    Py_DECREF(cr);
    return result;
	} else {
    SWIG_Error(SWIG_TypeError, "asSXMatrix: unrecognised type. Should have been caught by typemap(typecheck)");
    return false;
  }
	return true;
}

template <>
bool meta< CasADi::Matrix<CasADi::SX> >::couldbe(PyObject * p) {
  if (is_array(p)) { // Numpy arrays will be cast to dense Matrix<SX>
    if (array_type(p)==NPY_OBJECT)
      return true;
  } else if (meta< CasADi::Matrix<double> >::couldbe(p)) {
    return true;
  }
  
  return meta< CasADi::Matrix<CasADi::SX> >::isa(p) || meta< CasADi::SX >::couldbe(p) || meta< CasADi::Matrix<double> >::couldbe(p) || PyObject_HasAttrString(p,"__SXMatrix__");
}

meta_vector(std::vector<CasADi::SX>);
meta_vector(CasADi::SX);
meta_vector(CasADi::Matrix< CasADi::SX >);
meta_vector(std::vector< CasADi::Matrix< CasADi::SX > >);

/// CasADi::MX
template<> char meta< CasADi::MX >::expected_message[] = "Expecting (MX, numberarray)";

template <>
bool meta< CasADi::MX >::couldbe(PyObject * p) {
  return meta< CasADi::MX >::isa(p) || meta< CasADi::Matrix<double> >::couldbe(p) || PyObject_HasAttrString(p,"__MX__");
}

template <>
int meta< CasADi::MX >::as(PyObject * p,CasADi::MX &m) {
  NATIVERETURN(CasADi::MX,m)
  if(meta< CasADi::Matrix<double> >::couldbe(p)) {
    CasADi::DMatrix mt;
    bool result=meta< CasADi::Matrix<double> >::as(p,mt);
    if (!result)
      return false;
    m = CasADi::MX(mt);
    return true;
  } else if (PyObject_HasAttrString(p,"__MX__")) {
    char name[] = "__MX__";
    PyObject *cr = PyObject_CallMethod(p, name,0);
    if (!cr) { return false; }
    int result = meta< CasADi::MX >::as(cr,m);
    Py_DECREF(cr);
    return result;
  }
  return false;
}

meta_vector(CasADi::MX);
meta_vector(std::vector< CasADi::MX >);
%}


