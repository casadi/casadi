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


// Lower value means wil be checked first
#define PRECEDENCE_IVector 92
#define PRECEDENCE_PAIR_SLICE_SLICE 93
// Why are SLICE and IndexVector the same precedence?
// To circumvent an issue with typemap precedence.
// Originally, we had slice precedence < IndexList precedence, but this caused the following order:
//    indexed(Slice,Slice)
//    indexed(Slice,Martix<int>)
//    indexed(IndexList,IndexList)
//    indexed(IndexList,Martix<int>)
// While we intend it to be:
//    indexed(Slice,Slice)
//    indexed(IndexList,IndexList)
//    indexed(Slice,Martix<int>)
//    indexed(IndexList,Martix<int>)
#define PRECEDENCE_SLICE 94
#define PRECEDENCE_IndexVector 94
#define PRECEDENCE_PAIR_IVector_IVector 96

#define PRECEDENCE_IMatrix 97
#define PRECEDENCE_IMatrixVector 98
#define PRECEDENCE_IMatrixVectorVector 98

#define PRECEDENCE_DVector 99

#define PRECEDENCE_DMatrix 100
#define PRECEDENCE_DMatrixVector 101
#define PRECEDENCE_DMatrixVectorVector 101

#define PRECEDENCE_SX 103
#define PRECEDENCE_SXVector 103
#define PRECEDENCE_SXVectorVector 103
#define PRECEDENCE_MX 104
#define PRECEDENCE_MXVector 105
#define PRECEDENCE_MXVectorVector 106

#define PRECEDENCE_CREATOR 150
#define PRECEDENCE_DERIVATIVEGENERATOR 21
#define PRECEDENCE_CUSTOMEVALUATE 21
#define PRECEDENCE_CALLBACK 21

#define PRECEDENCE_GENERICTYPE 22
#define PRECEDENCE_DICTIONARY 21

%template(SparsityVector) std::vector< casadi::Sparsity > ;
%template(SparsityVectorVector) std::vector< std::vector< casadi::Sparsity> > ;
%template(SXVector) std::vector<casadi::Matrix<casadi::SXElement> > ;
%template(SXVectorVector) std::vector< std::vector<casadi::Matrix<casadi::SXElement> > > ;
%template(MXVector) std::vector<casadi::MX>;
%template(MXVectorVector) std::vector< std::vector<casadi::MX> >;
%template(IMatrixVector) std::vector<casadi::Matrix<int> > ;
%template(DMatrixVector) std::vector<casadi::Matrix<double> > ;
%template(DMatrixVectorVector) std::vector< std::vector<casadi::Matrix<double> > > ;
%template(IMatrixVectorVector) std::vector< std::vector<casadi::Matrix<int> > > ;


%fragment("to"{int}, "header", fragment="fwd") {
  int to_int(GUESTOBJECT *p, void *mv, int offs) {
    int *m = static_cast<int*>(mv);
    if (m) m += offs;

    // Check if built-in int
    if (SWIG_IsOK(SWIG_AsVal_int(p, m))) return true;

    // Scalar IMatrix
    if (is_a(p, type_IMatrix())) {
      casadi::Matrix<int> *temp;
      SWIG_ConvertPtr(p, (void **) &temp, type_IMatrix(), 0 );
      if (temp->isScalar(true)) {
        if (m) *m = temp->toScalar();
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
  int to_double(GUESTOBJECT *p, void *mv, int offs) {
    double *m = static_cast<double*>(mv);
    if (m) m += offs;

    // Check if built-in double
    if (SWIG_IsOK(SWIG_AsVal_double(p, m))) return true;

    // Scalar DMatrix
    if (is_a(p, type_DMatrix())) {
      casadi::Matrix<double> *ptr;
      SWIG_ConvertPtr(p, (void **) &ptr, type_DMatrix(), 0 );
      if (ptr->isScalar(true)) {
        if (m) *m = ptr->toScalar();
        return true;
      }
      return false;
    }

    // Scalar IMatrix
    if (is_a(p, type_IMatrix())) {
      casadi::Matrix<int> *ptr;
      SWIG_ConvertPtr(p, (void **) &ptr, type_IMatrix(), 0 );
      if (ptr->isScalar(true)) {
        if (m) *m = ptr->toScalar();
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
    } else if (a.isDoubleVector()) {
      p = swig::from( a.toDoubleVector());
    }  else if (a.isStringVector()) {
      p = swig::from(a.toStringVector());
    } else if (a.isDictionary()) {
      p = from_Dictionary(a.toDictionary());
    } else if (a.isNull()) {
      p = Py_None;
    }
    return p;
  }
}

%typemap(out, fragment="from"{GenericType}) casadi::GenericType {
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

%fragment("from"{Dictionary}, "header", fragment="fwd") {
  GUESTOBJECT * from_Dictionary(const casadi::GenericType::Dictionary &a) {
    PyObject *p = PyDict_New();
    casadi::GenericType::Dictionary::const_iterator end = a.end();
    for (casadi::GenericType::Dictionary::const_iterator it = a.begin(); it != end; ++it) {
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
%typemap(out, fragment="from"{Dictionary}) const casadi::GenericType::Dictionary&  {
  if(!($result = from_Dictionary(*$1))) SWIG_exception_fail(SWIG_TypeError,"GenericType not yet implemented");
}
#endif // SWIGPYTHON

%fragment("to"{string}, "header", fragment="fwd") {
  int to_string(GUESTOBJECT *p, void *mv, int offs) {
    std::string *m = static_cast<std::string*>(mv);
    if (m) m += offs;
#ifdef SWIGPYTHON
    std::string *mp = 0;
    if (!is_null(p) && SWIG_ConvertPtr(p, (void **) &mp, $descriptor(std::string *), 0) != -1) {
      if (m) *m=*mp;
      return true;
    }
    if (!PyString_Check(p)) return false;
    if (m) m->clear();
    if (m) m->append(PyString_AsString(p));
    return true;
#endif // SWIGPYTHON
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
    Function DerivativeGeneratorPythonInternal::call(Function& fcn, int nfwd, int nadj, void* user_data) {
      casadi_assert(p_!=0);
      PyObject * nfwd_py = PyInt_FromLong(nfwd);
      PyObject * nadj_py = PyInt_FromLong(nadj);
      PyObject * fcn_py = SWIG_NewPointerObj((new Function(static_cast< const Function& >(fcn))),
                                             $descriptor(casadi::Function *), SWIG_POINTER_OWN |  0 );
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
    } else if (to_Dictionary(p, 0) || to_DerivativeGenerator(p, 0) || to_Callback(p, 0)) {
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
%fragment("to"{Dictionary}, "header", fragment="fwd") {
  int to_Dictionary(GUESTOBJECT *p, void *mv, int offs) {
    casadi::GenericType::Dictionary *m = static_cast<casadi::GenericType::Dictionary*>(mv);
    if (m) m += offs;
#ifndef SWIGPYTHON
    return false;
#else // SWIGPYTHON
    casadi::GenericType::Dictionary *mp = 0;
    if (p != Py_None && SWIG_ConvertPtr(p, (void **) &mp, $descriptor(casadi::GenericType::Dictionary *), 0) != -1) {
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
%casadi_typemaps_constref(Dictionary, PRECEDENCE_DICTIONARY, casadi::GenericType::Dictionary)
#endif

%fragment("to"{SX}, "header", fragment="fwd") {
  int to_SX(GUESTOBJECT *p, void *mv, int offs) {
    casadi::SX *m = static_cast<casadi::SX*>(mv);
    if (m) m += offs;
    // Use operator=
    if (TRY_COPY(p, casadi::SX, type_SX(), m)) return true;
    if (TRY_COPY(p, casadi::DMatrix, type_DMatrix(), m)) return true;
    if (TRY_COPY(p, casadi::IMatrix, type_IMatrix(), m)) return true;
    // Try first converting to a temporary DMatrix
    {
      casadi::DMatrix mt;
      if(to_DMatrix(p, m ? &mt : 0)) {
        if (m) *m = mt;
        return true;
      }
    }
#ifdef SWIGPYTHON
    // Numpy arrays will be cast to dense SX
    if (is_array(p)) {
      if (array_type(p) != NPY_OBJECT) return false;
      if (array_numdims(p)>2 || array_numdims(p)<1) return false;
      PyArrayIterObject* it = (PyArrayIterObject*)PyArray_IterNew(p);
      std::vector<casadi::SXElement> v;
      if (m) v.resize(it->size);
      std::vector<casadi::SXElement>::iterator v_it = v.begin();
      casadi::SX tmp;
      PyObject *pe;
      while (it->index < it->size) { 
        pe = *((PyObject**) PyArray_ITER_DATA(it));
        if (!to_SX(pe, &tmp)) return false;
        if (!tmp.isScalar()) return false;
        if (m) *v_it++ = tmp.toScalar();
        PyArray_ITER_NEXT(it);
      }
      Py_DECREF(it);
      int nrows = array_size(p,0); // 1D array is cast into column vector
      int ncols  = array_numdims(p)==2 ? array_size(p,1) : 1;
      if (m) {
        *m = casadi::SX::zeros(nrows, ncols);
        m->set(v, casadi::DENSETRANS);
      }
      return true;
    }
    // Object has __SX__ method
    if (PyObject_HasAttrString(p,"__SX__")) {
      char name[] = "__SX__";
      PyObject *cr = PyObject_CallMethod(p, name, 0);
      if (!cr) return false;
      int result = to_SX(cr, m);
      Py_DECREF(cr);
      return result;
    }
#endif // SWIGPYTHON
    return false;
  }
 }
%casadi_typemaps_constref(SX, PRECEDENCE_SX, casadi::Matrix<casadi::SXElement>)
%casadi_typemaps_genericmatrix(SX, PRECEDENCE_SX, casadi::Matrix<casadi::SXElement>)
%casadi_typemaps_vector(SX, PRECEDENCE_SXVector, casadi::Matrix<casadi::SXElement>)
%casadi_typemaps_vector2(SX, PRECEDENCE_SXVectorVector, casadi::Matrix<casadi::SXElement>)

%fragment("to"{MX}, "header", fragment="fwd") {
  int to_MX(GUESTOBJECT *p, void *mv, int offs) {
    MX *m = static_cast<MX*>(mv);
    if (m) m += offs;
    // Use operator=
    if (TRY_COPY(p, casadi::MX, type_MX(), m)) return true;
    if (TRY_COPY(p, casadi::DMatrix, type_DMatrix(), m)) return true;
    // Try first converting to a temporary DMatrix
    {
      casadi::DMatrix mt;
      if(to_DMatrix(p, m ? &mt : 0)) {
        if (m) *m = mt;
        return true;
      }
    }
#ifdef SWIGPYTHON
    if (PyObject_HasAttrString(p,"__MX__")) {
      char name[] = "__MX__";
      PyObject *cr = PyObject_CallMethod(p, name,0);
      if (!cr) return false;
      int result = to_MX(cr, m);
      Py_DECREF(cr);
      return result;
    }
#endif // SWIGPYTHON
    // Failure if reached this point
    return false;
  }
 }
%casadi_typemaps_constref(MX, PRECEDENCE_MX, casadi::MX)
%casadi_typemaps_genericmatrix(MX, PRECEDENCE_MX, casadi::MX)

%fragment("to"{DMatrix}, "header", fragment="fwd,make_vector", fragment="to"{double}) {
  int to_DMatrix(GUESTOBJECT *p, void *mv, int offs) {
    casadi::DMatrix *m = static_cast<casadi::DMatrix*>(mv);
    if (m) m += offs;
    // Use operator=
    if (TRY_COPY(p, casadi::DMatrix, type_DMatrix(), m)) return true;
    if (TRY_COPY(p, casadi::IMatrix, type_IMatrix(), m)) return true;
    // If double
    {
      double t;
      if (to_double(p, &t)) {
        if (m) *m = t;
        return true;
      }
    }
#ifdef SWIGPYTHON
    // Object has __DMatrix__ method
    if (PyObject_HasAttrString(p,"__DMatrix__")) {
      char name[] = "__DMatrix__";
      PyObject *cr = PyObject_CallMethod(p, name, 0);
      if (!cr) return false;
      int result = to_DMatrix(cr, m);
      Py_DECREF(cr);
      return result;
    }
    // Numpy arrays will be cast to dense Matrix<double>
    if (is_array(p)) { 
      if (array_numdims(p)==0) {
        double d;
        int result = to_double(p, &d);
        if (!result) return result;
        if (m) *m = casadi::Matrix<double>(d);
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
        *m = casadi::Matrix<double>::zeros(nrows,ncols);
        m->set(d, casadi::DENSETRANS);
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
      
        if (m) *m = casadi::Matrix<double>(casadi::Sparsity(nrows,ncols,colindv,rowv), v);
      
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
      int result = to_DMatrix(cr, m);
      Py_DECREF(cr);
      return result;
    }
    {
      std::vector <double> t;
      int res = make_vector(p, &t, to_double);
      if (t.size()>0) {
        if (m) *m = casadi::Matrix<double>(t,t.size(),1);
      } else {
        if (m) *m = casadi::Matrix<double>(t,t.size(),0);
      }
      return res;
    }
#endif // SWIGPYTHON
#ifdef SWIGMATLAB
    // MATLAB dense matrix
    if (mxIsDouble(p)) {
      if (m) {
        // Get dimensions
        size_t nrow = mxGetM(p);
        size_t ncol = mxGetN(p);
        double* data = static_cast<double*>(mxGetData(p));
        *m = casadi::DMatrix::zeros(nrow,ncol);
        m->set(data);
      }
      return true;
    }
#endif // SWIGMATLAB
    return false;
  }
}
%casadi_typemaps_constref(DMatrix, PRECEDENCE_DMatrix, casadi::Matrix<double>)
%casadi_typemaps_genericmatrix(DMatrix, PRECEDENCE_DMatrix, casadi::Matrix<double>)
%casadi_typemaps_vector(MX, PRECEDENCE_MXVector, casadi::MX)

#ifdef SWIGPYTHON
%fragment("to"{IMatrix}, "header", fragment="fwd,make_vector") {
  int to_IMatrix(GUESTOBJECT *p, void *mv, int offs) {
    casadi::IMatrix *m = static_cast<casadi::IMatrix*>(mv);
    // Use operator=
    if (m) m += offs;
    if (TRY_COPY(p, casadi::IMatrix, type_IMatrix(), m)) return true;
    // If scalar
    {
      int t;
      if (to_int(p, &t)) {
        if (m) *m = t;
        return true;
      }
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
      int size=nrows*ncols; // number of elements in the dense matrix
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
          *m = m->trans();                  
        }
        // Free memory
        if (array_is_new_object) Py_DECREF(array);
        return true;
      }
    
      int* d=(int*) array_data(array);
      std::vector<int> v(d,d+size);
    
      if (m) {
        *m = casadi::Matrix<int>::zeros(nrows,ncols);
        m->set(d,casadi::DENSETRANS);
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
      if (m) *m = casadi::Matrix<int>(t,t.size(),1);
      return res;
    }
    return true;
#endif // SWIGPYTHON
#ifdef SWIGMATLAB
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
#endif // SWIGPYTHON

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
