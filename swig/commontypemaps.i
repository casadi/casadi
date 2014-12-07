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

#ifdef SWIGPYTHON
%fragment("to"{int}, "header") {
  int to_int(GUESTOBJECT *p, void *mv, int offs=0) {
    int *m = static_cast<int*>(mv);
    if (m) m += offs;
    return meta< int >::toCpp(p, static_cast<int *>(mv), $descriptor(int *));
  }
 }
%casadi_typemaps(int, SWIG_TYPECHECK_INTEGER, int)

%fragment("to"{double}, "header") {
  int to_double(GUESTOBJECT *p, void *mv, int offs=0) {
    double *m = static_cast<double*>(mv);
    if (m) m += offs;
    return meta< double >::toCpp(p, static_cast<double *>(mv), $descriptor(double *));
  }
 }
%casadi_typemaps(double, SWIG_TYPECHECK_DOUBLE, double)
%casadi_typemaps_constref(double, SWIG_TYPECHECK_DOUBLE, double)
#endif //SWIGPYTHON

#ifdef SWIGPYTHON
%typemap(out) casadi::GenericType {
  if(!($result = fromCpp($1)))
    SWIG_exception_fail(SWIG_TypeError,"GenericType not yet implemented");
}

%typemap(out) std::vector< casadi::GenericType > {
  PyObject* ret = PyList_New(0);
  std::vector< casadi::GenericType > & in = $1;
  for (int k=0 ; k < in.size(); ++k) {
    PyObject* rete;
    if (!(rete = fromCpp(in[k])))
      SWIG_exception_fail(SWIG_TypeError,"GenericType not yet implemented");
    PyList_Append(ret, rete);
  }
  $result = ret;
}

%typemap(out) const casadi::GenericType::Dictionary&  {
  if(!($result = fromCpp(*$1))) {
    SWIG_exception_fail(SWIG_TypeError,"GenericType not yet implemented");
  }
}
#endif // SWIGPYTHON

%fragment("to"{GenericType}, "header") {
  // Forward declaration
  int to_Dictionary(GUESTOBJECT *p, void *mv, int offs=0);

  int to_GenericType(GUESTOBJECT *p, void *mv, int offs=0) {
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
    } else if (meta< std::string >::toCpp(p, 0, *meta< std::string >::name)) {
      std::string temp;
      if (!meta< std::string >::toCpp(p, &temp, *meta< std::string >::name)) return false;
      if (m) *m = casadi::GenericType(temp);
    } else if (meta< std::vector<int> >::toCpp(p, 0, *meta< std::vector<int> >::name)) {
      std::vector<int> temp;
      if (!meta< std::vector<int> >::toCpp(p, &temp, *meta< std::vector<int> >::name)) return false;
      if (m) *m = casadi::GenericType(temp);
    } else if (meta< std::vector<double> >::toCpp(p, 0, *meta< std::vector<double> >::name)) {
      std::vector<double> temp;
      if (!meta< std::vector<double> >::toCpp(p, &temp, *meta< std::vector<double> >::name)) return false;
      if (m) *m = casadi::GenericType(temp);
    } else if (meta< std::vector<std::string> >::toCpp(p, 0, *meta< std::vector<std::string> >::name)) {
      std::vector<std::string> temp;
      if (!meta< std::vector<std::string> >::toCpp(p, &temp, *meta< std::vector<std::string> >::name)) return false;
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
    
    } else if (meta< casadi::Function >::toCpp(p, 0, *meta< casadi::Function >::name)) {
      casadi::Function temp;
      if (!meta< casadi::Function >::toCpp(p, &temp, *meta< casadi::Function >::name)) return false;
      if (m) *m = casadi::GenericType(temp);
    } else if (to_Dictionary(p, 0)
               || meta< casadi::DerivativeGenerator >::toCpp(p, 0, *meta< casadi::DerivativeGenerator >::name)
               || meta< casadi::Callback >::toCpp(p, 0, *meta< casadi::Callback >::name)) {
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
%fragment("to"{Dictionary}, "header", fragment="to"{GenericType}) {
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

%my_creator_typemap(PRECEDENCE_CREATOR, casadi::implicitFunctionCreator);
%my_creator_typemap(PRECEDENCE_CREATOR, casadi::linearSolverCreator);

#ifdef SWIGPYTHON
%fragment("to"{DerivativeGenerator}, "header") {
  int to_DerivativeGenerator(GUESTOBJECT *p, void *mv, int offs=0) {
    casadi::DerivativeGenerator *m = static_cast<casadi::DerivativeGenerator*>(mv);
    if (m) m += offs;
    return meta< casadi::DerivativeGenerator >::toCpp(p, m, $descriptor(casadi::DerivativeGenerator *));
  }
 }
%casadi_typemaps_constref(DerivativeGenerator, PRECEDENCE_DERIVATIVEGENERATOR, casadi::DerivativeGenerator)

%fragment("to"{CustomEvaluate}, "header") {
  int to_CustomEvaluate(GUESTOBJECT *p, void *mv, int offs=0) {
    casadi::CustomEvaluate *m = static_cast<casadi::CustomEvaluate*>(mv);
    if (m) m += offs;
    return meta< casadi::CustomEvaluate >::toCpp(p, m, $descriptor(casadi::CustomEvaluate *));
  }
 }
%casadi_typemaps_constref(CustomEvaluate, PRECEDENCE_CUSTOMEVALUATE, casadi::CustomEvaluate)

%fragment("to"{Callback}, "header") {
  int to_Callback(GUESTOBJECT *p, void *mv, int offs=0) {
    casadi::Callback *m = static_cast<casadi::Callback*>(mv);
    if (m) m += offs;
    return meta< casadi::Callback >::toCpp(p, m, $descriptor(casadi::Callback *));
  }
 }
%casadi_typemaps_constref(Callback, PRECEDENCE_CALLBACK, casadi::Callback)
#endif

%fragment("to"{DVector}, "header") {
  int to_DVector(GUESTOBJECT *p, void *mv, int offs=0) {
    std::vector<double> *m = static_cast<std::vector<double> *>(mv);
    if (m) m += offs;
    return meta< std::vector<double> >::toCpp(p, m, $descriptor(std::vector<double> *));
  }
 }
%casadi_typemaps_constref(DVector, PRECEDENCE_DVector, std::vector<double>)

%fragment("to"{IVector}, "header") {
  int to_IVector(GUESTOBJECT *p, void *mv, int offs=0) {
    std::vector<int> *m = static_cast<std::vector<int> *>(mv);
    if (m) m += offs;
    return meta< std::vector<int> >::toCpp(p, m, $descriptor(std::vector<int> *));
  }
 }
%casadi_typemaps_constref(IVector, PRECEDENCE_IVector, std::vector<int>)

%fragment("to"{SX}, "header") {
  int to_SX(GUESTOBJECT *p, void *mv, int offs=0) {
    casadi::SX *m = static_cast<casadi::SX*>(mv);
    if (m) m += offs;
    return meta< casadi::SX >::toCpp(p, m, $descriptor(casadi::Matrix<casadi::SXElement> *));
  }
 }
%casadi_typemaps_constref(SX, PRECEDENCE_SX, casadi::Matrix<casadi::SXElement>)
%my_genericmatrix_const_typemap(PRECEDENCE_SX,casadi::Matrix<casadi::SXElement>)
%casadi_typemaps_vector(SXVector, PRECEDENCE_SXVector, casadi::Matrix<casadi::SXElement>)
%casadi_typemaps_vector(SXVectorVector, PRECEDENCE_SXVectorVector, std::vector< casadi::Matrix<casadi::SXElement> >)

%fragment("to"{MX}, "header") {
  int to_MX(GUESTOBJECT *p, void *mv, int offs=0) {
    MX *m = static_cast<MX*>(mv);
    if (m) m += offs;
    return meta< casadi::MX >::toCpp(p, m, $descriptor(casadi::MX *));
  }
 }
%casadi_typemaps_constref(MX, PRECEDENCE_MX, casadi::MX)
%my_genericmatrix_const_typemap(PRECEDENCE_MX,casadi::MX)

%fragment("to"{DMatrix}, "header") {
  int to_DMatrix(GUESTOBJECT *p, void *mv, int offs=0) {
    casadi::DMatrix *m = static_cast<casadi::DMatrix*>(mv);
    if (m) m += offs;
    return meta< casadi::DMatrix >::toCpp(p, m, $descriptor(casadi::Matrix<double> *));
  }
 }
%casadi_typemaps_constref(DMatrix, PRECEDENCE_DMatrix, casadi::Matrix<double>)
%my_genericmatrix_const_typemap(PRECEDENCE_DMatrix,casadi::Matrix<double>)
%casadi_typemaps_vector(MXVector, PRECEDENCE_MXVector, casadi::MX)

#ifdef SWIGPYTHON
%fragment("to"{IMatrix}, "header") {
  int to_IMatrix(GUESTOBJECT *p, void *mv, int offs=0) {
    casadi::IMatrix *m = static_cast<casadi::IMatrix*>(mv);
    if (m) m += offs;
    return meta< casadi::IMatrix >::toCpp(p, m, $descriptor(casadi::Matrix<int> *));
  }
 }
%casadi_typemaps_constref(IMatrix, PRECEDENCE_IMatrix, casadi::Matrix<int>)
%my_genericmatrix_const_typemap(PRECEDENCE_IMatrix,casadi::Matrix<int>)

%casadi_typemaps_vector(MXVectorVector, PRECEDENCE_MXVectorVector, std::vector<casadi::MX>)
%casadi_typemaps_vector(DMatrixVector, PRECEDENCE_DMatrixVector, casadi::Matrix<double>)
%casadi_typemaps_vector(IMatrixVector, PRECEDENCE_IMatrixVector, casadi::Matrix<int>)
%casadi_typemaps_vector(DMatrixVectorVector, PRECEDENCE_DMatrixVectorVector, std::vector< casadi::Matrix<double> >)
%casadi_typemaps_vector(IMatrixVectorVector, PRECEDENCE_IMatrixVectorVector, std::vector< casadi::Matrix<int> >)
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
