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
 return PyInt_Check(p) || PyLong_Check(p) || PyBool_Check(p) ;
}

/// double
template<> char meta< double >::expected_message[] = "Expecting double";

template <>
int meta< double >::as(PyObject * p, double &m) {
  NATIVERETURN(double,m)
  if (PyInt_Check(p) || PyBool_Check(p) || PyFloat_Check(p)) {
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
  }
}
   
template <> bool meta< double >::couldbe(PyObject * p) {
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

  return PyInt_Check(p) || PyBool_Check(p) || PyFloat_Check(p) ;
}


// Forward declarations
template<> int meta< CasADi::GenericType::Dictionary >::as(PyObject * p,CasADi::GenericType::Dictionary &s);
template<> bool meta< CasADi::GenericType::Dictionary >::couldbe(PyObject * p);
template<> bool meta< CasADi::GenericType::Dictionary >::toPython(CasADi::GenericType::Dictionary &a, PyObject *&p);


/// CasADi::GenericType
template<> char meta< CasADi::GenericType >::expected_message[] = "Expecting number, string, vector(number)";

template <>
int meta< CasADi::GenericType >::as(PyObject * p,CasADi::GenericType &s) {
  NATIVERETURN(CasADi::GenericType, s)
  if (PyBool_Check(p)) {
    s=CasADi::GenericType((bool) PyInt_AsLong(p));
  } else if (PyInt_Check(p)) {
    s=CasADi::GenericType((int) PyInt_AsLong(p));
  } else if (PyFloat_Check(p)) {
    s=CasADi::GenericType(PyFloat_AsDouble(p));
  } else if (PyString_Check(p)) {
    s=CasADi::GenericType(std::string(PyString_AsString(p)));
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
  } else if (PyType_Check(p) && PyObject_HasAttrString(p,"creator")) {
    PyObject *c = PyObject_GetAttrString(p,"creator");
    int result = meta< CasADi::GenericType >::as(c,s);
    Py_DECREF(c);
    return result;
  } else {
    PyObject* pPyObjectModuleName = PyString_FromString("casadi");
    if (!pPyObjectModuleName) { PyErr_Clear(); return false; }
    PyObject* pObjectModule = PyImport_Import(pPyObjectModuleName);
    if (!pObjectModule) { PyErr_Clear(); Py_DECREF(pPyObjectModuleName); return false; }
    PyObject* pObjectDict = PyModule_GetDict(pObjectModule);
    if (!pObjectDict) { PyErr_Clear(); Py_DECREF(pObjectModule); Py_DECREF(pPyObjectModuleName); return false; }
    PyObject* pClass = PyDict_GetItemString(pObjectDict,  "GenericType");
    if (!pClass) { PyErr_Clear(); Py_DECREF(pObjectModule); Py_DECREF(pPyObjectModuleName); return false; }
    
    PyObject* args = PyTuple_New(1);
    PyTuple_SetItem(args,0,p);
    
    PyObject* g = PyObject_CallObject(pClass,args);
    
    Py_DECREF(args);
    Py_DECREF(pObjectModule);
    Py_DECREF(pPyObjectModuleName);
    
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
  return meta< CasADi::GenericType >::isa(p) || PyBool_Check(p) ||  PyInt_Check(p) || PyFloat_Check(p) || PyString_Check(p) || meta< std::vector<double> >::couldbe(p) || ( PyType_Check(p) && PyObject_HasAttrString(p,"creator")) || PyType_Check(p) || meta< CasADi::GenericType::Dictionary >::couldbe(p);

  
  }

template <>
bool meta< CasADi::GenericType >::toPython(CasADi::GenericType &a, PyObject *&p) {
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
bool meta< CasADi::GenericType::Dictionary >::toPython(CasADi::GenericType::Dictionary &a, PyObject *&p) {
  p = PyDict_New();
  //CasADi::Dictionary::const_iterator end = a.end(); 
  CasADi::GenericType::Dictionary::iterator end = a.end();
  for (CasADi::GenericType::Dictionary::iterator it = a.begin(); it != end; ++it)
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
  } else if (meta< CasADi::Matrix< CasADi::SX > >::isa(p)) {
    CasADi::Matrix< CasADi::SX > m;
    meta< CasADi::Matrix< CasADi::SX > >::as(p,m);
    if (m.numel()==1 && m.size()==1) {
      s = m.at(0);
      return true;
    }
  } else {
    return false;
  }
  return true;
}

template <>
bool meta< CasADi::SX >::couldbe(PyObject * p) {
  if (meta< CasADi::Matrix< CasADi::SX > >::isa(p)) {
    CasADi::Matrix< CasADi::SX > m;
    meta< CasADi::Matrix< CasADi::SX > >::as(p,m);
    if (m.numel()==1 && m.size()==1)
      return true;
  }
  return (meta< CasADi::SX >::isa(p) || meta< double >::couldbe(p));
}


/// CasADi::Matrix<double>
template<> char meta< CasADi::Matrix<double> >::expected_message[] = "Expecting numpy.array2D, numpy.matrix, csr_matrix, DMatrix";

template <>
int meta< CasADi::Matrix<double> >::as(PyObject * p,CasADi::Matrix<double> &m) {
  NATIVERETURN(CasADi::Matrix<double>,m)
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

    double* d=(double*) array->data;
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
    int size=array_size(array,0); // number on non-zeros
    double* d=(double*) array->data;
    std::vector<double> v(d,d+size);

    // Get the dimensions of the csr_matrix
    PyObject * shape = PyObject_GetAttrString( p, "shape"); // need's to be decref'ed
    int nrows=PyInt_AsLong(PyTuple_GetItem(shape,0));
    int ncols=PyInt_AsLong(PyTuple_GetItem(shape,1));
		
    // Construct the 'col' vector needed for initialising the correct sparsity
    PyObject * col = PyObject_GetAttrString(p,"indices"); // need's to be decref'ed
    if (!(is_array(col) && array_numdims(col)==1 && array_type(col)==NPY_INT))
      SWIG_Error_return(SWIG_TypeError, "asMatrixDouble: data.indices should be numpy array");
    int* cold=(int*) array_data(col);
    std::vector<int> colv(cold,cold+size);
    
    // Construct the 'rowind' vector needed for initialising the correct sparsity
    PyObject * rowind = PyObject_GetAttrString(p,"indptr"); // need's to be decref'ed
    if (!(is_array(rowind) && array_numdims(rowind)==1 && array_type(rowind)==NPY_INT))
      SWIG_Error_return(SWIG_TypeError, "asMatrixDouble: data.indptr should be numpy array");
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
    m = CasADi::Matrix<double>(t,t.size(),1);
    return res;
  } else {
    SWIG_Error(SWIG_TypeError, "asDMatrix: unrecognised type. Should have been caught by typemap(typecheck)");
    return false;
  }
  return true;
}

// Disallow 1D numpy arrays. Allowing them may introduce conflicts with other typemaps or overloaded methods
template <>
bool meta< CasADi::Matrix<double> >::couldbe(PyObject * p) {
  return meta< double >::couldbe(p) || ((is_array(p) && array_numdims(p)==2) && array_type(p)!=NPY_OBJECT|| PyObjectHasClassName(p,"csr_matrix") || PyObjectHasClassName(p,"DMatrix")) || meta< double >::couldbe_sequence(p);
}

meta_vector(CasADi::Matrix<double>)

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
	} else if(meta< CasADi::SX >::couldbe_sequence(p)) {
    std::vector<CasADi::SX> sxv;
    int result = meta< CasADi::SX >::as_vector(p,sxv);
    if (result) {
      m = CasADi::SXMatrix(sxv);
    } else {
      return false;
    }
  } else if(meta< std::vector<CasADi::SX> >::couldbe_sequence(p)) {
    std::vector< std::vector<CasADi::SX> > sxv;
    int result = meta< std::vector<CasADi::SX> >::as_vector(p,sxv);
    if (result) {
      m = CasADi::SXMatrix(sxv);
    } else {
      return false;
    }
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
  
  return meta< CasADi::Matrix<CasADi::SX> >::isa(p) || meta< CasADi::SX >::couldbe(p) || meta< CasADi::Matrix<double> >::couldbe(p) || meta< CasADi::SX >::couldbe_sequence(p) || meta< std::vector< CasADi::SX > >::couldbe_sequence(p);
}

/// std::vector< CasADi::Matrix<CasADi::SX> >
#ifdef SWIGPYTHON
template<> char meta< std::vector< CasADi::Matrix<CasADi::SX> > >::expected_message[] = "Expecting sequence(numpy.ndarray(SX/number) , SXMatrix, SX, number, sequence(SX/number))";

template <>
int meta< std::vector< CasADi::Matrix<CasADi::SX> > >::as(PyObject * p,std::vector< CasADi::Matrix<CasADi::SX> > &m) {
  NATIVERETURN(std::vector< CasADi::Matrix<CasADi::SX> >,m)
  //if(isSXMatrix(p)) {
  //  CasADi::SXMatrix *mp;
  //  int result = getSXMatrix(p,mp);
  //  if (!result)
  //    return false;
  //  m.push_back(*mp);
  //}
  if( PyIsSequence(p)) {
    return meta< CasADi::Matrix<CasADi::SX> >::as_vector(p,m);
  } else if (PyDict_Check(p)) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    int num=0;
    while (PyDict_Next(p, &pos, &key, &value)) {
      if (PyInt_Check(key)) {
        if (PyInt_AsLong(key)+1>num)
          num=PyInt_AsLong(key)+1;
      } else if (PyString_Check(key) && strcmp(PyString_AsString(key),"NUM")==0) {
        if (PyInt_Check(value))
          num=PyInt_AsLong(value);
      } else {
        return false;
      }
    }
    m.resize(num);
    pos=0;
    while (PyDict_Next(p, &pos, &key, &value)) {
      if (PyInt_Check(key)) {
        bool result=meta< CasADi::Matrix< CasADi::SX > >::as(value,m[PyInt_AsLong(key)]);
        if (!result)
          return false;
      }
    }
    return true;
  } else {
    return false;
  }
  return true;
}

template <>
bool meta< std::vector< CasADi::Matrix<CasADi::SX> > >::couldbe(PyObject * p) {
  if( meta< CasADi::Matrix< CasADi::SX > >::couldbe_sequence(p)) {
    return true;
  } else if (PyDict_Check(p)) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (PyDict_Next(p, &pos, &key, &value)) {
      if (!((PyInt_Check(key) || (PyString_Check(key) && strcmp(PyString_AsString(key),"NUM")==0)) && meta< CasADi::Matrix<CasADi::SX> >::couldbe(value)))
        return false;
    }
    return true;
  }
  return meta< std::vector< CasADi::Matrix< CasADi::SX > > >::isa(p);
}

#endif //SWIGPYTHON


meta_vector(std::vector<CasADi::SX>);
meta_vector(CasADi::SX);

/// CasADi::Slice
template<> char meta< CasADi::Slice >::expected_message[] = "Expecting Slice or number";
template <>
int meta< CasADi::Slice >::as(PyObject * p,CasADi::Slice &m) {
  NATIVERETURN(CasADi::Slice,m)

  if (PyInt_Check(p)) {
    m.start_ = PyInt_AsLong(p);
    m.stop_ = m.start_+1;
    if (m.stop_==0) m.stop_ = std::numeric_limits<int>::max();
    return true;
  } else if (PySlice_Check(p)) {
    PySliceObject *r = (PySliceObject*)(p);
    if(r->start!=Py_None) m.start_ = PyInt_AsLong(r->start);
    m.stop_  = (r->stop ==Py_None) ? std::numeric_limits<int>::max() : PyInt_AsLong(r->stop) ;
    if(r->step !=Py_None) m.step_  = PyInt_AsLong(r->step);
    return true;
  } else {
    return false;
  }

}

template <>
bool meta<  CasADi::Slice >::couldbe(PyObject * p) {
  return meta< CasADi::Slice >::isa(p) || PyInt_Check(p) || PySlice_Check(p);
}

/// CasADi::IndexList
template<> char meta< CasADi::IndexList >::expected_message[] = "Expecting Slice or number or list of ints";
template <>
int meta< CasADi::IndexList >::as(PyObject * p,CasADi::IndexList &m) {
  
  if (meta< int >::couldbe(p)) {
    m.type = CasADi::IndexList::INT;
    meta< int >::as(p,m.i);
  } else if (meta< std::vector<int> >::couldbe(p)) {
    m.type = CasADi::IndexList::IVECTOR;
    return meta< std::vector<int> >::as(p,m.iv);
  } else if (meta< CasADi::Slice>::couldbe(p)) {
    m.type = CasADi::IndexList::SLICE;
    return meta< CasADi::Slice >::as(p,m.slice);
  } else {
    return false;
  }
  return true;
}


template <>
bool meta<  CasADi::IndexList >::couldbe(PyObject * p) {
  return meta< CasADi::Slice >::couldbe(p) || meta< std::vector<int> >::couldbe(p) || meta< int >::couldbe(p);
}

/// CasADi::MX
template<> char meta< CasADi::MX >::expected_message[] = "Expecting (MX, numberarray)";

template <>
bool meta< CasADi::MX >::couldbe(PyObject * p) {
  return (meta< CasADi::MX >::isa(p) || meta< CasADi::Matrix<double> >::couldbe(p) );
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
  }
  return false;
}

/// std::vector< CasADi::MX >
template<> char meta< std::vector< CasADi::MX > >::expected_message[] = "Expecting sequence(MX, number)";
template <>
int meta< std::vector< CasADi::MX > >::as(PyObject * p,std::vector< CasADi::MX > &m) {
  NATIVERETURN(std::vector< CasADi::MX >,m)
  if( PyIsSequence(p) ) {
    return meta< CasADi::MX >::as_vector(p,m);
  } else if (PyDict_Check(p)) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    int num=0;
    while (PyDict_Next(p, &pos, &key, &value)) {
      if (PyInt_Check(key)) {
        if (PyInt_AsLong(key)+1>num)
          num=PyInt_AsLong(key)+1;
      } else if (PyString_Check(key) && strcmp(PyString_AsString(key),"NUM")==0) {
        if (PyInt_Check(value))
          num=PyInt_AsLong(value);
      } else {
        return false;
      }
    }
    m.resize(num);
    pos=0;
    while (PyDict_Next(p, &pos, &key, &value)) {
      if (PyInt_Check(key)) {
        bool result=meta< CasADi::MX >::as(value,m[PyInt_AsLong(key)]);
        if (!result)
          return false;
      }
    }
    return true;
  } else {
    return false;
  }
return true;
}

template <>
bool meta< std::vector< CasADi::MX > >::couldbe(PyObject * p) {
  if(meta< CasADi::MX >::couldbe_sequence(p)) {
    return true;
  } else if (PyDict_Check(p)) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (PyDict_Next(p, &pos, &key, &value)) {
      if (!((PyInt_Check(key) || (PyString_Check(key) && strcmp(PyString_AsString(key),"NUM")==0)) && meta< CasADi::MX >::couldbe(value)))
        return false;
    }
    return true;
  }
  return meta< std::vector<CasADi::MX> >::isa(p);
}
%}

