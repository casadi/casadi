%{
#include "casadi/mx/mx.hpp"
#include "casadi/mx/mx_tools.hpp"
%}

%include "casadi/mx/mx.hpp"

#ifdef SWIGPYTHON

%inline %{

int getMXVector_ptr(PyObject * p, std::vector<CasADi::MX> * & m) {
  void *pd = 0 ;
  int res = SWIG_ConvertPtr(p, &pd,SWIGTYPE_p_std__vectorT_CasADi__MX_std__allocatorT_CasADi__MX_t_t, 0 );
  if (!SWIG_IsOK(res)) {
    return false;
  }
  m = reinterpret_cast< std::vector<CasADi::MX> * >(pd);
  return true;
}

int getMX_ptr(PyObject * p, CasADi::MX * & m) {
  void *pd = 0 ;
  int res = SWIG_ConvertPtr(p, &pd,SWIGTYPE_p_CasADi__MX, 0 );
  if (!SWIG_IsOK(res)) {
    return false;
  }
  m = reinterpret_cast< CasADi::MX* >(pd);
  return true;
}

bool couldbeMX(PyObject* p) {
  return (isMX(p) || couldbePyNumber(p) || couldbeDMatrix(p) );
}

bool asMX(PyObject*p, CasADi::MX &m ) {
  if (isMX(p)) {
    CasADi::MX * mx;
    int result = getMX_ptr(p,mx);
    if (!result)
      return false;
    m=*mx;
  } else if(couldbeDMatrix(p)) {
    if (isDMatrix(p)) { 
      CasADi::DMatrix * mt;
      int result = getDMatrix_ptr(p,mt);
		  if (!result)
			  return false;
      m = CasADi::MX(*mt);
	  } else {  
	    CasADi::DMatrix mt;
      bool result=asDMatrix(p,mt);
      if (!result)
        return false;
      m = CasADi::MX(mt);
    }
  } else if (couldbePyNumber(p)) {
    double res;
    int result = getPyNumber(p,&res);
    if (!result)
      return false;
    m=CasADi::MX(res);
  } else {
    return false;
  }
  return true;
}

bool asMXVector(PyObject*p, std::vector<CasADi::MX> &m ) {
  if(PySequence_Check(p) &&! isSXMatrix(p) && !isMX(p)) {
    PyObject *it = PyObject_GetIter(p);
    PyObject *pe;
    m.resize(PySequence_Size(p));
    int i=0;
    while (pe = PyIter_Next(it)) {                                // Iterate over the sequence inside the sequence
      bool result=asMX(pe,m[i++]);
      if (!result) {
        Py_DECREF(pe);Py_DECREF(it);
        return false;
      }
      Py_DECREF(pe);
    }
    Py_DECREF(it);
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
        bool result=asMX(value,m[PyInt_AsLong(key)]);
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


bool couldbeMXVector(PyObject* p) {
  if(PySequence_Check(p) &&! isSXMatrix(p) && !isMX(p)) {
    PyObject *it = PyObject_GetIter(p);
    PyObject *pe;
    int i=0;
    while (pe = PyIter_Next(it)) {                                // Iterate over the sequence inside the sequence
      if (!couldbeMX(pe)) {
        Py_DECREF(pe);Py_DECREF(it);return false;
      }
      Py_DECREF(pe);
    }
    Py_DECREF(it);
    return true;
  } else if (PyDict_Check(p)) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (PyDict_Next(p, &pos, &key, &value)) {
      if (!((PyInt_Check(key) || (PyString_Check(key) && strcmp(PyString_AsString(key),"NUM")==0)) && couldbeMX(value)))
        return false;
    }
    return true;
  }
  return isMXVector(p);
}

int getDMatrixVector_ptr(PyObject * p, std::vector<CasADi::DMatrix> * & m) {
  void *pd = 0 ;
  int res = SWIG_ConvertPtr(p, &pd,SWIGTYPE_p_std__vectorT_CasADi__MatrixT_double_t_std__allocatorT_CasADi__MatrixT_double_t_t_t, 0 );
  if (!SWIG_IsOK(res)) {
    return false;
  }
  m = reinterpret_cast< std::vector<CasADi::DMatrix> * >(pd);
  return true;
}

bool couldbeDMatrixVector(PyObject* p) {
  if(PySequence_Check(p) &&! isSXMatrix(p) && !isMX(p)) {
    PyObject *it = PyObject_GetIter(p);
    PyObject *pe;
    int i=0;
    while (pe = PyIter_Next(it)) {                                // Iterate over the sequence inside the sequence
      if (!couldbeDMatrix(pe)) {
        Py_DECREF(pe);Py_DECREF(it);return false;
      }
      Py_DECREF(pe);
    }
    Py_DECREF(it);
    return true;
  }
  return isDMatrixVector(p);
}

bool asDMatrixVector(PyObject*p, std::vector<CasADi::Matrix<double> > &m ) {
  if(PySequence_Check(p) &&! isSXMatrix(p) && !isMX(p)) {
    PyObject *it = PyObject_GetIter(p);
    PyObject *pe;
    m.resize(PySequence_Size(p));
    int i=0;
    while (pe = PyIter_Next(it)) {                                // Iterate over the sequence inside the sequence
      bool result=asDMatrix(pe,m[i++]);
      if (!result) {
        Py_DECREF(pe);Py_DECREF(it);
        return false;
      }
      Py_DECREF(pe);
    }
    Py_DECREF(it);
  } else {
    return false;
  }
  return true;
}

%}

namespace CasADi{

%typemap(in) const std::vector< MX > &  (std::vector<MX> v) {
  if (isMXVector($input)) { // MXVector object get passed on as-is, and fast.
    int result = getMXVector_ptr($input,$1);
    if (!result)
      SWIG_exception_fail(SWIG_TypeError,"MX cast failed");
  } else {  
    bool result=asMXVector($input,v);
    if (!result)
      SWIG_exception_fail(SWIG_TypeError,"Expecting sequence(MX, number)");
    $1 = &v;
  }
}

%typemap(typecheck,precedence=PRECEDENCE_MXVector) const std::vector< MX > & { $1 = couldbeMXVector($input); }
%typemap(freearg) const std::vector< MX > & {}



%typemap(in) const std::vector< Matrix<double> > &  (std::vector<CasADi::DMatrix> v) {
  if (isDMatrixVector($input)) { // DMatrixVector object get passed on as-is, and fast.
    int result = getDMatrixVector_ptr($input,$1);
    if (!result)
      SWIG_exception_fail(SWIG_TypeError,"DMatrixVector cast failed");
  } else {  
    bool result=asDMatrixVector($input,v);
    if (!result)
      SWIG_exception_fail(SWIG_TypeError,"Expecting sequence(DMatrix)");
    $1 = &v;
  }
}

%typemap(typecheck,precedence=PRECEDENCE_DMatrixVector) const std::vector< Matrix<double> > &  { $1 = couldbeDMatrixVector($input); }
%typemap(freearg) const std::vector< Matrix<double> > &  {}

}

#endif // SWIGPYTHON

%include "casadi/mx/mx_tools.hpp"

#ifdef SWIGPYTHON
namespace CasADi{
  %extend MX{
    %pythoncode %{
    def __getitem__(self,s):
      if isinstance(s,int):
        if s < 0:
          s = s + self.size()
      elif isinstance(s,list):
        return self.getNZ(s)
      elif isinstance(s,slice):
        s = (s,[0])
      if isinstance(s,tuple):
        if len(s)!=2:
          raise Exception("get/setitem can only do 1D or 2D indexing")
        s = list(s)
        if isinstance(s[0],int) and isinstance(s[1],int):
          for k in range(2):
            if s[k]<0:
              s[k]=s[k]+self.shape[k]
        else:
          for k in range(2):
            if isinstance(s[k],slice):
              J = s[k].indices(self.shape[k])
              s[k] = range(J[0],J[1],J[2])
            elif isinstance(s[k],int):
              if s[k]<0:
                s[k]=s[k]+self.shape[k]
              s[k] = [s[k]]
      return self.getitem(s)
    %}

    %pythoncode %{
    def __setitem__(self,s,val):
      if isinstance(s,int):
        if s < 0:
          s = s + self.size()
      elif isinstance(s,list):
        self.setNZ(s,val)
        return
      elif isinstance(s,slice):
        s = (s,[0])
      if isinstance(s,tuple):
        if len(s)!=2:
          raise Exception("get/setitem can only do 1D or 2D indexing")
        s = list(s)
        if isinstance(s[0],int) and isinstance(s[1],int):
          for k in range(2):
            if s[k]<0:
              s[k]=s[k]+self.shape[k]
        else:
          for k in range(2):
            if isinstance(s[k],slice):
              J = s[k].indices(self.shape[k])
              s[k] = range(J[0],J[1],J[2])
            elif isinstance(s[k],int):
              if s[k]<0:
                s[k]=s[k]+self.shape[k]
              s[k] = [s[k]]
      self.setitem(s,val)
    %}


    %pythoncode %{
        @property
        def shape(self):
            return (self.size1(),self.size2())
    %}
    
    %pythoncode %{
        @property
        def T(self):
            return trans(self)
    %}

  };
} // namespace CasADi
#endif // SWIGPYTHON

