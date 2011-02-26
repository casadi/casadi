%{
#include "casadi/matrix/crs_sparsity.hpp"
#include "casadi/matrix/matrix.hpp"
#include "casadi/sx/sx.hpp"
#include "casadi/sx/sx_tools.hpp"
%}

%include "typemaps.i"
%include "casadi/matrix/crs_sparsity.hpp"
%include "casadi/matrix/matrix.hpp"
%include "casadi/sx/sx.hpp"

%extend std::vector<CasADi::SX>{
  std::string __repr__(){ return CasADi::getRepresentation(*$self); }
  std::string __str__(){ return CasADi::getDescription(*$self); }
};

%extend std::vector<CasADi::Matrix< CasADi::SX> >{
  std::string __repr__(){ return CasADi::getRepresentation(*$self); }
  std::string __str__(){ return CasADi::getDescription(*$self); }
};

%extend std::vector<std::vector< CasADi::SX> >{
  std::string __repr__(){ return CasADi::getRepresentation(*$self); }
  std::string __str__(){ return CasADi::getDescription(*$self); }
};

namespace CasADi {
  %extend Matrix<double> {
  

    
    %pythoncode %{
    def __getitem__(self,s):
      if isinstance(s,int):
        if s < 0:
          s = s + self.size()
      elif isinstance(s,tuple):
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
      else:
        raise Exception("get/setitem expecting a tuple or int. Got %s of type %s" % (str(s),str(type(s))))
      return self.getitem(s)
    %}
    
    %pythoncode %{
    def __setitem__(self,s,val):
      if isinstance(s,int):
        if s < 0:
          s = s + self.size()
      elif isinstance(s,tuple):
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
      else:
        raise Exception("get/setitem expecting a tuple or int. Got %s of type %s" % (str(s),str(type(s))))
      self.setitem(s,val)
    %}

  };

%extend SX {
  %pythoncode %{
    def __lt__(self,other):
      return _casadi.__lt__(self,other)
    def __le__(self,other):
      return _casadi.__le__(self,other)
    def __eq__(self,other):
      return _casadi.__eq__(self,other)
    def __ne__(self,other):
      return _casadi.__ne__(self,other)
    def __gt__(self,other):
      return _casadi.__gt__(self,other)
    def __ge__(self,other):
      return _casadi.__ge__(self,other)
  %}
};
  
%extend Matrix<SX>{
    // The constructor has to be added since SX::operator Matrix<SX does not work
    // Matrix<SX>(const SX&){ *$self
    
    
    %pythoncode %{
      def __lt__(self,other):
        return _casadi.__lt__(self,other)
      def __le__(self,other):
        return _casadi.__le__(self,other)
      def __eq__(self,other):
        return _casadi.__eq__(self,other)
      def __ne__(self,other):
        return _casadi.__ne__(self,other)
      def __gt__(self,other):
        return _casadi.__gt__(self,other)
      def __ge__(self,other):
        return _casadi.__ge__(self,other)
    %}
    
    // These methods must be added since the implicit type cast does not work
    Matrix<SX> __pow__ (double b) const{ return $self->__pow__(CasADi::Matrix<CasADi::SX>(b));}
    Matrix<SX> __rpow__(double b) const{ return CasADi::Matrix<CasADi::SX>(b).__pow__(*$self);}
    Matrix<SX> __add__ (double b) const{ return *$self + CasADi::Matrix<CasADi::SX>(b);}
    Matrix<SX> __radd__(double b) const{ return CasADi::Matrix<CasADi::SX>(b) + *$self;}
    Matrix<SX> __sub__ (double b) const{ return *$self - CasADi::Matrix<CasADi::SX>(b);}
    Matrix<SX> __rsub__(double b) const{ return CasADi::Matrix<CasADi::SX>(b) - *$self;}
    Matrix<SX> __mul__ (double b) const{ return *$self * CasADi::Matrix<CasADi::SX>(b);}
    Matrix<SX> __rmul__(double b) const{ return CasADi::Matrix<CasADi::SX>(b) * *$self;}
    Matrix<SX> __div__ (double b) const{ return *$self / CasADi::Matrix<CasADi::SX>(b);}
    Matrix<SX> __rdiv__(double b) const{ return CasADi::Matrix<CasADi::SX>(b) / *$self;}

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
    
    %pythoncode %{
    def toArray(self):
      import numpy as n
      r = n.array((),dtype=object)
      r.resize(self.size1(),self.size2())
      for i in range(self.size1()):  # loop over rows
        for el in range(self.rowind(i),self.rowind(i+1)): # loop over the non-zero elements
          j=self.col(el)  # column
          r[i,j] = self[el] # add the non-zero element

      return r
    %}
    
    
  %pythoncode %{
  def __array_wrap__(self,*args,**kwargs):
    fun=getattr(self, args[1][0].__name__)
    return fun()
  %}

  %pythoncode %{
    def __array__(self,*args,**kwargs):
      import numpy as n
      if len(args) > 1 and isinstance(args[1],tuple) and isinstance(args[1][0],n.ufunc):
        return n.array([])
      else:
        return self.toArray()
  %}
        
    %pythoncode %{
    def toMatrix(self):
      import numpy as n
      return n.matrix(self.toArray())
    %}
    
    %pythoncode %{
    def __getitem__(self,s):
      if isinstance(s,int):
        if s < 0:
          s = s + self.size()
      elif isinstance(s,tuple):
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
      else:
        raise Exception("get/setitem expecting a tuple or int. Got %s of type %s" % (str(s),str(type(s))))
      return self.getitem(s)
    %}

    %pythoncode %{
    def __setitem__(self,s,val):
      if isinstance(s,int):
        if s < 0:
          s = s + self.size()
      elif isinstance(s,tuple):
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
      else:
        raise Exception("get/setitem expecting a tuple or int. Got %s of type %s" % (str(s),str(type(s))))
      self.setitem(s,val)
    %}
};
} // namespace CasADi


#ifdef WITH_NUMPY
#include <numpy/arrayobject.h>
#endif // WITH_NUMPY

// Template instantiations
%template(vector_PyObject)    std::vector<PyObject*>;


%inline %{

bool istype(PyObject *p, swig_type_info *type) {
  int res = SWIG_ConvertPtr(p, 0, type, 0);
  return SWIG_CheckState(res);
}

bool isSXMatrix(PyObject * p) {
  return istype(p,SWIGTYPE_p_CasADi__MatrixT_CasADi__SX_t);
}

bool isSXMatrixVector(PyObject * p) {
  return istype(p,SWIGTYPE_p_std__vectorT_CasADi__MatrixT_CasADi__SX_t_std__allocatorT_CasADi__MatrixT_CasADi__SX_t_t_t);
}

bool isSX(PyObject * p) {
  return istype(p,SWIGTYPE_p_CasADi__SX);
}

bool isMX(PyObject * p) {
  return istype(p,SWIGTYPE_p_CasADi__MX);
}

int getSXMatrix(PyObject * p, CasADi::SXMatrix * & m) {
  void *pd = 0 ;
  int res = SWIG_ConvertPtr(p, &pd,SWIGTYPE_p_CasADi__MatrixT_CasADi__SX_t, 0 );
  if (!SWIG_IsOK(res)) {
    return false;
  }
  m = reinterpret_cast< CasADi::SXMatrix* >(pd);
  return true;
}

int getSX_ptr(PyObject * p, CasADi::SX * & m) {
  void *pd = 0 ;
  int res = SWIG_ConvertPtr(p, &pd,SWIGTYPE_p_CasADi__SX, 0 );
  if (!SWIG_IsOK(res)) {
    return false;
  }
  m = reinterpret_cast< CasADi::SX * >(pd);
  return true;
}

int getSXMatrixVector_ptr(PyObject * p, std::vector< CasADi::Matrix<CasADi::SX> > * & m) {
  void *pd = 0 ;
  int res = SWIG_ConvertPtr(p, &pd,SWIGTYPE_p_std__vectorT_CasADi__MatrixT_CasADi__SX_t_std__allocatorT_CasADi__MatrixT_CasADi__SX_t_t_t, 0 );
  if (!SWIG_IsOK(res)) {
    return false;
  }
  m = reinterpret_cast< std::vector< CasADi::Matrix<CasADi::SX> > *  >(pd);
  return true;
}



int getPyNumber(PyObject * p, double * m) {
  PyObject *r = PyNumber_Float(p);
  int ret = !(r==NULL);
  if (ret)
    *m = PyFloat_AsDouble(r);
  Py_DECREF(r);
  return ret;
}


/**
\param p   input  SX or number
\param s    output SX
\return true if succesful
*/
bool asSX(PyObject* p, CasADi::SX & s) {
  if (isSX(p)) {
    CasADi::SX * sx;
    int result = getSX_ptr(p,sx);
    if (!result)
      return false;
    s=*sx;
  } else if (PyNumber_Check(p)) {
    double res;
    int result = getPyNumber(p,&res);
    if (!result)
      return false;
    s=CasADi::SX(res);
  } else {
    return false;
  }
  return true;
}

bool couldbeSX(PyObject* p) {
  return (isSX(p) || PyNumber_Check(p));
}

/**
\param p   input   numpy.ndarray(SX) or SXMatrix or SX or number or sequence(SX)
\param s    output SX
\return true if succesful
*/
bool asSXMatrix(PyObject*p, CasADi::SXMatrix &m ) {
if (is_array(p)) { // Numpy arrays will be cast to dense Matrix<SX>
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
      bool result=asSX(pe,v[i++]);
      if (!result)
        return false;
		  PyArray_ITER_NEXT(it);
		}
    Py_DECREF(it);
		m = CasADi::Matrix< CasADi::SX >(v, nrows, ncols);
	} else if (isSXMatrix(p)) { // SXMatrix object get passed on as-is.
    CasADi::SXMatrix *mp;
    int result = getSXMatrix(p,mp);
		if (!result)
			return false;
    m=*mp;
	} else if (isSX(p)) {
    CasADi::SX * sx;
    int result = getSX_ptr(p,sx);
		if (!result)
			return false;
    m = CasADi::Matrix< CasADi::SX >(*sx);
  } else if (PyNumber_Check(p))  {
    double res;
    int result = getPyNumber(p,&res);
		if (!result)
			return false;
    m = CasADi::Matrix< CasADi::SX >(res);
  } else if(PySequence_Check(p) &&! isSXMatrix(p) && !isMX(p)) {
    PyObject *it = PyObject_GetIter(p);
    PyObject *pe;
    std::vector<CasADi::SX> sxv(PySequence_Size(p));
    int i=0;
    while (pe = PyIter_Next(it)) {                                // Iterate over the sequence inside the sequence
      bool result=asSX(pe,sxv[i++]);
      if (!result) {
        Py_DECREF(pe);Py_DECREF(it);
        return false;
      }
      Py_DECREF(pe);
    }
    m = CasADi::SXMatrix(sxv);
    Py_DECREF(it);
  } else {
    SWIG_Error(SWIG_TypeError, "asSXMatrix: unrecognised type. Should have been caught by typemap(typecheck)");
    return false;
  }
	return true;
}


bool couldbeSXMatrix(PyObject* p) {
  if (is_array(p)) { // Numpy arrays will be cast to dense Matrix<SX>
		if (array_type(p)==NPY_OBJECT)
			return true;
	} else if(PySequence_Check(p) &&! isSXMatrix(p) && !isMX(p)) {
    PyObject *it = PyObject_GetIter(p);
    PyObject *pe;
    int i=0;
    while (pe = PyIter_Next(it)) {                                // Iterate over the sequence inside the sequence
      if (!couldbeSX(pe)) {
        Py_DECREF(pe);Py_DECREF(it);return false;
      }
      Py_DECREF(pe);
    }
    Py_DECREF(it);
    return true;
  } 
  return isSXMatrix(p) || couldbeSX(p) ;
}

/**
\param p   input   sequence(anything that SXMatrix accepts)
\param s    output SX
\return true if succesful
*/
bool asVSXMatrix(PyObject*p, std::vector<CasADi::SXMatrix> &m ) {
  //if(isSXMatrix(p)) {
  //  CasADi::SXMatrix *mp;
  //  int result = getSXMatrix(p,mp);
  //  if (!result)
  //    return false;
  //  m.push_back(*mp);
  //}
  if(PySequence_Check(p) &&! isSXMatrix(p) && !isMX(p)) {
    PyObject *it = PyObject_GetIter(p);
    PyObject *pe;
    m.resize(PySequence_Size(p));
    int i=0;
    while (pe = PyIter_Next(it)) {                                // Iterate over the sequence inside the sequence
      bool result=asSXMatrix(pe,m[i++]);
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

bool couldbeSXMatrixVector(PyObject* p) {
  if(PySequence_Check(p) &&! isSXMatrix(p) && !isMX(p)) {
    PyObject *it = PyObject_GetIter(p);
    PyObject *pe;
    int i=0;
    while (pe = PyIter_Next(it)) {                                // Iterate over the sequence inside the sequence
      if (!couldbeSXMatrix(pe)) {
        Py_DECREF(pe);Py_DECREF(it);return false;
      }
      Py_DECREF(pe);
    }
    Py_DECREF(it);
    return true;
  }
  return isSXMatrixVector(p);
}

%}

%template(SXMatrixVector)       std::vector<CasADi::Matrix<CasADi::SX> > ;
%template(SXMatrixVectorVector) std::vector< std::vector<CasADi::Matrix<CasADi::SX> > > ;

namespace CasADi{
/*
Attempts to form its argument into a std::vector<SXMatrix> form
Accepts: sequence(numpy.ndarray(SX/number) , SXMatrix, SX, number, sequence(SX/number))

matching on SXMatrix is prohibited as per wish of Joel
*/
%typemap(in) const std::vector< Matrix<SX> > &  (std::vector<CasADi::SXMatrix> v) {
  if (isSXMatrixVector($input)) { // SXMatrixVector object get passed on as-is, and fast.
    int result = getSXMatrixVector_ptr($input,$1);
		if (!result)
			SWIG_exception_fail(SWIG_TypeError,"SXMatrix cast failed");
	} else {  
    bool result=asVSXMatrix($input,v);
    if (!result)
      SWIG_exception_fail(SWIG_TypeError,"Expecting sequence(numpy.ndarray(SX/number) , SXMatrix, SX, number, sequence(SX/number))");
    $1 = &v;
  }
}

%typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) const std::vector< Matrix<SX> > & {
    if (couldbeSXMatrixVector($input)) {
      $1 = 1;
    } else {
      $1=0;
    }
}

%typemap(freearg) const std::vector< Matrix<SX> > & {
}



// numpy.ndarray(SX/number) , SXMatrix, SX, number, sequence(SX/number)
%typemap(in) const Matrix<SX> & (Matrix<SX> m) {
  if (isSXMatrix($input)) { // SXMatrix object get passed on as-is, and fast.
    int result = getSXMatrix($input,$1);
		if (!result)
			SWIG_exception_fail(SWIG_TypeError,"SXMatrix cast failed");
	} else {
    bool result=asSXMatrix($input,m);
    if (!result)
      SWIG_exception_fail(SWIG_TypeError,"Expecting one of: numpy.ndarray(SX/number) , SXMatrix, SX, number, sequence(SX/number)");
    $1 = &m;
  }
}


%typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) const Matrix<SX> & {
  if (couldbeSXMatrix($input)) {
    $1=1;
  } else {
    $1=0;
  }
}

%typemap(freearg) const Matrix<SX>  & {
}



}

%template(SXVector)             std::vector<CasADi::SX>;
%template(SXVectorVector)       std::vector<std::vector<CasADi::SX> > ;
%template(SXVectorVectorVector) std::vector< std::vector<std::vector<CasADi::SX> > > ;
%template(SXMatrix)             CasADi::Matrix<CasADi::SX>;


// NOTE: mx.py fails terminate called after throwing an instance of 'std::bad_alloc'
