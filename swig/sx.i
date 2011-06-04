%{
#include "casadi/matrix/crs_sparsity.hpp"
#include "casadi/matrix/slice.hpp"
#include "casadi/matrix/matrix.hpp"
#include "casadi/sx/sx.hpp"
#include "casadi/sx/sx_tools.hpp"
%}


#ifdef SWIGOCTAVE
%rename(__el_mul__) __mul__;
%rename(__el_div__) __div__;
%rename(__mul__) prod;
%rename(__transpose__) trans;
#endif // SWIGOCTAVE

%include "typemaps.i"
%include "casadi/matrix/crs_sparsity.hpp"
%include "casadi/matrix/slice.hpp"
%include "casadi/matrix/matrix.hpp"
%include "casadi/sx/sx.hpp"

#ifdef SWIGPYTHON
%extend CasADi::CRSSparsity{
    %pythoncode %{
        @property
        def shape(self):
            return (self.size1(),self.size2())
    %}
};
#endif // SWIGPYTHON

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

#ifdef SWIGPYTHON
%pythoncode %{
try:
  import numpy
  constpow=numpy.frompyfunc(lambda x,y: x.constpow(y),2,1)
except:
  pass
%}
#endif // SWIGPYTHON

namespace CasADi {

#ifdef SWIGPYTHON
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
    # define binopsT(T) \
    Matrix<SX> __pow__ (T) const{ return $self->__pow__(CasADi::Matrix<CasADi::SX>(b));} \
    Matrix<SX> __rpow__(T) const{ return CasADi::Matrix<CasADi::SX>(b).__pow__(*$self);} \
    Matrix<SX> __add__ (T) const{ return *$self + CasADi::Matrix<CasADi::SX>(b);} \
    Matrix<SX> __radd__(T) const{ return CasADi::Matrix<CasADi::SX>(b) + *$self;} \
    Matrix<SX> __sub__ (T) const{ return *$self - CasADi::Matrix<CasADi::SX>(b);} \
    Matrix<SX> __rsub__(T) const{ return CasADi::Matrix<CasADi::SX>(b) - *$self;} \
    Matrix<SX> __mul__ (T) const{ return *$self * CasADi::Matrix<CasADi::SX>(b);} \
    Matrix<SX> __rmul__(T) const{ return CasADi::Matrix<CasADi::SX>(b) * *$self;} \
    Matrix<SX> __div__ (T) const{ return *$self / CasADi::Matrix<CasADi::SX>(b);} \
    Matrix<SX> __rdiv__(T) const{ return CasADi::Matrix<CasADi::SX>(b) / *$self;}

    binopsT(double b)
    binopsT(const Matrix<double>& b)
    
    #undef binopsT
    

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
  __array_priority__ = 1000
  %}
    
  %pythoncode %{
  def __array_wrap__(self,out_arr,context=None):
        name = context[0].__name__
        conversion = {"multiply": "mul", "divide": "div", "subtract":"sub","power":"pow"}
        if name in conversion:
          name = conversion[name]
        if len(context[1])==2 and context[1][1] is self:
          name = 'r' + name
        if not(hasattr(self,name)):
          name = '__' + name + '__'
        fun=getattr(self, name)
        return fun(*context[1][0:-1])
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
};

#endif // SWIGPYTHON

} // namespace CasADi

#ifdef SWIGPYTHON
#ifdef WITH_NUMPY
#include <numpy/arrayobject.h>
#endif // WITH_NUMPY

// Template instantiations
%template(vector_PyObject)    std::vector<PyObject*>;


%inline %{

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

bool isMXVector(PyObject * p) {
  return istype(p,SWIGTYPE_p_std__vectorT_CasADi__MX_std__allocatorT_CasADi__MX_t_t);
}

bool isDMatrixVector(PyObject * p) {
  return istype(p,SWIGTYPE_p_std__vectorT_CasADi__MatrixT_double_t_std__allocatorT_CasADi__MatrixT_double_t_t_t);
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

bool couldbePyNumber(PyObject * p) {
  return PyInt_Check(p) || PyBool_Check(p) || PyFloat_Check(p);
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
  } else if (couldbePyNumber(p)) {
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
  return (isSX(p) || couldbePyNumber(p));
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
  } else if (couldbePyNumber(p))  {
    double res;
    int result = getPyNumber(p,&res);
		if (!result)
			return false;
    m = CasADi::Matrix< CasADi::SX >(res);
  } else if(couldbeDMatrix(p)) {
    if (isDMatrix(p)) { 
      CasADi::DMatrix * mt;
      int result = getDMatrix_ptr(p,mt);
		  if (!result)
			  return false;
      m = CasADi::SXMatrix(*mt);
	  } else {  
	    CasADi::DMatrix mt;
      bool result=asDMatrix(p,mt);
      if (!result)
        return false;
      m = CasADi::SXMatrix(mt);
    }
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
  } else if (couldbeDMatrix(p)) {
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
        bool result=asSXMatrix(value,m[PyInt_AsLong(key)]);
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
  } else if (PyDict_Check(p)) {
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (PyDict_Next(p, &pos, &key, &value)) {
      if (!((PyInt_Check(key) || (PyString_Check(key) && strcmp(PyString_AsString(key),"NUM")==0)) && couldbeSXMatrix(value)))
        return false;
    }
    return true;
  }
  return isSXMatrixVector(p);
}

%}

#endif // SWIGPYTHON

#ifdef SWIGOCTAVE
%inline %{
  
bool isSXMatrix(const octave_value& p){
  return istype(p,SWIGTYPE_p_CasADi__MatrixT_CasADi__SX_t);
}

bool couldbeSXMatrix(const octave_value& p){
  return p.is_real_matrix();
}

bool getSXMatrix(const octave_value& p, CasADi::SXMatrix * & m){
  void *pd = 0;
  int res = SWIG_ConvertPtr(p, &pd,SWIGTYPE_p_CasADi__MatrixT_CasADi__SX_t, 0 );
  m = reinterpret_cast< CasADi::SXMatrix * >(pd);
  return SWIG_IsOK(res);
}

bool asSXMatrix(const octave_value& p, CasADi::SXMatrix &m){
  if(p.is_real_matrix()){
    Matrix mat = p.matrix_value();
    m = CasADi::SXMatrix(mat.rows(),mat.cols(),0);
    for(int i=0; i<mat.rows(); ++i){
      for(int j=0; j<mat.cols(); ++j){
        m(i,j) = mat(i,j);
      }
    }
  }
    
  return true;
}

bool isSXMatrixVector(const octave_value& p){
  return istype(p,SWIGTYPE_p_std__vectorT_CasADi__MatrixT_CasADi__SX_t_std__allocatorT_CasADi__MatrixT_CasADi__SX_t_t_t);
}

int getSXMatrixVector_ptr(const octave_value& p, std::vector< CasADi::Matrix<CasADi::SX> > * & m) {
  void *pd = 0 ;
  int res = SWIG_ConvertPtr(p, &pd,SWIGTYPE_p_std__vectorT_CasADi__MatrixT_CasADi__SX_t_std__allocatorT_CasADi__MatrixT_CasADi__SX_t_t_t, 0 );
  m = reinterpret_cast< std::vector< CasADi::Matrix<CasADi::SX> > *  >(pd);
  return SWIG_IsOK(res);
}


bool couldbeSXMatrixVector(const octave_value& p){
  return p.is_cell();
}

bool asVSXMatrix(const octave_value& p, std::vector<CasADi::SXMatrix> &m ) {
  int nrow = p.rows();
  int ncol = p.columns();
  if(nrow != 1) return false;
  m.resize(ncol);
  
  for(int i=0; i<ncol; ++i){
    // Get the octave object
    const octave_value& obj_i = p.cell_value()(i);
    
    // Check if it is an SXMatrix
    if(isSXMatrix(obj_i)){
      CasADi::SXMatrix *m_i;
      bool ret = getSXMatrix(obj_i, m_i);
      if(!ret) return false;
      m[i] = *m_i;
    } else if(couldbeSXMatrix(obj_i)){
      bool ret = asSXMatrix(obj_i,m[i]);
      if(!ret) return false;
    }
  }
  return true;
}


%}
#endif // SWIGOCTAVE

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
    if (!result) SWIG_exception_fail(SWIG_TypeError,"SXMatrix cast failed");
  } else {  
    bool result=asVSXMatrix($input,v);
    if (!result) SWIG_exception_fail(SWIG_TypeError,"Expecting sequence(numpy.ndarray(SX/number) , SXMatrix, SX, number, sequence(SX/number))");
    $1 = &v;
  }
}

%typemap(typecheck,precedence=PRECEDENCE_SXMatrixVector) const std::vector< Matrix<SX> > & { $1 = couldbeSXMatrixVector($input); }
%typemap(freearg) const std::vector< Matrix<SX> > & {}

%fragment("generic_typemap", "header") {

template<class T>
class meta {
  public:
    /// Check if Python object is of type T
    static bool is(PyObject * p) {
      return istype(p,*meta<T>::name);
    };
    /// Convert Python object to pointer of type T
    static bool get_ptr(PyObject * p,T*& m) {
      void *pd = 0 ;
      int res = SWIG_ConvertPtr(p, &pd,*meta<T>::name, 0 );
      if (!SWIG_IsOK(res)) {
        return false;
      }
      m = reinterpret_cast< T*  >(pd);
      return true;
    };
    /// Convert Python object to type T
    static int as(PyObject * p,T&);
    /// Check if Python object could ultimately be converted to type T
    static bool couldbe(PyObject * p);
    static swig_type_info** name;
    static char expected_message[];
};

// $descriptor is quite useless: it only work in a typemap context.
template<> swig_type_info** meta< CasADi::Matrix<CasADi::SX> >::name = &SWIGTYPE_p_CasADi__MatrixT_CasADi__SX_t;

template<> char meta< CasADi::Matrix<CasADi::SX> >::expected_message[] = "Expecting one of: numpy.ndarray(SX/number) , SXMatrix, SX, number, sequence(SX/number)";


template <>
int meta< CasADi::Matrix<CasADi::SX> >::as(PyObject * p,CasADi::Matrix<CasADi::SX> &m) {return asSXMatrix(p,m);}

template <>
bool meta< CasADi::Matrix<CasADi::SX> >::couldbe(PyObject * p) {return couldbeSXMatrix(p);}

}


%define %my_generic_const_typemap(Type,Precedence) 
%typemap(in, fragment="generic_typemap") const Type & (Type m) {
  if (meta< Type >::is($input)) { // Type object get passed on as-is, and fast.
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

%typemap(typecheck,precedence=Precedence, fragment="generic_typemap") const Type & { $1 = meta< Type >::couldbe($input); }
%typemap(freearg) const Type  & {}

%enddef


%my_generic_const_typemap(CasADi::Matrix<CasADi::SX>,PRECEDENCE_SXMatrix);

} // namespace CasADi




%template(SXVector)             std::vector<CasADi::SX>;
%template(SXVectorVector)       std::vector<std::vector<CasADi::SX> > ;
%template(SXVectorVectorVector) std::vector< std::vector<std::vector<CasADi::SX> > > ;
%template(SXMatrix)             CasADi::Matrix<CasADi::SX>;


// NOTE: mx.py fails terminate called after throwing an instance of 'std::bad_alloc'
