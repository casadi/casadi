%{
#include "casadi/matrix/crs_sparsity.hpp"
#include "casadi/matrix/slice.hpp"
#include "casadi/matrix/matrix.hpp"
#include "casadi/matrix/matrix_tools.hpp"
#include "casadi/sx/sx.hpp"
#include "casadi/sx/sx_tools.hpp"


%}

// http://www.gnu.org/software/octave/doc/interpreter/Operator-Overloading.html#Operator-Overloading
// http://mentat.za.net/DaCodaAlFine.pdf
// "dispatch binary operator" can be found in octrun.swg: look for dispatch_binary_op in swig generated

%include "typemaps.i"
%include "casadi/matrix/crs_sparsity.hpp"
%include "casadi/matrix/slice.hpp"
%include "casadi/matrix/matrix.hpp"
%include "casadi/sx/sx.hpp"

%inline %{
template<> swig_type_info** meta< CasADi::SX >::name = &SWIGTYPE_p_CasADi__SX;
template<> swig_type_info** meta< CasADi::Matrix<CasADi::SX> >::name = &SWIGTYPE_p_CasADi__MatrixT_CasADi__SX_t;
template<> swig_type_info** meta< std::vector< CasADi::Matrix<CasADi::SX> > >::name = &SWIGTYPE_p_std__vectorT_CasADi__MatrixT_CasADi__SX_t_std__allocatorT_CasADi__MatrixT_CasADi__SX_t_t_t;
template<> swig_type_info** meta< std::vector< CasADi::SX > >::name = &SWIGTYPE_p_std__vectorT_CasADi__SX_std__allocatorT_CasADi__SX_t_t;
template<> swig_type_info** meta< std::vector< std::vector< CasADi::SX > > >::name = &SWIGTYPE_p_std__vectorT_std__vectorT_CasADi__SX_std__allocatorT_CasADi__SX_t_t_std__allocatorT_std__vectorT_CasADi__SX_std__allocatorT_CasADi__SX_t_t_t_t;

%}

#ifdef SWIGPYTHON
%extend CasADi::CRSSparsity{
    %pythoncode %{
        @property
        def shape(self):
            return (self.size1(),self.size2())
            
        def __array__(self,*args,**kwargs):
            return DMatrix(self,1).toArray()
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
  def constpow(x,y):
    if hasattr(x,'__constpow__'):
      return x.__constpow__(y)
  constpow=numpy.frompyfunc(constpow,2,1)
except:
  pass
%}
#endif // SWIGPYTHON

namespace CasADi {


%extend SX {
#ifdef SWIGPYTHON
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
  
  
  %pythoncode %{
  __array_priority__ = 1000.0
  
  def __array_wrap__(self,out_arr,context=None):
    if context is None:
      return out_arr
    name = context[0].__name__
    args = list(context[1])
    
    selfM = SXMatrix(self)
    if "vectorized" in name:
      name = name[:-len(" (vectorized)")]

    conversion = {"multiply": "mul", "divide": "div", "subtract":"sub","power":"pow"}
    if name in conversion:
      name = conversion[name]
    if len(context[1])==2 and context[1][1] is self:
      name = 'r' + name
      args.reverse()
    if not(hasattr(selfM,name)):
      name = '__' + name + '__'
    fun=getattr(selfM, name)
    return fun(*args[1:])

  def toArray(self):
    import numpy as n
    r = n.array((),dtype=object)
    r.resize(1,1)
    r[0,0] = self
    return r
      
  def __array__(self,*args,**kwargs):
    import numpy as n
    if len(args) > 1 and isinstance(args[1],tuple) and isinstance(args[1][0],n.ufunc):
      return n.array([1])
    else:
      return self.toArray()
  %}
  
  

  #endif // SWIGPYTHON
  
  #ifdef SWIGOCTAVE
  std::vector<int> __dims__() const {
    std::vector<int> ret(2);
    ret[0] = 1;
    ret[1] = 1;
    return ret;
  }
  
  #endif // SWIGOCTAVE
  
  binopsFull(const Matrix<double>& b,CasADi::Matrix<CasADi::SX>,CasADi::Matrix<CasADi::SX>,CasADi::Matrix<CasADi::SX>)
  binopsFull(const CasADi::Matrix<CasADi::SX> & b,,CasADi::Matrix<CasADi::SX>,CasADi::Matrix<CasADi::SX>)
    
};
  
%extend Matrix<SX>{
    // The constructor has to be added since SX::operator Matrix<SX does not work
    // Matrix<SX>(const SX&){ *$self
    #ifdef SWIGPYTHON
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
    #endif

    binopsFull(double b,CasADi::Matrix<CasADi::SX>,,CasADi::Matrix<CasADi::SX>)
    binopsFull(const CasADi::Matrix<double>& b,CasADi::Matrix<CasADi::SX>,,CasADi::Matrix<CasADi::SX>)
    #ifdef SWIGOCTAVE
    binopsFull(const CasADi::SX& b,CasADi::Matrix<CasADi::SX>,,CasADi::Matrix<CasADi::SX>)
    #endif // SWIGOCTAVE
  
    %python_matrix_convertors
    %python_matrix_helpers(CasADi::Matrix<CasADi::SX>)
       
    #ifdef SWIGPYTHON

    
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
  __array_priority__ = 1001.0
  %}
    
  %pythoncode %{
  def __array_wrap__(self,out_arr,context=None):
    name = context[0].__name__
    args = list(context[1])

    conversion = {"multiply": "mul", "divide": "div", "subtract":"sub","power":"pow"}
    if name in conversion:
      name = conversion[name]
    if len(context[1])==2 and context[1][1] is self:
      name = 'r' + name
      args.reverse()
    if not(hasattr(self,name)):
      name = '__' + name + '__'
    fun=getattr(self, name)
    return fun(*args[1:])
  %}

  %pythoncode %{
    def __array__(self,*args,**kwargs):
      import numpy as n
      if len(args) > 1 and isinstance(args[1],tuple) and isinstance(args[1][0],n.ufunc):
        return n.array([1])
      else:
        return self.toArray()
  %}
  #endif // SWIGPYTHON    
};



} // namespace CasADi

#ifdef SWIGPYTHON
#ifdef WITH_NUMPY
#include <numpy/arrayobject.h>
#endif // WITH_NUMPY

// Template instantiations
%template(vector_PyObject)    std::vector<PyObject*>;


#endif // SWIGPYTHON


%template(SXMatrixVector)       std::vector<CasADi::Matrix<CasADi::SX> > ;
%template(SXMatrixVectorVector) std::vector< std::vector<CasADi::Matrix<CasADi::SX> > > ;



%inline %{
template<> char meta< CasADi::SX >::expected_message[] = "Expecting SX or number";

// Explicit intialization of these two member functions, so we can use them in meta< CasADi::SX >
template<> int meta< CasADi::Matrix<CasADi::SX> >::as(GUESTOBJECT,CasADi::Matrix<CasADi::SX> &);
template<> bool meta< CasADi::Matrix<CasADi::SX> >::couldbe(GUESTOBJECT);

%}

/// CasADi::SX
#ifdef SWIGPYTHON
%inline %{


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

%}
#endif //SWIGPYTHON


/// CasADi::SX
#ifdef SWIGOCTAVE
%inline %{

template <>
int meta< CasADi::SX >::as(const octave_value& p,CasADi::SX &s) {
  NATIVERETURN(CasADi::SX, s)
  if ((p.is_real_scalar() && p.is_numeric_type())) {
    s=CasADi::SX(p.double_value());
    return true;
  } else if (meta< CasADi::Matrix< CasADi::SX > >::isa(p)) {
    CasADi::Matrix< CasADi::SX > m;
    meta< CasADi::Matrix< CasADi::SX > >::as(p,m);
    if (m.numel()==1 && m.size()==1) {
      s = m.at(0);
      return true;
    }
  }
  return false;
}

template <>
bool meta< CasADi::SX >::couldbe(const octave_value& p) {
  if (meta< CasADi::Matrix< CasADi::SX > >::isa(p)) {
    CasADi::Matrix< CasADi::SX > m;
    meta<CasADi::Matrix< CasADi::SX > >::as(p,m);
    if (m.numel()==1 && m.size()==1)
      return true;
  }
  return (meta< CasADi::SX >::isa(p) || (p.is_real_scalar() && p.is_numeric_type()) );
}

%}
#endif //SWIGOCTAVE

%inline %{
template<> char meta< CasADi::Matrix<CasADi::SX> >::expected_message[] = "Expecting one of: numpy.ndarray(SX/number) , SXMatrix, SX, number, sequence(SX/number)";
%}

/// CasADi::Matrix<CasADi::SX>
#ifdef SWIGPYTHON
%inline %{
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

%}
#endif //SWIGPYTHON

/// CasADi::Matrix<CasADi::SX>
#ifdef SWIGOCTAVE
%inline %{

template <>
int meta< CasADi::Matrix<CasADi::SX> >::as(const octave_value& p,CasADi::Matrix<CasADi::SX> &m) {
  NATIVERETURN(CasADi::Matrix<CasADi::SX>, m)
  NATIVERETURN(CasADi::Matrix<double>, m)
  NATIVERETURN(CasADi::SX, m)
  if((p.is_real_matrix() && p.is_numeric_type())){
    Matrix mat = p.matrix_value();
    m = CasADi::SXMatrix(mat.rows(),mat.cols(),0);
    for(int i=0; i<mat.rows(); ++i){
      for(int j=0; j<mat.cols(); ++j){
        m(i,j) = mat(i,j);
      }
    }
    return true;
  } 
  if ((p.is_real_scalar() && p.is_numeric_type())) {
    m = CasADi::SX(p.double_value());
    return true;
  }
  return false;
}

template <> bool meta< CasADi::Matrix<CasADi::SX> >::couldbe(const octave_value& p) { return meta< CasADi::Matrix<CasADi::SX> >::isa(p) || meta< CasADi::SX >::isa(p) || meta< CasADi::Matrix<double> >::isa(p)  || (p.is_real_matrix() && p.is_numeric_type()) || (p.is_real_scalar() && p.is_numeric_type());}
%}
#endif //SWIGOCTAVE

/// std::vector< CasADi::Matrix<CasADi::SX> >
#ifdef SWIGPYTHON
%inline %{
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

%}
#endif //SWIGPYTHON


%my_generic_const_typemap(PRECEDENCE_SX,CasADi::SX);
%my_generic_const_typemap(PRECEDENCE_SXMatrix,CasADi::Matrix<CasADi::SX>);
%my_generic_const_typemap(PRECEDENCE_SXMatrixVector,std::vector< CasADi::Matrix<CasADi::SX> >);

#ifdef SWIGOCTAVE
%meta_vector(CasADi::Matrix<CasADi::SX>);
#endif //SWIGOCTAVE
%meta_vector(std::vector<CasADi::SX>);
%meta_vector(CasADi::SX);

%template(SXVector)             std::vector<CasADi::SX>;
%template(SXVectorVector)       std::vector<std::vector<CasADi::SX> > ;
%template(SXVectorVectorVector) std::vector< std::vector<std::vector<CasADi::SX> > > ;
%template(SXMatrix)             CasADi::Matrix<CasADi::SX>;

