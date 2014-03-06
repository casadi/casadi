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

#ifdef SWIGOCTAVE
%rename(__el_mul__) __mul__;
%rename(__el_div__) __div__;
%rename(__rel_mul__) __rmul__;
%rename(__rel_div__) __rdiv__;
%rename(__el_pow__) __pow__;
%rename(__rel_pow__) __rpow__;
%rename(__mul__) mul;
%rename(__rmul__) rmul;
%rename(__transpose__) trans;
%rename(__div__) __mrdivide__;
%rename(__rdiv__) __rmrdivide__;
%rename(__pow__) __mpower__;
%rename(__rpow__) __rmpower__;
#endif // SWIGOCTAVE

%{
#include "symbolic/matrix/sparsity.hpp"
#include "symbolic/matrix/matrix.hpp"
#include <sstream>
#include "symbolic/casadi_exception.hpp"

// to allow for typechecking
#include "symbolic/sx/sx_element.hpp"

// to typecheck for MX
#include "symbolic/mx/mx.hpp"
// to have prod available
#include "symbolic/mx/mx_tools.hpp"
%}

#ifdef SWIGPYTHON
%pythoncode %{
def prod(self,*args):
    raise Exception("'prod' is not supported anymore in CasADi. Use 'mul' to do matrix multiplication.")
def dot(self,*args):
    raise Exception("'dot' is not supported anymore in CasADi. Use 'mul' to do matrix multiplication.")
%}

%define %matrix_convertors
%pythoncode %{
        
    def toMatrix(self):
        import numpy as n
        return n.matrix(self.toArray())

    def __iter__(self):
      return self.data().__iter__()
        
%}
%enddef 
%define %matrix_helpers(Type)
%pythoncode %{
    @property
    def shape(self):
        return (self.size1(),self.size2())
        
    def reshape(self,arg):
        return _casadi_global.reshape(self,arg)
        
    @property
    def T(self):
        return _casadi_global.trans(self)
        
    def __getitem__(self,s):
        if isinstance(s,tuple) and len(s)==2:
          return self.__Cgetitem__(s[0],s[1])  
        return self.__Cgetitem__(s)

    def __setitem__(self,s,val):
        if isinstance(s,tuple) and len(s)==2:
          return self.__Csetitem__(s[0],s[1],val)  
        return self.__Csetitem__(s,val)
        
    def prod(self,*args):
        raise Exception("'prod' is not supported anymore in CasADi. Use 'mul' to do matrix multiplication.")
     
%}
%enddef 
    
%define %python_array_wrappers(arraypriority)
%pythoncode %{

  __array_priority__ = arraypriority

  def __array_wrap__(self,out_arr,context=None):
    if context is None:
      return out_arr
    name = context[0].__name__
    args = list(context[1])
    
    if len(context[1])==3:
      raise Exception("Error with %s. Looks like you are using an assignment operator, such as 'a+=b' where 'a' is a numpy type. This is not supported, and cannot be supported without changing numpy." % name)

    if "vectorized" in name:
        name = name[:-len(" (vectorized)")]
    
    conversion = {"multiply": "mul", "divide": "div", "true_divide": "div", "subtract":"sub","power":"pow","greater_equal":"ge","less_equal": "le", "less": "lt", "greater": "gt"}
    if name in conversion:
      name = conversion[name]
    if len(context[1])==2 and context[1][1] is self and not(context[1][0] is self):
      name = 'r' + name
      args.reverse()
    if not(hasattr(self,name)) or ('mul' in name):
      name = '__' + name + '__'
    fun=getattr(self, name)
    return fun(*args[1:])
     
     
  def __array__(self,*args,**kwargs):
    import numpy as n
    if len(args) > 1 and isinstance(args[1],tuple) and isinstance(args[1][0],n.ufunc) and isinstance(args[1][0],n.ufunc) and len(args[1])>1 and args[1][0].nin==len(args[1][1]):
      if len(args[1][1])==3:
        raise Exception("Error with %s. Looks like you are using an assignment operator, such as 'a+=b'. This is not supported when 'a' is a numpy type, and cannot be supported without changing numpy itself. Either upgrade a to a CasADi type first, or use 'a = a + b'. " % args[1][0].__name__)
      return n.array([n.nan])
    else:
      if hasattr(self,'__array_custom__'):
        return self.__array_custom__(*args,**kwargs)
      else:
        return self.toArray()
      
%}
%enddef
#endif // SWIGPYTHON

#ifdef SWIGOCTAVE
%define %matrix_helpers(Type)

  Type __hermitian__() const { return trans((*$self)); }
  
  std::vector<int> __dims__() const {
    std::vector<int> ret(2);
    ret[0] = $self->size1();
    ret[1] = $self->size2();
    return ret;
  }
  
%enddef
#endif //SWIGOCTAVE

#ifdef SWIGXML
%define %matrix_helpers(Type)
%enddef
#endif

#ifndef SWIGPYTHON
%define %matrix_convertors
%enddef
#endif

#ifdef SWIGOCTAVE
%rename(__paren__) indexed_one_based;
%rename(__paren__) indexed;
%rename(__paren_asgn__) indexed_one_based_assignment;
%rename(__paren_asgn__) indexed_assignment;
%rename(__vertcat__) vertcat;
%rename(__horzcat__) horzcat;
#endif
#ifdef SWIGPYTHON
%rename(__Cgetitem__) indexed_zero_based;
%rename(__Cgetitem__) indexed;
%rename(__Csetitem__) indexed_zero_based_assignment;
%rename(__Csetitem__) indexed_assignment;
#endif


#ifdef SWIGOCTAVE
namespace CasADi{
%extend Matrix<double> {
/// Create a 2D contiguous NP_DOUBLE numpy.ndarray

octave_value toSparse() {
  int nz = (*$self).size(), nr = (*$self).size1(), nc = (*$self).size2();
  
  Array<int> Ar(nz);
  Array<int> Ac(nz);
  
  std::vector<int> vc = (*$self).sparsity().getCol();
  Array<double> mydata(nz);
  const std::vector<double> &cdata = (*$self).data();
  
  for(int k=0; k<nz; ++k){
    Ac(k)=vc[k];
    Ar(k)=(*$self).sparsity().row()[k];
    mydata(k)=cdata[k];
  }
  
  return octave_value(SparseMatrix(mydata,Ar,Ac,nr,nc));
}

binopsrFull(CasADi::Matrix<double>)
binopsFull(const CasADi::Matrix<CasADi::SXElement> & b,,CasADi::Matrix<CasADi::SXElement>,CasADi::Matrix<CasADi::SXElement>)
binopsFull(const CasADi::MX & b,,CasADi::MX,CasADi::MX)

}; // extend Matrix<double>
} // namespace CasADi
#endif // SWIGOCTAVE

namespace CasADi{
  %extend Matrix<double> {

    void assign(const CasADi::Matrix<double>&rhs) { (*$self)=rhs; }
    %matrix_convertors
    %matrix_helpers(CasADi::Matrix<double>)

  }
  %extend Matrix<int> {

    void assign(const CasADi::Matrix<int>&rhs) { (*$self)=rhs; }
    %matrix_convertors
    %matrix_helpers(CasADi::Matrix<int>)

  }
}

#ifdef SWIGPYTHON
namespace CasADi{
%extend Matrix<double> {
/// Create a 2D contiguous NP_DOUBLE numpy.ndarray

#ifdef WITH_NUMPY
PyObject* arrayView() {
  if ($self->size()!=$self->numel()) 
    throw  CasADi::CasadiException("Matrix<double>::arrayview() can only construct arrayviews for dense DMatrices.");
  npy_intp dims[2];
  dims[0] = $self->size2();
  dims[1] = $self->size1();
  std::vector<double> &v = $self->data();
  PyArrayObject* temp = (PyArrayObject*) PyArray_New(&PyArray_Type, 2, dims, NPY_DOUBLE, NULL, &v[0], 0, NPY_ARRAY_CARRAY, NULL);
  PyObject* ret = PyArray_Transpose(temp,NULL);
  Py_DECREF(temp); 
  return ret;
}
#endif // WITH_NUMPY
    
%pythoncode %{
  def toArray(self,shared=False):
    import numpy as n
    if shared:
      if self.size()!=self.numel():
        raise Expection("toArray(shared=True) only possible for dense arrays.")
      return self.arrayView()
    else:
      r = n.zeros((self.size1(),self.size2()))
      self.get(r)
    return r
%}

%python_array_wrappers(999.0)

// The following code has some trickery to fool numpy ufunc.
// Normally, because of the presence of __array__, an ufunctor like nump.sqrt
// will unleash its activity on the output of __array__
// However, we wish DMatrix to remain a DMatrix
// So when we receive a call from a functor, we return a dummy empty array
// and return the real result during the postprocessing (__array_wrap__) of the functor.
%pythoncode %{
  def __array_custom__(self,*args,**kwargs):
    if "dtype" in kwargs and not(isinstance(kwargs["dtype"],n.double)):
      return n.array(self.toArray(),dtype=kwargs["dtype"])
    else:
      return self.toArray()
%}

%pythoncode %{
  def toCsc_matrix(self):
    import numpy as n
    from scipy.sparse import csc_matrix
    return csc_matrix( (list(self.data()),self.row(),self.colind()), shape = self.shape, dtype=n.double )

  def tocsc(self):
    return self.toCsc_matrix()

%}

%pythoncode %{
  def __float__(self):
    if self.numel()!=1:
      raise Exception("Only a scalar can be cast to a float")
    if self.size()==0:
      return 0.0
    return self.toScalar()
%}

%pythoncode %{
  def __int__(self):
    if self.numel()!=1:
      raise Exception("Only a scalar can be cast to an int")
    if self.size()==0:
      return 0
    return int(self.toScalar())
%}

%pythoncode %{
  def __nonzero__(self):
    if self.numel()!=1:
      raise Exception("Only a scalar can be cast to a float")
    if self.size()==0:
      return 0
    return self.toScalar()!=0
%}

%pythoncode %{
  def __abs__(self):
    return abs(self.__float__())
%}

binopsrFull(CasADi::Matrix<double>)
binopsFull(const CasADi::Matrix<CasADi::SXElement> & b,,CasADi::Matrix<CasADi::SXElement>,CasADi::Matrix<CasADi::SXElement>)
binopsFull(const CasADi::MX & b,,CasADi::MX,CasADi::MX)

}; // extend Matrix<double>

%extend Matrix<int> {

  %python_array_wrappers(998.0)

  %pythoncode %{
    def toArray(self):
      import numpy as n
      r = n.zeros((self.size1(),self.size2()))
      for i in range(r.shape[0]):
        for j in range(r.shape[1]):
          r[i,j] = self.elem(i,j)
      return r
  %}
  
  %pythoncode %{
    def __float__(self):
      if self.numel()!=1:
        raise Exception("Only a scalar can be cast to a float")
      if self.size()==0:
        return 0.0
      return float(self.toScalar())
  %}

  %pythoncode %{
    def __int__(self):
      if self.numel()!=1:
        raise Exception("Only a scalar can be cast to an int")
      if self.size()==0:
        return 0
      return self.toScalar()
  %}

  binopsrFull(CasADi::Matrix<int>)
  binopsFull(const CasADi::Matrix<CasADi::SXElement> & b,,CasADi::Matrix<CasADi::SXElement>,CasADi::Matrix<CasADi::SXElement>)
  binopsFull(const CasADi::Matrix<double> & b,,CasADi::Matrix<double>,CasADi::Matrix<double>)
  binopsFull(const CasADi::MX & b,,CasADi::MX,CasADi::MX)
  %pythoncode %{
    def __abs__(self):
      return int(self.__int__())
  %}
} // extend Matrix<int>


// Logic for pickling
%extend Sparsity {

  %pythoncode %{
    def __setstate__(self, state):
        if state:
          self.__init__(state["nrow"],state["ncol"],state["colind"],state["row"])
        else:
          self.__init__()

    def __getstate__(self):
        if self.isNull(): return {}
        return {"nrow": self.size1(), "ncol": self.size2(), "colind": numpy.array(self.colind(),dtype=int), "row": numpy.array(self.row(),dtype=int)}
  %}
  
}

%extend Matrix<int> {

  %pythoncode %{
    def __setstate__(self, state):
        sp = Sparsity.__new__(Sparsity)
        sp.__setstate__(state["sparsity"])
        self.__init__(sp,state["data"])

    def __getstate__(self):
        return {"sparsity" : self.sparsity().__getstate__(), "data": numpy.array(self.data(),dtype=int)}
  %}
  
}

%extend Matrix<double> {

  %pythoncode %{
    def __setstate__(self, state):
        sp = Sparsity.__new__(Sparsity)
        sp.__setstate__(state["sparsity"])
        self.__init__(sp,state["data"])

    def __getstate__(self):
        return {"sparsity" : self.sparsity().__getstate__(), "data": numpy.array(self.data(),dtype=float)}
  %}
  
}


} // namespace CasADi
#endif // SWIGPYTHON

  


namespace CasADi{

%{
#ifdef SWIGPYTHON
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
    m.start_ = (r->start == Py_None) ? std::numeric_limits<int>::min() : PyInt_AsLong(r->start);
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
#endif //SWIGPYTHON

#ifdef SWIGOCTAVE
/// CasADi::Slice
template<> char meta< CasADi::Slice >::expected_message[] = "Expecting Slice or number";

template <>
int meta< CasADi::Slice >::as(const octave_value& p,CasADi::Slice &m) {
  if (p.is_range()) {
    Range r = p.range_value();
    m.start_ = r.base()-1;
    m.stop_ = r.limit();
    m.step_ = r.inc();
  } else if (p.is_magic_colon()) {
    m.start_ = 0;
    m.stop_ = std::numeric_limits<int>::max();
  } else if (p.is_numeric_type()) {
    m.start_ = p.int_value()-1;
    m.stop_ = m.start_+1;
  } else {
    return false;
  }
  return true;
}

template <>
bool meta<  CasADi::Slice >::couldbe(const octave_value& p) {
  return p.is_range() || p.is_magic_colon()|| (p.is_real_scalar() && p.is_numeric_type());
}


/// CasADi::IndexList
template<> char meta< CasADi::IndexList >::expected_message[] = "Expecting Slice or number or list of ints";

template <>
int meta< CasADi::IndexList >::as(const octave_value& p,CasADi::IndexList &m) {
  if ((p.is_real_scalar() && p.is_numeric_type())) {
    m.type = CasADi::IndexList::INT;
    m.i = p.int_value()-1;
  } else if (meta< std::vector<int> >::couldbe(p)) {
    m.type = CasADi::IndexList::IVECTOR;
    bool result = meta< std::vector<int> >::as(p,m.iv);
    if (!result) return false;
    for (int k=0; k < m.iv.size();k++) m.iv[k]--;
  } else if (meta< CasADi::Slice>::couldbe(p)) {
    m.type = CasADi::IndexList::SLICE;
    return meta< CasADi::Slice >::as(p,m.slice);
  } else {
    return false;
  }
  return true;
}


template <>
bool meta<  CasADi::IndexList >::couldbe(const octave_value& p) {
  return meta< CasADi::Slice >::couldbe(p) || meta< std::vector<int> >::couldbe(p) || (p.is_real_scalar() && p.is_numeric_type());
}
#endif // SWIGOCTAVE
%}

%{
template<> swig_type_info** meta< CasADi::Slice >::name = &SWIGTYPE_p_CasADi__Slice;
template<> swig_type_info** meta< CasADi::IndexList >::name = &SWIGTYPE_p_CasADi__IndexList;
%}

%my_generic_const_typemap(PRECEDENCE_SLICE,CasADi::Slice);
%my_generic_const_typemap(PRECEDENCE_IndexVector,CasADi::IndexList);


#ifdef SWIGPYTHON

#ifdef WITH_NUMPY
/**

Accepts: 2D numpy.ndarray, numpy.matrix (contiguous, native byte order, datatype double)   - DENSE
         1D numpy.ndarray, numpy.matrix (contiguous, native byte order, datatype double)   - SPARSE
         2D scipy.csc_matrix
*/

%typemap(in,numinputs=1) (double * val,int len,int stride1, int stride2,SparsityType sp)  {
	PyObject* p = $input;
	$3 = 0;
	$4 = 0;
	if (is_array(p)) {
			if (!(array_is_native(p) && array_type(p)==NPY_DOUBLE))
			  SWIG_exception_fail(SWIG_TypeError, "Array should be native & of datatype double");
			  
	    if (!(array_is_contiguous(p))) {
	      if (PyArray_CHKFLAGS((PyArrayObject *) p,NPY_ALIGNED)) {
	        $3 = PyArray_STRIDE((PyArrayObject *) p,0)/sizeof(double);
	        $4 = PyArray_STRIDE((PyArrayObject *) p,1)/sizeof(double);
	      } else {
			   SWIG_exception_fail(SWIG_TypeError, "Array should be contiguous or aligned");
	      }
	    }
	    
			if (array_numdims(p)==2) {
				if (!(array_size(p,0)==arg1->size1() && array_size(p,1)==arg1->size2()) ) {
				  std::stringstream s;
				  s << "SWIG::typemap(in) (double *val,int len,SparsityType sp) " << std::endl;
				  s << "Array is not of correct shape.";
				  s << "Expecting shape (" << arg1->size1() << "," << arg1->size2() << ")" << ", but got shape (" << array_size(p,0) << "," << array_size(p,1) <<") instead.";
          const std::string tmp(s.str());
          const char* cstr = tmp.c_str();
			    SWIG_exception_fail(SWIG_TypeError,  cstr);
			  }
			  $5 = CasADi::DENSETRANS;
			  $2 = array_size(p,0)*array_size(p,1);
			  $1 = (double*) array_data(p);
			} else if (array_numdims(p)==1) {
				if (!(array_size(p,0)==arg1->size()) ) {
				  std::stringstream s;
				  s << "SWIG::typemap(in) (double *val,int len,SparsityType sp) " << std::endl;
				  s << "Array is not of correct size. Should match number of non-zero elements.";
				  s << "Expecting " << array_size(p,0) << " non-zeros, but got " << arg1->size() <<" instead.";
          const std::string tmp(s.str());
          const char* cstr = tmp.c_str();
			    SWIG_exception_fail(SWIG_TypeError,  cstr);
			  }
			  $5 = CasADi::SPARSE;
			  $2 = array_size(p,0);
			  $1 = (double*) array_data(p);
			} else {
			  SWIG_exception_fail(SWIG_TypeError, "Expecting 1D or 2D numpy.ndarray");
			}
	} else if (PyObjectHasClassName(p,"csc_matrix")) {
			$5 = CasADi::SPARSE;
			PyObject * narray=PyObject_GetAttrString( p, "data"); // narray needs to be decref'ed
			if (!(array_is_contiguous(narray) && array_is_native(narray) && array_type(narray)==NPY_DOUBLE))
			  SWIG_exception_fail(SWIG_TypeError, "csc_matrix should be contiguous, native & of datatype double");
			$2 = array_size(narray,0);
			if (!(array_size(narray,0)==arg1->size() ) ) {
					std::stringstream s;
				  s << "SWIG::typemap(in) (double *val,int len,SparsityType sp) " << std::endl;
				  s << "csc_matrix does not have correct number of non-zero elements.";
				  s << "Expecting " << arg1->size() << " non-zeros, but got " << array_size(narray,0) << " instead.";
          const std::string tmp(s.str());
          const char* cstr = tmp.c_str();
		      Py_DECREF(narray);
			    SWIG_exception_fail(SWIG_TypeError,  cstr);
			}
			$1 = (double*) array_data(narray);
			Py_DECREF(narray);
	} else {
			SWIG_exception_fail(SWIG_TypeError, "Unrecognised object");
	}
	
}

/**
Accepts: 2D numpy.ndarray, numpy.matrix (any setting of contiguous, native byte order, datatype)  - DENSE
         1D numpy.ndarray, numpy.matrix (any setting of contiguous, native byte order, datatype double) - SPARSE
         2D scipy.csc_matrix (any setting of contiguous, native byte order, datatype double) 
*/
%typemap(in,numinputs=1) (const double *val,int len,SparsityType sp) (PyArrayObject* array, int array_is_new_object=0)  {
	PyObject* p = $input;
	if (is_array(p)) {
			array = obj_to_array_contiguous_allow_conversion(p,NPY_DOUBLE,&array_is_new_object);
			if (array_numdims(array)==2) {
				if (!(array_size(array,0)==arg1->size1() && array_size(array,1)==arg1->size2()) ) {
				  std::stringstream s;
				  s << "SWIG::typemap(in) (const double *val,int len,SparsityType sp) " << std::endl;
				  s << "Array is not of correct shape.";
				  s << "Expecting shape (" << arg1->size1() << "," << arg1->size2() << ")" << ", but got shape (" << array_size(array,0) << "," << array_size(array,1) <<") instead.";
          const std::string tmp(s.str());
          const char* cstr = tmp.c_str();
			    SWIG_exception_fail(SWIG_TypeError,  cstr);
			  }
			  $3 = CasADi::DENSETRANS;
			  $2 = array_size(array,0)*array_size(array,1);
			  $1 = (double*) array_data(array);
			} else if (array_numdims(array)==1) {
				if (!(array_size(array,0)==arg1->size()) ) {
				  std::stringstream s;
				  s << "SWIG::typemap(in) (const double *val,int len,SparsityType sp) " << std::endl;
				  s << "Array is not of correct size. Should match number of non-zero elements.";
				  s << "Expecting " << arg1->size() << " non-zeros, but got " << array_size(array,0) << " instead.";
          const std::string tmp(s.str());
          const char* cstr = tmp.c_str();
			    SWIG_exception_fail(SWIG_TypeError,  cstr);
			  }
			  $3 = CasADi::SPARSE;
			  $2 = array_size(array,0);
			  $1 = (double*) array_data(array);
			} else {
			  SWIG_exception_fail(SWIG_TypeError, "Expecting 1D or 2D numpy.ndarray");
			}
	} else if (PyObjectHasClassName(p,"csc_matrix")) {
			$3 = CasADi::SPARSE;
			PyObject * narray=PyObject_GetAttrString( p, "data"); // narray needs to be decref'ed
			$2 = array_size(narray,0);
			if (!(array_size(narray,0)==arg1->size() ) ) {
					std::stringstream s;
				  s << "SWIG::typemap(in) (const double *val,int len,SparsityType sp) " << std::endl;
				  s << "csc_matrix does not have correct number of non-zero elements.";
				  s << "Expecting " << arg1->size() << " non-zeros, but got " << array_size(narray,0) << " instead.";
          const std::string tmp(s.str());
          const char* cstr = tmp.c_str();
          Py_DECREF(narray);
			    SWIG_exception_fail(SWIG_TypeError,  cstr);
			}
			array = obj_to_array_contiguous_allow_conversion(narray,NPY_DOUBLE,&array_is_new_object);
			$1 = (double*) array_data(array);
			Py_DECREF(narray);
	} else {
			SWIG_exception_fail(SWIG_TypeError, "Unrecognised object");
	}
	
}

%typemap(freearg) (const double *val,int len,SparsityType sp) {
    if (array_is_new_object$argnum && array$argnum) { Py_DECREF(array$argnum); }
}


%typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) (double * val,int len,int stride1, int stride2,SparsityType sp) {
  PyObject* p = $input;
  if (((is_array(p) && array_numdims(p) < 3)  && array_type(p)!=NPY_OBJECT)|| PyObjectHasClassName(p,"csc_matrix")) {
    $1=1;
  } else {
    $1=0;
  }
}

%typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) (const double * val,int len,SparsityType sp) {
  PyObject* p = $input;
  if (((is_array(p) && array_numdims(p) < 3)  && array_type(p)!=NPY_OBJECT)|| PyObjectHasClassName(p,"csc_matrix")) {
    $1=1;
  } else {
    $1=0;
  }
}
#endif // WITH_NUMPY
#endif // SWIGPYTHON


} // namespace CasADi
